# ─────────────────────────────────────────────────────────────────────────────
# CSV input/output: self-describing, exact, and read by column NAME
# ─────────────────────────────────────────────────────────────────────────────
#
# Every file written here has the same layout:
#
#     # QuantumNoiseSimulator CSV v1          ← format line (checked on load)
#     # kind = delta_sweep                     ← what the file holds (checked on load)
#     # tg = 1.0                               ← any other metadata, one "key = value" per line
#     delta,sigma,eps,eps_se,...               ← ONE header row of column names
#     0.02,0.1929,0.00331,0.000102,...         ← data, one row per record
#
# Two rules make the round trip safe:
#   1. Numbers are written with Julia's shortest exact representation
#      (`string(x::Float64)`), so reading a file back gives bit-for-bit the
#      same Float64. The tests check this with `==`, not `≈`.
#   2. Loaders find columns by NAME, never by position. The old `load_propagator`
#      assumed columns (param, time, …), but `save_propagator` wrote
#      (time, …), so every value was read from the wrong column.

const _CSV_FORMAT = "QuantumNoiseSimulator CSV v1"

_fmt(x::Real) = string(Float64(x))

"Write metadata, a header row and a numeric matrix (rows = records)."
function _write_table(path::AbstractString, kind::AbstractString, meta, cols::Vector{String},
                      data::AbstractMatrix{<:Real})
    size(data, 2) == length(cols) ||
        throw(DimensionMismatch("$(length(cols)) column names for $(size(data, 2)) columns"))
    all(c -> !occursin(',', c), cols) || throw(ArgumentError("column names must not contain commas"))
    dir = dirname(path)
    isempty(dir) || mkpath(dir)
    open(path, "w") do io
        println(io, "# ", _CSV_FORMAT)
        println(io, "# kind = ", kind)
        for (k, v) in meta
            s = string(v)
            occursin('\n', s) && throw(ArgumentError("metadata value for $k contains a newline"))
            println(io, "# ", k, " = ", s)
        end
        println(io, join(cols, ','))
        for i in axes(data, 1)
            println(io, join((_fmt(x) for x in view(data, i, :)), ','))
        end
    end
    return path
end

"""
    read_table(path) -> (meta, cols, data)

Read any file written by this package: `meta::Dict{String,String}` (including
`"kind"`), the column names, and the numeric data as a `Matrix{Float64}` whose
columns are in file order. Use it to inspect a file whose kind you don't know yet.
"""
function read_table(path::AbstractString)
    lines = [rstrip(l) for l in readlines(path)]         # rstrip also drops Windows "\r"
    isempty(lines) && throw(ArgumentError("$path is empty"))
    strip(lines[1]) == "# " * _CSV_FORMAT ||
        throw(ArgumentError("$path is not a $_CSV_FORMAT file (first line: $(repr(lines[1])))"))
    meta = Dict{String,String}()
    i = 2
    while i <= length(lines) && startswith(lines[i], "#")
        body = strip(lines[i][2:end])
        k = findfirst(" = ", body)
        k === nothing || (meta[strip(body[1:first(k)-1])] = strip(body[last(k)+1:end]))
        i += 1
    end
    i <= length(lines) || throw(ArgumentError("$path has no header row"))
    cols = [String(strip(c)) for c in split(lines[i], ',')]
    rows = [parse.(Float64, split(l, ',')) for l in lines[i+1:end] if !isempty(strip(l))]
    all(r -> length(r) == length(cols), rows) ||
        throw(ArgumentError("$path: a data row does not have $(length(cols)) values"))
    data = isempty(rows) ? zeros(0, length(cols)) : permutedims(reduce(hcat, rows))
    return (meta=meta, cols=cols, data=data)
end

function _read_kind(path, kind)
    t = read_table(path)
    get(t.meta, "kind", "") == kind ||
        throw(ArgumentError("$path holds kind = $(get(t.meta, "kind", "?")), expected $kind"))
    return t
end

"Column `name` of a table, or `nothing` if the file does not have it."
function _col(t, name)
    j = findfirst(==(name), t.cols)
    return j === nothing ? nothing : t.data[:, j]
end

function _col!(t, name)
    c = _col(t, name)
    c === nothing && throw(ArgumentError("missing required column \"$name\""))
    return c
end

_colornan(t, name) = (c = _col(t, name); c === nothing ? fill(NaN, size(t.data, 1)) : c)

# ── Noise ensembles ──────────────────────────────────────────────────────────

"""
    save_ensemble(path, ens; meta=()) -> path

Columns `t, xi_1, …, xi_M`, one row per time point. This is the exact noise
that drove a simulation, so someone else can reuse it (or validate it) later.
Pass extra metadata as pairs, e.g. `meta = ("model" => string(model), "seed" => 1)`.
"""
function save_ensemble(path::AbstractString, ens::NoiseEnsemble; meta=())
    M = ntrajectories(ens)
    cols = ["t"; ["xi_$j" for j in 1:M]]
    info = ("n_trajectories" => M, "dt" => timestep(ens), meta...)
    return _write_table(path, "ensemble", info, cols, hcat(times(ens), samples(ens)))
end

"`load_ensemble(path) -> NoiseEnsemble`, bit-for-bit what was saved."
function load_ensemble(path::AbstractString)
    t = _read_kind(path, "ensemble")
    tt = _col!(t, "t")
    M = count(c -> startswith(c, "xi_"), t.cols)
    X = Matrix{Float64}(undef, length(tt), M)
    for j in 1:M
        X[:, j] = _col!(t, "xi_$j")
    end
    return NoiseEnsemble(tt, X)
end

# ── Averaged channels ⟨S(t)⟩ ─────────────────────────────────────────────────

const _S_COLS = [s for i in 1:4 for j in 1:4 for s in ("S$(i)$(j)_re", "S$(i)$(j)_im")]

"""
    save_channel(path, r::GateResult; all_times=false) -> path
    save_channel(path, times, S; meta=()) -> path

The averaged superoperator ⟨S(t)⟩ (Eq. 3.21), which maps vec_dm(ρ(0)) to
vec_dm(⟨ρ(t)⟩). Columns are `time`, then `Sij_re, Sij_im` for row i and column j.
`all_times=false` writes only t = t_g. The second form saves any list of 4×4 maps,
for example the output of `tcl2_evolution`.
"""
function save_channel(path::AbstractString, times_::AbstractVector{<:Real},
                      S::AbstractVector{<:AbstractMatrix}; meta=())
    length(times_) == length(S) || throw(DimensionMismatch("times and S differ in length"))
    all(A -> size(A) == (4, 4), S) || throw(DimensionMismatch("every S must be 4×4"))
    data = zeros(length(S), 33)
    for (n, A) in enumerate(S)
        data[n, 1] = times_[n]
        k = 2
        for i in 1:4, j in 1:4                      # same order as _S_COLS
            data[n, k], data[n, k+1] = real(A[i, j]), imag(A[i, j])
            k += 2
        end
    end
    return _write_table(path, "channel", meta, ["time"; _S_COLS], data)
end

function save_channel(path::AbstractString, r::GateResult; all_times::Bool=false)
    meta = ("theta" => r.gate.θ, "tg" => r.gate.tg, "solver" => r.solver,
            "n_trajectories" => ntrajectories(r.ensemble), "dt" => timestep(r.ensemble))
    idx = all_times ? eachindex(r.times) : [lastindex(r.times)]
    return save_channel(path, r.times[idx], r.S[idx]; meta)
end

"""
    load_channel(path) -> (times, S, meta)

Read a file written by [`save_channel`](@ref). `S[k]` is a 4×4 `Matrix{ComplexF64}`.
"""
function load_channel(path::AbstractString)
    t = _read_kind(path, "channel")
    cols = [complex.(_col!(t, "S$(i)$(j)_re"), _col!(t, "S$(i)$(j)_im")) for i in 1:4, j in 1:4]
    S = [[cols[i, j][k] for i in 1:4, j in 1:4] for k in 1:size(t.data, 1)]
    return (times=_col!(t, "time"), S=S, meta=t.meta)
end

# ── δ sweeps: the reference curves, and the hand-in format for comparisons ───

const _STATE_NAMES = (:zp, :zm, :xp, :xm, :yp, :ym)

function _sweep_cols()
    cols = ["delta", "sigma", "eps", "eps_se"]
    for s in _STATE_NAMES
        for a in 1:2, b in 1:2, part in ("re", "im")
            push!(cols, "$(s)_rho$(a)$(b)_$(part)")
        end
        for a in 1:2, b in 1:2, part in ("re", "im")
            push!(cols, "$(s)_rho$(a)$(b)_se_$(part)")
        end
        append!(cols, ["$(s)_lambda1_re", "$(s)_lambda1_im", "$(s)_lambda2_re", "$(s)_lambda2_im",
                       "$(s)_lambda1_se", "$(s)_lambda2_se"])
    end
    return cols
end

_field(x, name, default) = hasproperty(x, name) ? getproperty(x, name) : default

function _sweep_row(row)
    v = Float64[row.δ, _field(row, :σ, NaN), row.ε, _field(row, :ε_se, NaN)]
    nanρ = fill(complex(NaN, NaN), 2, 2)
    states = _field(row, :states, nothing)
    for s in _STATE_NAMES
        st = states === nothing ? nothing : _field(states, s, nothing)
        ρ = st === nothing ? nanρ : _field(st, :ρ, nanρ)
        se = st === nothing ? nanρ : _field(st, :se, nanρ)
        λ = st === nothing ? fill(complex(NaN, NaN), 2) : _field(st, :λ, fill(complex(NaN, NaN), 2))
        λse = st === nothing ? [NaN, NaN] : _field(st, :λ_se, [NaN, NaN])
        for a in 1:2, b in 1:2
            append!(v, (real(ρ[a, b]), imag(ρ[a, b])))
        end
        for a in 1:2, b in 1:2
            append!(v, (real(se[a, b]), imag(se[a, b])))
        end
        append!(v, (real(λ[1]), imag(λ[1]), real(λ[2]), imag(λ[2]), λse[1], λse[2]))
    end
    return v
end

"""
    save_delta_sweep(path, rows; meta=()) -> path

Write the output of [`delta_sweep`](@ref) (or [`tcl2_delta_sweep`](@ref)) as one row
per δ. The columns are:

- `delta, sigma, eps, eps_se`
- for each state s ∈ zp, zm, xp, xm, yp, ym:
  - `s_rhoab_re`, `s_rhoab_im`: ⟨ρ_s(t_g)⟩ elements
  - `s_rhoab_se_re`, `s_rhoab_se_im`: their standard errors
  - `s_lambda1_re`, `s_lambda1_im`, `s_lambda2_re`, `s_lambda2_im`: the eigenvalues
  - `s_lambda1_se`, `s_lambda2_se`: SE of the eigenvalues' real parts

Anything a row does not have (for example SEs for a master-equation result) is written as NaN.
"""
function save_delta_sweep(path::AbstractString, rows; meta=())
    cols = _sweep_cols()
    data = permutedims(reduce(hcat, [_sweep_row(r) for r in rows]))
    return _write_table(path, "delta_sweep", meta, cols, data)
end

"""
    load_delta_sweep(path) -> Vector of rows

Read a δ-sweep file back into the same structure `delta_sweep` returns:
`(δ, σ, ε, ε_se, states = (zp = (ρ, se, λ, λ_se), …))`.

This is also the HAND-IN FORMAT for a master-equation result. Only the columns
`delta` and `eps` are required. Every other column is optional, and a missing column
loads as NaN. So a teammate can send `delta,eps,xp_lambda1_re,xp_lambda2_re,…` with the
format line `# QuantumNoiseSimulator CSV v1` and the line `# kind = delta_sweep` on top,
and it loads and plots next to the brute-force curve.
"""
function load_delta_sweep(path::AbstractString)
    t = _read_kind(path, "delta_sweep")
    δ, ε = _col!(t, "delta"), _col!(t, "eps")
    σ, εse = _colornan(t, "sigma"), _colornan(t, "eps_se")
    cplx(name, k) = complex(_colornan(t, name * "_re")[k], _colornan(t, name * "_im")[k])
    return map(eachindex(δ)) do k
        states = NamedTuple{_STATE_NAMES}(map(_STATE_NAMES) do s
            ρ = [cplx("$(s)_rho$(a)$(b)", k) for a in 1:2, b in 1:2]
            se = [cplx("$(s)_rho$(a)$(b)_se", k) for a in 1:2, b in 1:2]
            λ = [cplx("$(s)_lambda1", k), cplx("$(s)_lambda2", k)]
            λse = [_colornan(t, "$(s)_lambda1_se")[k], _colornan(t, "$(s)_lambda2_se")[k]]
            (ρ=ρ, se=se, λ=λ, λ_se=λse)
        end)
        (δ=δ[k], σ=σ[k], ε=ε[k], ε_se=εse[k], states=states)
    end
end