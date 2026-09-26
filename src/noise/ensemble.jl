# ─────────────────────────────────────────────────────────────────────────────
# NoiseEnsemble — the one container type every generator returns
# ─────────────────────────────────────────────────────────────────────────────

"""
    NoiseEnsemble(t, X)

`M` independent noise trajectories sampled on one uniform time grid.

- `t :: Vector{Float64}`: the time grid, of length `n_t` (uniformly spaced, increasing)
- `X :: Matrix{Float64}`: an `n_t × M` matrix where **column `j` is trajectory `j`**

Why a matrix: Julia arrays are column-major, so each trajectory is contiguous in
memory. Per-trajectory FFTs and autocovariances are then fast, and
ensemble statistics at a fixed time are just row operations (`X[k, :]`).

# Accessors
`times(e)`, `samples(e)`, `ntimes(e)`, `ntrajectories(e)`, `timestep(e)`.
`e[j]` returns trajectory `j` as a view, and `for x in e` iterates over trajectories.
`subensemble(e, idx)` selects trajectories without regenerating anything.
"""
struct NoiseEnsemble
    t::Vector{Float64}
    X::Matrix{Float64}

    function NoiseEnsemble(t::AbstractVector{<:Real}, X::AbstractMatrix{<:Real})
        n_t = length(t)
        n_t == size(X, 1) || throw(DimensionMismatch(
            "length(t) = $n_t but size(X, 1) = $(size(X, 1)); rows of X must be time points"))
        n_t >= 2 || throw(ArgumentError("a NoiseEnsemble needs at least 2 time points"))
        size(X, 2) >= 1 || throw(ArgumentError("a NoiseEnsemble needs at least 1 trajectory"))
        dt = (t[end] - t[1]) / (n_t - 1)
        dt > 0 || throw(ArgumentError("time grid must be increasing"))
        for k in 2:n_t
            isapprox(t[k] - t[k-1], dt; rtol=1e-6) ||
                throw(ArgumentError("time grid must be uniform (step $k differs from mean step $dt)"))
        end
        return new(Vector{Float64}(t), Matrix{Float64}(X))
    end
end

"Time grid of the ensemble."
times(e::NoiseEnsemble) = e.t

"The `n_t × M` sample matrix (column `j` = trajectory `j`)."
samples(e::NoiseEnsemble) = e.X

"Number of time points `n_t`."
ntimes(e::NoiseEnsemble) = size(e.X, 1)

"Number of trajectories `M`."
ntrajectories(e::NoiseEnsemble) = size(e.X, 2)

"Sampling interval Δt (computed from the whole grid, so it is robust to rounding)."
timestep(e::NoiseEnsemble) = (e.t[end] - e.t[1]) / (ntimes(e) - 1)

"""
    subensemble(e, idx) -> NoiseEnsemble

The trajectories `idx` (a range or vector of column indices) of `e`, as a new
ensemble on the same time grid. The data are copied; nothing is regenerated.
"""
subensemble(e::NoiseEnsemble, idx::AbstractVector{<:Integer}) = NoiseEnsemble(e.t, e.X[:, idx])

Base.length(e::NoiseEnsemble) = ntrajectories(e)
Base.getindex(e::NoiseEnsemble, j::Integer) = view(e.X, :, j)
Base.iterate(e::NoiseEnsemble, j::Int=1) = j > length(e) ? nothing : (e[j], j + 1)

function Base.show(io::IO, e::NoiseEnsemble)
    print(io, "NoiseEnsemble(", ntrajectories(e), " trajectories × ", ntimes(e),
          " time points, dt = ", timestep(e), ")")
end
