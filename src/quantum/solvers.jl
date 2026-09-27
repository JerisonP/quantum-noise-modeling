# ─────────────────────────────────────────────────────────────────────────────
# Brute force: solve the equation for every noise trajectory, then average
# ─────────────────────────────────────────────────────────────────────────────
#
# For ONE noise realisation ξ(t), the qubit obeys the von Neumann equation
#     dρ/dt = −i[H(t), ρ],     H(t) = f_x(t)/2 σx + ξ(t)/2 σz      (Eq. 3.12).
# As a 4×4 superoperator acting on vec_dm(ρ):
#     dS/dt = −i(1 ⊗ H − Hᵀ ⊗ 1) S,   S(0) = 1,    so ρ(t) = S(t)[ρ(0)].
# Averaging over realisations gives the exact noise-averaged evolution
#     ⟨ρ(t)⟩ = ⟨S(t)⟩[ρ(0)],      ⟨S(t)⟩ = L₀(t)⟨L_I(t)⟩          (Eq. 3.21).
# Every master equation for ⟨ρ⟩ (the SMNE of Eq. 2.57, the PLME of Eq. 1.37)
# approximates this. The brute-force average has no expansion in δ, so it has
# only two errors, and both are measured:
#   • statistical: finite number of trajectories, reported as standard errors;
#   • time step: checked with `timestep_convergence`.
#
# Between samples, ξ is linearly interpolated (the old code did the same).
#
# TWO INDEPENDENT SOLVERS. The tests require them to agree trajectory by trajectory.
#
# solver = :rk4    Solves the equation exactly as written above: classical RK4 on
#                  dS/dt in the frame of Eq. 3.12. It uses only `hamiltonian`, with
#                  no interaction picture and no closed forms. This is the old code's
#                  approach (superoperator ODE), with a fixed step so it can be
#                  convergence-tested.
#
# solver = :magnus (default) Interaction picture U = U_q,0 · U_I, with
#                  i dU_I/dt = V_I(t) U_I and V_I = ξ/2 · n(t)·σ. Each step uses the
#                  4th-order Magnus method at the two Gauss points
#                  t_a,b = t_mid ∓ h/(2√3) (Blanes, Casas, Oteo & Ros, Phys. Rep. 470,
#                  151 (2009)). For V_a = a·σ and V_b = b·σ, using
#                  [a·σ, b·σ] = 2i(a×b)·σ:
#                      U_step = exp(−i v·σ),   v = (h/2)(a + b) − (√3/6) h² (a × b).
#                  Every step is an exact SU(2) matrix, so the result stays unitary
#                  to rounding. It is faster than :rk4 and is the default.
#
# Both are 4th order: halving the step cuts the error 16× (tested).

_cross(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])

"""
    su2_exp(v) -> exp(−i v·σ)

Closed form for a real 3-vector v: cos|v|·1 − i sin|v|·(v̂·σ).
"""
function su2_exp(v)
    n = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
    iszero(n) && return copy(I2)
    s = sin(n) / n
    return cos(n) * I2 - im * s * (v[1] * σx + v[2] * σy + v[3] * σz)
end

"One Magnus-4 step of U_I on [t1, t2], with ξ varying linearly from ξ1 to ξ2."
function _magnus4_step(g::CosineGate, t1, t2, ξ1, ξ2)
    h = t2 - t1
    c = h / (2 * sqrt(3))
    tm = (t1 + t2) / 2
    ta, tb = tm - c, tm + c
    ξa = ξ1 + (ξ2 - ξ1) * (ta - t1) / h
    ξb = ξ1 + (ξ2 - ξ1) * (tb - t1) / h
    a = (ξa / 2) .* noise_axis(g, ta)            # V_I(t_a) = a·σ
    b = (ξb / 2) .* noise_axis(g, tb)
    v = (h / 2) .* (a .+ b) .- (sqrt(3) / 6 * h^2) .* _cross(a, b)
    return su2_exp(v)
end

"Liouvillian of the commutator: vec(−i[H, ρ]) = −i(1 ⊗ H − Hᵀ ⊗ 1) vec(ρ)."
_liouvillian(H::AbstractMatrix) = -im .* (kron(I2, H) .- kron(transpose(H), I2))

"One classical RK4 step of dS/dt = L(t) S on [t1, t2], with ξ varying linearly."
function _rk4_step(g::CosineGate, t1, t2, ξ1, ξ2, S)
    h = t2 - t1
    L1 = _liouvillian(hamiltonian(g, ξ1, t1))
    Lm = _liouvillian(hamiltonian(g, (ξ1 + ξ2) / 2, t1 + h / 2))
    L2 = _liouvillian(hamiltonian(g, ξ2, t2))
    k1 = L1 * S
    k2 = Lm * (S + (h / 2) * k1)
    k3 = Lm * (S + (h / 2) * k2)
    k4 = L2 * (S + h * k3)
    return S + (h / 6) * (k1 + 2k2 + 2k3 + k4)
end

"S(t_k) for one trajectory at the sorted grid indices `save_idx`."
function _propagate(g::CosineGate, t, ξ, save_idx, solver::Symbol, substeps::Integer)
    solver in (:magnus, :rk4) || throw(ArgumentError("solver must be :magnus or :rk4, got :$solver"))
    substeps >= 1 || throw(ArgumentError("substeps must be ≥ 1"))
    U = copy(I2)                  # :magnus state (interaction picture)
    S = kron(I2, I2)              # :rk4 state (frame of Eq. 3.12)
    out = Vector{Matrix{ComplexF64}}(undef, length(save_idx))
    j = 1
    for k in 1:length(t)
        if k > 1
            for s in 1:substeps
                f1, f2 = (s - 1) / substeps, s / substeps
                t1 = t[k-1] + f1 * (t[k] - t[k-1])
                t2 = t[k-1] + f2 * (t[k] - t[k-1])
                ξ1 = ξ[k-1] + f1 * (ξ[k] - ξ[k-1])
                ξ2 = ξ[k-1] + f2 * (ξ[k] - ξ[k-1])
                if solver === :magnus
                    U = _magnus4_step(g, t1, t2, ξ1, ξ2) * U
                else
                    S = _rk4_step(g, t1, t2, ξ1, ξ2, S)
                end
            end
        end
        if j <= length(save_idx) && save_idx[j] == k
            out[j] = solver === :magnus ? superoperator(ideal_propagator(g, t[k]) * U) : copy(S)
            j += 1
        end
    end
    return out
end

"""
    propagate_trajectory(g, t, ξ; solver=:magnus, substeps=1) -> Vector of S(t_k)

The 4×4 superoperator S(t_k) at every grid time, for ONE noise trajectory sampled
as `ξ[k] = ξ(t[k])`. It is expressed in the frame of Eq. 3.12, so ρ(t_k) = S(t_k)[ρ(0)].
`substeps` splits every grid interval (ξ is still linearly interpolated). It is
used for convergence checks and to make the two solvers agree to high precision.
"""
function propagate_trajectory(g::CosineGate, t::AbstractVector{<:Real}, ξ::AbstractVector{<:Real};
                              solver::Symbol=:magnus, substeps::Integer=1)
    length(t) == length(ξ) || throw(DimensionMismatch("t and ξ differ in length"))
    return _propagate(g, t, ξ, 1:length(t), solver, substeps)
end

"""
    GateResult

Output of [`simulate_gate`](@ref):

- `gate`: the `CosineGate`
- `times`: the saved times
- `S`: the noise-averaged superoperator ⟨S(t)⟩ at each saved time (Eq. 3.21), so
  `apply_channel(S[i], ρ0)` is ⟨ρ(times[i])⟩
- `S_final`: every trajectory's own S(t_g). All standard errors come from these.
- `ensemble`: the exact noise that was used, so it can also be validated
- `solver`: `:magnus` or `:rk4`
"""
struct GateResult
    gate::CosineGate
    times::Vector{Float64}
    S::Vector{Matrix{ComplexF64}}
    S_final::Vector{Matrix{ComplexF64}}
    ensemble::NoiseEnsemble
    solver::Symbol
end

"""
    simulate_gate(g, ens; solver=:magnus, save_every=1, substeps=1) -> GateResult
    simulate_gate(g, model, n_samples; dt=g.tg/1000, rng, solver, save_every, substeps)

Brute-force solution: solve Eq. 3.12 for every trajectory of the noise ensemble
`ens`, whose time grid must span [0, tg], and average. The second form generates
the ensemble from `model` first. Keeping generation separate means the *same*
ensemble can be passed to `validate_noise_model`.
"""
function simulate_gate(g::CosineGate, ens::NoiseEnsemble; solver::Symbol=:magnus,
                       save_every::Integer=1, substeps::Integer=1)
    t = times(ens)
    isapprox(t[1], 0; atol=1e-12 * g.tg) && isapprox(t[end], g.tg; rtol=1e-9) ||
        throw(ArgumentError("noise time grid must span [0, tg] = [0, $(g.tg)], got [$(t[1]), $(t[end])]"))
    save_every >= 1 || throw(ArgumentError("save_every must be ≥ 1"))
    idx = collect(1:save_every:length(t))
    idx[end] == length(t) || push!(idx, length(t))
    Ssum = [zeros(ComplexF64, 4, 4) for _ in idx]
    finals = Vector{Matrix{ComplexF64}}(undef, ntrajectories(ens))
    for (j, ξ) in enumerate(ens)
        Ss = _propagate(g, t, ξ, idx, solver, substeps)
        for i in eachindex(idx)
            Ssum[i] .+= Ss[i]
        end
        finals[j] = Ss[end]
    end
    M = ntrajectories(ens)
    return GateResult(g, t[idx], [A ./ M for A in Ssum], finals, ens, solver)
end

function simulate_gate(g::CosineGate, model::AbstractNoiseModel, n_samples::Integer;
                       dt::Real=g.tg / 1000, rng::AbstractRNG=default_rng(), kwargs...)
    return simulate_gate(g, generate_ensemble(model, (0.0, g.tg), dt, n_samples; rng); kwargs...)
end

"""
    timestep_convergence(g, ens; factors=(1, 2, 4, 8), solver=:magnus) -> Vector of rows

Is the time step small enough? The gate is re-run on the SAME trajectories, keeping
every f-th sample. For an exact generator such as OU, every f-th sample of an exact
path is itself an exact path on the coarser grid. Differences between rows are
therefore time-step error only, not Monte-Carlo noise, and they are paired
trajectory by trajectory.

Each row contains `factor`, `dt`, `ε` and `ε_se` (Eq. 3.25), and
`Δε = ε(f·dt) − ε(finest)` with its paired standard error `Δε_se`.
The time step is converged when |Δε| at f = 2 is far below `ε_se`.

Note: the time-step error comes mostly from representing a rough noise path by
straight lines, not from the integrator. `dt` must resolve τ_c as well as t_g.
"""
function timestep_convergence(g::CosineGate, ens::NoiseEnsemble; factors=(1, 2, 4, 8),
                              solver::Symbol=:magnus, substeps::Integer=1)
    fs = sort(collect(Int, factors))
    n = ntimes(ens) - 1
    all(f -> f >= 1 && n % f == 0, fs) ||
        throw(ArgumentError("every factor must divide the number of steps ($n)"))
    t, X = times(ens), samples(ens)
    base = Float64[]
    rows = map(fs) do f
        sub = NoiseEnsemble(t[1:f:end], X[1:f:end, :])
        e = trajectory_errors(simulate_gate(g, sub; solver, save_every=max(1, ntimes(sub) - 1), substeps))
        isempty(base) && append!(base, e)
        d = e .- base
        M = length(e)
        (factor=f, dt=f * timestep(ens), ε=mean(e), ε_se=std(e) / sqrt(M),
         Δε=mean(d), Δε_se=std(d) / sqrt(M))
    end
    return rows
end