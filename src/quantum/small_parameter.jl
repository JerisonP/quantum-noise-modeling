# ─────────────────────────────────────────────────────────────────────────────
# The teammates' small parameter δ, and brute-force reference curves against it
# ─────────────────────────────────────────────────────────────────────────────

"""
    small_noise_parameter(σ, τ_c, θ, t_g) -> δ                 (thesis Eq. 3.20)
    small_noise_parameter(g::CosineGate, m::OUNoiseModel)

    δ = σ² τ_c / √(1 + θ² τ_c² / t_g²)

The small parameter chosen to characterise when the perturbative master equation
is valid (thesis §3.2.2). The thesis figures plot everything against δ at fixed
t_g/τ_c, varying δ through σ²τ_c. Here σ and τ_c are those of the OU process
(Eq. 3.10), i.e. `OUNoiseModel(σ, τ_c)`.

Units: ξ is an angular frequency, so σ²τ_c has units of 1/time. δ becomes a pure
number only once a time unit is fixed, so use the same time unit as whoever you
are comparing against.
"""
small_noise_parameter(σ::Real, τc::Real, θ::Real, tg::Real) = σ^2 * τc / sqrt(1 + θ^2 * τc^2 / tg^2)
small_noise_parameter(g::CosineGate, m::OUNoiseModel) = small_noise_parameter(m.σ, m.τ_c, g.θ, g.tg)

"The OU σ that gives small parameter `δ`, i.e. Eq. 3.20 solved for σ."
sigma_for_delta(δ::Real, τc::Real, θ::Real, tg::Real) = sqrt(δ * sqrt(1 + θ^2 * τc^2 / tg^2) / τc)

"""
    delta_sweep(g, τ_c, δs, n_samples; dt, rng, solver=:magnus, nbatches=20, substeps=1)

Brute-force reference values for the thesis figures, at each δ in `δs`:

- `δ`, and `σ` from Eq. 3.20
- `ε`, `ε_se`: ⟨ε⟩_ξ (Eq. 3.25), as in Figs. 3.20–3.22
- `states`: for each cardinal state (zp, zm, xp, xm, yp, ym), `ρ`, `se` (⟨ρ_j(t_g)⟩ and its
  element-wise SE, Eq. 3.21), and `λ`, `λ_se` (its eigenvalues, as in Figs. 3.2–3.19)

The noise is OU (Eq. 3.10), zero mean, and stationary from t = 0. One unit-σ ensemble
is generated and multiplied by σ for every δ. This is exact, because the OU equation
is linear in σ (ξ_σ = σ·ξ₁). Using the same random numbers makes the curve smooth in
δ, and it makes differences between δ values far more precise than independent runs
would.

Default time step: `dt = min(t_g, τ_c)/1000`, so that both the gate and the noise
correlation time are resolved. Check it with [`timestep_convergence`](@ref).
"""
function delta_sweep(g::CosineGate, τc::Real, δs, n_samples::Integer;
                     dt::Real=_default_dt(g, τc), rng::AbstractRNG=default_rng(),
                     solver::Symbol=:magnus, nbatches::Integer=20, substeps::Integer=1)
    unit = generate_ensemble(OUNoiseModel(1.0, τc), (0.0, g.tg), dt, n_samples; rng)
    return map(collect(δs)) do δ
        σ = sigma_for_delta(δ, τc, g.θ, g.tg)
        ens = NoiseEnsemble(times(unit), σ .* samples(unit))
        r = simulate_gate(g, ens; solver, save_every=ntimes(ens) - 1, substeps)
        err = average_error(r)
        states = map(cardinal_states()) do ρ0
            fs = final_state(r, ρ0)
            ev = final_state_eigenvalues(r, ρ0; nbatches)
            (ρ=fs.ρ, se=fs.se, λ=ev.λ, λ_se=ev.se)
        end
        (δ=δ, σ=σ, ε=err.ε, ε_se=err.se, states=states)
    end
end

"`min(t_g, τ_c)/1000`, rounded so that t_g is an exact whole number of steps."
_default_dt(g::CosineGate, τc::Real) = g.tg / ceil(Int, 1000 * g.tg / min(g.tg, τc) - 1e-9)