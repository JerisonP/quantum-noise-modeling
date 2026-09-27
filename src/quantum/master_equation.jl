# ─────────────────────────────────────────────────────────────────────────────
# The 2nd-order master equation: an approximate equation, validated against brute force
# ─────────────────────────────────────────────────────────────────────────────
#
# THIS IS NOT THE REFERENCE. The brute-force average (`simulate_gate`) is the
# reference. This file solves the kind of equation that brute force is used to check.
#
# At 2nd order, the thesis's SMNE (Eq. 2.57 with Ξ⁽⁴⁾ = 0) and the 2nd-order PLME
# are identical (thesis §3.2.3). For Gaussian noise, both are the time-local
# 2nd-order equation
#
#     d⟨ρ_I⟩/dt = −∫₀ᵗ C(t − s) [V(t), [V(s), ⟨ρ_I⟩]] ds,     V(t) = ½ n(t)·σ,
#
# where C(τ) = ⟨ξ(t)ξ(t+τ)⟩ (compare PLME Eqs. 1.34–1.36 with g = 1, Â(t) = V(t)).
# Why it is in the package:
#   1. It reproduces the teammates' 2nd-order curves from the same gate and noise
#      objects. That confirms the conventions (σ, time unit, δ) before comparing.
#   2. It gives the brute-force pipeline a test it must pass. At θ = 0 every V(t)
#      commutes, so this equation is EXACT for Gaussian noise.
#   3. Comparing it with brute force shows where 2nd order breaks down as δ grows.
#
# Numerics: V(t) is evaluated on a grid of spacing h/2, the memory integral
# ∫₀ᵗ C(t−s)V(s)ds uses the trapezoid rule on that grid, and the equation is stepped
# with classical RK4 of step h (the midpoints use the half-grid). The error is O(h²),
# limited by the trapezoid rule.

"""
    tcl2_evolution(g, model; nsteps=2000) -> (times, S)

Solve the 2nd-order time-local master equation above for the averaged superoperator.
It returns `S[k]`, the map from ρ(0) to ⟨ρ(times[k])⟩ in the frame of Eq. 3.12 (the
same convention as `GateResult.S`). Use it with `apply_channel`, `fidelity_map`,
`is_cptp` and `density_matrix_eigenvalues`.

The model must have a continuous-time autocovariance: `OUNoiseModel` or
`BandLimitedOneOverFNoiseModel`.
"""
function tcl2_evolution(g::CosineGate, model::Union{OUNoiseModel,BandLimitedOneOverFNoiseModel};
                        nsteps::Integer=2000)
    nsteps >= 1 || throw(ArgumentError("nsteps must be ≥ 1"))
    tt = range(0, g.tg; length=2nsteps + 1)          # half-step grid
    hh = step(tt)
    ad = map(tt) do x                                 # superoperator of ρ → [V(x), ρ]
        n = noise_axis(g, x)
        V = (n[2] * σy + n[3] * σz) / 2
        kron(I2, V) .- kron(transpose(V), I2)
    end
    C = [theoretical_autocovariance(model, l * hh) for l in 0:length(tt)-1]
    K = [zeros(ComplexF64, 4, 4) for _ in tt]         # generator K(t) = −ad(t)·∫₀ᵗ C(t−s) ad(s) ds
    for k in 2:length(tt)
        B = zeros(ComplexF64, 4, 4)
        for j in 1:k
            w = (j == 1 || j == k) ? hh / 2 : hh
            B .+= (w * C[k-j+1]) .* ad[j]
        end
        K[k] = -ad[k] * B
    end
    h = 2hh
    Λ = kron(I2, I2)                                  # interaction-picture map ⟨L_I(t)⟩
    ts = collect(tt[1:2:end])
    S = [superoperator(ideal_propagator(g, 0.0)) * Λ]
    for i in 1:nsteps
        K1, Km, K2 = K[2i-1], K[2i], K[2i+1]
        k1 = K1 * Λ
        k2 = Km * (Λ + (h / 2) * k1)
        k3 = Km * (Λ + (h / 2) * k2)
        k4 = K2 * (Λ + h * k3)
        Λ = Λ + (h / 6) * (k1 + 2k2 + 2k3 + k4)
        push!(S, superoperator(ideal_propagator(g, ts[i+1])) * Λ)
    end
    return (times=ts, S=S)
end