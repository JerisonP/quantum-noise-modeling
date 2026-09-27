# ─────────────────────────────────────────────────────────────────────────────
# The physical model: thesis Chapter 3, Eqs. 3.1–3.15
# ─────────────────────────────────────────────────────────────────────────────
#
# Lab frame (Eq. 3.9):
#     H_q(t) = ½[ω₀ + ξ(t)] σz + f_x(t) cos(ω_d t) σx
# The noise ξ(t) shifts the qubit splitting ω₀ (detuning noise), so it sits next
# to ω₀ on σz/2.
#
# Frame rotating at ω_d (Eqs. 3.2–3.3), rotating-wave approximation, resonant
# drive ω_d = ω_q (Eq. 3.4). What is left is Eq. 3.12, the equation the
# teammates' master equations are derived from, and the one solved here:
#
#     H(t) = f_x(t)/2 · σx  +  ξ(t)/2 · σz          (ħ = 1)
#            ─────────────     ──────────
#            H_q,0 (Eq. 3.4)   V(t), the noise
#
# ξ(t) is the OU process of Eq. 3.10, i.e. OUNoiseModel(σ, τ_c), or any other
# AbstractNoiseModel.
#
# Drive envelope. The thesis's χ(t) (Eq. 3.15) fixes it:
#     θ(t) = ∫₀ᵗ f_x = θ t/t_g − χ(t),   χ(t) = (θ/2π) sin(2πt/t_g)     (Eq. 3.7)
#  ⇔  f_x(t) = (θ/t_g)(1 − cos 2πt/t_g).
# The envelope switches on and off smoothly, and θ(t_g) = θ. With θ = π/2, the
# ideal gate is U_Ideal = (1 − iσx)/√2 (Eq. 3.8).
#
# Noiseless propagator (Eq. 3.6):  U_q,0(t) = cos(θ(t)/2)·1 − i sin(θ(t)/2)·σx.
#
# Interaction picture. Since i∂U_q,0 = H_q,0 U_q,0 (Eq. 3.5), the noise in the
# drive's frame is
#     V_I(t) = U_q,0† V U_q,0 = ξ(t)/2 · n(t)·σ,   n(t) = (0, sin θ(t), cos θ(t)).
# Eq. 3.13 prints U V U†, which flips the sign of n_y. See docs/04_quantum.md §7.

"""
    CosineGate(θ, tg)

The drive of thesis Ch. 3: a rotation by `θ` about x in gate time `tg`, with
envelope f_x(t) = (θ/tg)(1 − cos 2πt/tg). The thesis uses θ = π/2 (Eq. 3.8).
"""
struct CosineGate
    θ::Float64
    tg::Float64

    function CosineGate(θ::Real, tg::Real)
        tg > 0 || throw(ArgumentError("gate time tg must be > 0, got $tg"))
        return new(Float64(θ), Float64(tg))
    end
end

"Drive envelope f_x(t) = (θ/tg)(1 − cos 2πt/tg)."
drive(g::CosineGate, t::Real) = g.θ / g.tg * (1 - cos(2π * t / g.tg))

"Rotation angle θ(t) = ∫₀ᵗ f_x = θt/tg − χ(t), with χ(t) = (θ/2π) sin(2πt/tg) (Eqs. 3.7, 3.15)."
rotation_angle(g::CosineGate, t::Real) = g.θ * t / g.tg - g.θ / (2π) * sin(2π * t / g.tg)

"Noiseless propagator U_q,0(t) = cos(θ(t)/2)·1 − i sin(θ(t)/2)·σx (Eq. 3.6)."
function ideal_propagator(g::CosineGate, t::Real)
    a = rotation_angle(g, t)
    return cos(a / 2) * I2 - im * sin(a / 2) * σx
end

"The target gate U_Ideal = U_q,0(tg), which is (1 − iσx)/√2 for θ = π/2 (Eq. 3.8)."
ideal_gate(g::CosineGate) = ideal_propagator(g, g.tg)

"""
    hamiltonian(g, ξ, t) = f_x(t)/2·σx + ξ/2·σz

Thesis Eq. 3.12, the equation being solved. `ξ` is the value of the noise at time `t`.
"""
hamiltonian(g::CosineGate, ξ::Real, t::Real) = drive(g, t) / 2 * σx + ξ / 2 * σz

"""
    noise_axis(g, t) -> (0, sin θ(t), cos θ(t))

Unit vector n(t) with U_q,0†(t) σz U_q,0(t) = n(t)·σ.
"""
function noise_axis(g::CosineGate, t::Real)
    a = rotation_angle(g, t)
    return (0.0, sin(a), cos(a))
end

"""
    interaction_noise(g, ξ, t) = U_q,0† (ξ/2 σz) U_q,0 = ξ/2 · n(t)·σ

The noise Hamiltonian in the frame of the noiseless drive.
"""
function interaction_noise(g::CosineGate, ξ::Real, t::Real)
    n = noise_axis(g, t)
    return ξ / 2 * (n[2] * σy + n[3] * σz)
end