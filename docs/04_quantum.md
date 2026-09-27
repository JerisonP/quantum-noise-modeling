# Step 4: Brute-force solution of the thesis model

*What equation is solved, how, how we know the answer is right, and how to compare
a master equation against it.* Equation numbers refer to the thesis
(Ch. 3 unless stated otherwise).

## 1. The job of this layer

The teammates derive approximate, time-local master equations for the
noise-averaged density matrix ⟨ρ(t)⟩: the SMNE (Eq. 2.57) and the PLME (Eq. 1.37),
each at 2nd and 4th order. This layer computes the same ⟨ρ(t)⟩ **without any
approximation in the noise**:

1. generate many exact noise trajectories ξ(t) (Steps 1–3 check they are right);
2. for each trajectory, solve the Schrödinger/von Neumann equation directly;
3. average.

There is no expansion in δ, so nothing breaks down as the noise gets strong. The
only errors are statistical (reported as a standard error on every number) and
time-step error (measured by `timestep_convergence`). That is why the result can
serve as the reference that a master equation is checked against.

## 2. The model: exactly the thesis equation

| Thesis | Equation | Code |
|--------|----------|------|
| Eq. 3.9 | H_q = ½[ω₀ + ξ(t)]σz + f_x(t)cos(ω_d t)σx (lab frame) | (starting point) |
| Eqs. 3.3–3.4, 3.12 | rotating frame, RWA, ω_d = ω_q: **H(t) = f_x(t)/2·σx + ξ(t)/2·σz** | `hamiltonian(g, ξ, t)` |
| Eq. 3.7 + 3.15 | θ(t) = ∫f_x = θt/t_g − χ(t), χ = (θ/2π)sin(2πt/t_g) ⇒ f_x = (θ/t_g)(1 − cos 2πt/t_g) | `rotation_angle`, `drive` |
| Eq. 3.6 | U_q,0(t) = cos(θ(t)/2)·1 − i sin(θ(t)/2)·σx | `ideal_propagator` |
| Eq. 3.8 | θ(t_g) = π/2 ⇒ U_Ideal = (1 − iσx)/√2 | `ideal_gate(CosineGate(π/2, t_g))` |
| Eq. 3.10 | dξ = −ξ/τ_c dt + √(2σ²/τ_c) dW | `OUNoiseModel(σ, τ_c)` |
| Eq. 3.21 | ⟨ρ(t_g)⟩ = L₀(t_g)⟨L_I(t_g)⟩ρ_j(0) | `final_state(r, ρ_j)` |
| §3.2.3 | the states zp, zm, xp, xm, yp, ym | `cardinal_states()` |
| Eq. 3.24 | F̄ = (1/6)Σ_j Tr[U_Ideal ρ_j U_Ideal† ρ_j(t_g)] | `fidelity_map`, `average_fidelity` |
| Eq. 3.25 | ⟨ε⟩_ξ = 1 − ⟨F̄⟩_ξ | `average_error` |
| Eq. 3.20 | δ = σ²τ_c/√(1 + θ²τ_c²/t_g²) | `small_noise_parameter` |

**Where the noise enters.** ξ(t) shifts the qubit splitting ω₀, so it sits next to ω₀
on σz/2 (Eq. 3.9) and stays there through the frame change and the RWA (Eq. 3.12). The
drive is on σx. `test_gate.jl` pins H element by element to
`[ξ/2 f/2; f/2 −ξ/2]`. `test_solvers.jl` shows that moving the noise to σx or σy
changes the answer by more than 0.05, so the tests would catch a misplaced noise term.

**The noise** is the OU process of Eq. 3.10: zero mean, stationary from t = 0, with
⟨ξ(t)ξ(t′)⟩ = σ²e^{−|t−t′|/τ_c}. It is generated exactly (Step 1), and
Step 3's `validate_noise_model` can be run on the very ensemble that drives the qubit.

## 3. How it is solved: two independent solvers

For one trajectory, ρ(t) = S(t)[ρ(0)] with the 4×4 superoperator obeying
dS/dt = −i(1 ⊗ H − Hᵀ ⊗ 1)S. Between samples, ξ is linearly interpolated.

| `solver` | Method | Why |
|----------|--------|-----|
| `:rk4` | classical RK4 directly on dS/dt, using only `hamiltonian` | solves the equation exactly as written; no interaction picture, no closed forms (the old code's approach, with a fixed, testable step) |
| `:magnus` (default) | interaction picture with 4th-order Magnus steps: U_step = exp(−i v·σ), v = (h/2)(a+b) − (√3/6)h²(a×b) | exactly unitary at every step, and faster |

The tests require the two to agree on the same trajectory to 10⁻⁷ (strong,
rough noise, ξ ~ 3). They also require both to agree with a third, independent
exponential-midpoint solver, and both to show 4th-order convergence (error ÷16 when
h halves). The averaged map ⟨S(t)⟩ is exactly L₀⟨L_I⟩ of Eq. 3.21.

**Time step.** Most of the time-step error comes from drawing a rough noise path as
straight lines, not from the integrator. `timestep_convergence` re-runs the gate on
the *same* trajectories at dt, 2dt, 4dt, … and reports the paired change Δε. Paired
differences have no Monte-Carlo noise from resampling, so Δε is measured precisely.
dt is converged when Δε at 2dt is far below the SE of ε. A NumPy check with
t_g = 100τ_c found biases of +0.4% at dt = τ_c/2 and +5.6% at dt = τ_c: dt must
resolve τ_c as well as t_g. `delta_sweep` therefore defaults to dt = min(t_g, τ_c)/1000.

## 4. What is reported, with errors

| Quantity | Function | Standard error from |
|----------|----------|---------------------|
| ⟨ρ_j(t_g)⟩ for the six states | `final_state(r, ρ0)` → `(ρ, se)` | the spread of per-trajectory ρ |
| eigenvalues λ₁ ≥ λ₂ (Figs. 3.2–3.19) | `final_state_eigenvalues(r, ρ0)` | batch means (the eigenvalue is nonlinear in ρ) |
| ⟨F̄⟩_ξ, ⟨ε⟩_ξ (Figs. 3.20–3.22) | `average_fidelity(r)`, `average_error(r)` | the spread of per-trajectory ε |
| physicality | `is_cptp(S)`, `validate_density_matrix(ρ)` | exact (no SE needed) |

**Comparing a master-equation result with brute force.** Take its ⟨ρ_j(t_g)⟩ or ε
at the same (t_g, τ_c, σ). If it lies within about 3–5 SE of brute force, the two
agree at that precision. If not, the difference is real, and its size is the error of
the master equation. Basis-independent quantities (ρ_j(t_g), its eigenvalues, ε) are
the safest to compare, because they do not depend on how a superoperator was vectorised.

The exact dynamics always gives a Hermitian ρ with real eigenvalues in [0, 1] and a
CPTP map (tested). So an imaginary eigenvalue, as in the 4th-order SMNE curves,
cannot be a feature of the physics: it comes from truncating the expansion.

## 5. δ and the reference curves

`small_noise_parameter` is Eq. 3.20 exactly. `delta_sweep(g, τ_c, δs, M)` produces
brute-force values for the thesis figures. At each δ it sets σ from Eq. 3.20 (varying
σ²τ_c at fixed t_g/τ_c, as the thesis does) and returns ε ± SE, and for each
cardinal state ⟨ρ_j(t_g)⟩ ± SE and λ₁,₂ ± SE. One unit-σ ensemble is scaled by σ
for every δ. This is exact because OU is linear in σ, and it makes the curve smooth.

**Units.** ξ is a rate, so σ²τ_c has units of 1/time. δ is a pure number only
after a time unit is fixed. For t_g = τ_c the choice does not matter; for t_g = τ_c/10
and τ_c/100 it does. Use the same time unit as the master-equation code.

## 6. The 2nd-order master equation (the thing being validated)

`tcl2_evolution(g, model)` solves
d⟨ρ_I⟩/dt = −∫₀ᵗ C(t−s)[V(t),[V(s),⟨ρ_I⟩]] ds with V = ½n(t)·σ. At 2nd order, the SMNE
and the PLME are identical (thesis §3.2.3), and this is that equation. It is **not** the
reference. It is here to (a) reproduce the teammates' 2nd-order curves from the same
objects, which confirms conventions before comparing, and (b) test the pipeline:

- at θ = 0 it is exact for Gaussian noise, and it matches the closed form to 10⁻⁶;
- at weak noise (δ = 0.02) it agrees with brute force within the SE.

Preliminary NumPy result (to be confirmed with Julia at larger M): at t_g = τ_c,
2nd order gives an ε about 3–5% above brute force for δ = 0.1–0.5. This agrees with the
thesis's observation that 2nd order overestimates the effect of the noise.

## 7. Things found in the thesis while checking the equations

| Where | As printed | Effect |
|-------|-----------|--------|
| Eq. 3.13 | V_I = U V U† | The interaction picture for i∂U = HU (Eq. 3.5) is U†VU (as in PLME Eq. 1.30). U V U† flips the sign of n_y, which swaps the ±x and ±y results; F̄ is unchanged. **Worth checking which one the master-equation code uses.** The brute force does not use an interaction picture, so it is unaffected. |
| Eq. 3.14 | ℒ_I(t) | As printed, its eigenvalues at t = 0 are ±i, ±i. A commutator −i[V_I, ·] must have two zero eigenvalues (and real parts 0). It is probably a typo, but **worth checking against the code**. |
| Eq. 3.17 | −(1 − e^{…})t_gτ_c/(t_g − iθτ_c) | The overall sign is wrong. No effect on δ, which uses the modulus. |
| Eq. 3.18 | +2e^{−t/τ_c}…cos(tθ/t_g) | It should be −2e^{−t/τ_c}…cos. No effect on δ, since that term is dropped for Eq. 3.19. |
| Eq. 3.10 | −ξ/τ_c | Missing dt. |
| Eq. 3.11 | ⟨ξξ⟩ = e^{−\|t−t₁\|/τ_c} | Missing σ² (Eq. 3.10 gives σ²e^{−\|τ\|/τ_c}). |

The corrected Eqs. 3.17–3.19 are checked by numerical integration in
`test_small_parameter.jl`.

## 8. Tests

| Test | What it proves |
|------|----------------|
| `hamiltonian == [ξ/2 f/2; f/2 −ξ/2]`, `ideal_gate == (1 − iσx)/√2` | Eqs. 3.12 and 3.8, element by element |
| rotation angle = θt/t_g − χ(t), dθ/dt = f_x, i∂U₀ = H₀U₀ | Eqs. 3.5–3.7, 3.15 |
| `:rk4` = `:magnus` to 10⁻⁷ on rough strong noise | two independent solvers of the same equation |
| both = a third solver; noise on σx/σy ⇒ differs by > 0.05 | the noise is on σz, and the test can tell |
| error ratio 16 when h halves, both solvers | 4th order, so the step can be controlled |
| θ = 0 OU: ⟨ρ₀₁⟩ = ½exp(−σ²τ_c(t − τ_c(1 − e^{−t/τ_c}))) within 5 SE | brute force reproduces an exact result for the thesis noise |
| CPTP, unital, all ρ valid, real eigenvalues | the exact average is physical |
| Eq. 3.24 = trace formula; per-trajectory mean = F̄ of mean | fidelity computed two ways |
| `timestep_convergence` paired Δε grows with dt; Δε_se ≪ ε_se | the convergence tool works |
| δ values, inverse, Eqs. 3.17–3.19 by quadrature | Eq. 3.20 and its derivation |
| TCL2: θ = 0 exact; weak noise = brute force | the reference and the approximation agree where they must |

## 9. Usage

```julia
using QuantumNoiseSimulator, Random

g  = CosineGate(π/2, 1.0)                          # thesis gate, t_g = 1
τc = 1.0                                           # t_g = τ_c
m  = OUNoiseModel(sigma_for_delta(0.1, τc, g.θ, g.tg), τc)
ens = generate_ensemble(m, (0.0, g.tg), 1e-3, 4000; rng = Xoshiro(1))

validate_noise_model(m, ens)                       # Step 3: the noise is right…
r = simulate_gate(g, ens)                          # …then solve the qubit with the SAME noise
average_error(r)                                   # ⟨ε⟩_ξ ± SE   (Eq. 3.25)
final_state(r, cardinal_states().xp)               # ⟨ρ_xp(t_g)⟩ ± SE (Eq. 3.21)
final_state_eigenvalues(r, cardinal_states().xp)   # λ₁, λ₂ ± SE (Fig. 3.2)
timestep_convergence(g, ens)                       # is dt small enough?

rows = delta_sweep(g, τc, 0.0:0.05:0.5, 4000; rng = Xoshiro(2))   # Fig. 3.20 reference
me = tcl2_evolution(g, m)                          # 2nd-order master equation, same inputs
```