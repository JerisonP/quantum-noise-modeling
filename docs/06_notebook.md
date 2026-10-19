# Step 6: The notebook

`notebooks/Noise_Validation.ipynb` is the numerics appendix. Every number and figure in it is
computed by the tested library, and nothing is typed in by hand.

## How to run it

```
cd notebooks
julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'    # once; installs Plots
jupyter notebook Noise_Validation.ipynb                        # then "Run All"
```

The notebook has its own environment, `notebooks/Project.toml`, which points at the package in
the folder above (`[sources]`, Julia ≥ 1.11). On Julia 1.10, run
`Pkg.develop(path="..")` in that environment instead. Commit `notebooks/Manifest.toml` after the
first run: it pins the exact package versions the figures were made with.

| Switch (cell §0) | Trajectories per gate run | Run time (rough) |
|------------------|---------------------------|------------------|
| `FAST = true`    | 1 000                     | about 1 minute   |
| `FAST = false`   | 10 000                    | a few minutes    |

Outputs: `notebooks/figures/*.png` and `notebooks/results/*.csv`.

## What each section establishes

| Section | Claim | How it is shown |
|---------|-------|-----------------|
| 0 | parameters are the thesis's | τ_c = 1 (the time unit), θ = π/2, t_g ∈ {1, 0.1, 0.01}, δ ranges of Figs. 3.20–3.22 |
| 1 | 1/f (α = 1) generator is right, including its nonstationarity | `validate_noise_model` + exact Var ξ(t); slope β̂ vs the same fit on the exact E[S] for α = 0.5, 1, 1.5 |
| 2 | band-limited 1/f (Mathematica kernel) is right | `validate_noise_model` |
| 3 | OU (Eq. 3.10) is right, and exact at any Δt | validation, Gaussian marginal, τ̂_c, and a PASS at Δt = 0.5, 0.1, 0.01 |
| 4 | white noise is right | `validate_noise_model` |
| 6 | the solver solves Eq. 3.12 | H printed, U_Ideal = (1 − iσx)/√2, two independent solvers agree, measured 4th order |
| 7 | brute force reproduces an exact result | θ = 0 dephasing law, and the 2nd-order equation is exact there too |
| 8 | the numbers are converged and physical | noise driving the gate validated; Magnus = RK4 on the same trajectories; F̄ two ways; CPTP; valid ρ for all 6 states; paired time-step check at every t_g |
| 9 | brute-force reference curves | ⟨ε⟩ ± SE vs δ (Figs. 3.20–3.22) and eigenvalue panels for all 6 states (Figs. 3.2–3.19), with the 2nd-order master equation overlaid, and a table of relative differences and z-scores |
| 10 | teammates get the numbers | one CSV per t_g, with parameters and seed; the files reload bit-for-bit |
| 11 | summary | collected from the cells above |

## What to look for when it runs

* **Part I:** every report says PASS. With about 10 reports × 4 checks at a 10⁻³ false-alarm
  rate each, an occasional single FAIL can happen by chance. Rerun that section with another seed
  before suspecting the generator.
* **§8 time step:** at the production step (factor 8, dt = t_g/400), |Δε| should be well inside
  the ± 1 SE band.
* **§9 2nd order vs brute force:** at small δ the relative difference should go to 0. A difference
  that grows steadily with δ is the breakdown of 2nd order. At t_g = τ_c the notebook finds
  +0.3% to +1.7% for δ = 0.1–0.5.
* **Convention check:** at the right edge, the 2nd-order ε should be about 8.8 × 10⁻⁴ (t_g = τ_c/10,
  δ = 0.05) and 8.9 × 10⁻⁶ (t_g = τ_c/100, δ = 0.005). If the thesis figures show those values, the
  teammates' conventions match (τ_c = 1).

## External cross-validation (QuTiP)

`notebooks/export_for_qutip.jl` writes 300 OU noise paths (t_g = τ_c = 1, θ = π/2, δ = 0.5) and
the package's per-trajectory results. `notebooks/qutip_check.py` solves the same Hamiltonian,
H(t) = f_x(t)/2·σx + ξ(t)/2·σz, with QuTiP's own operators and ODE solver on the identical paths.
Result (QuTiP 5.3.1): largest per-trajectory |Δε| = 1.7 × 10⁻⁸; largest difference in ⟨ρ(t_g)⟩
= 4.5 × 10⁻⁹ (from |0⟩) and 4.2 × 10⁻⁹ (from |+⟩). Moving the noise to σx in the QuTiP script
changes the results by ~10⁻¹, so the check discriminates. The output is saved in
`notebooks/results/qutip_check.txt`. Rerun:

    julia --project=. export_for_qutip.jl
    uv run --python 3.12 --with qutip --with numpy python qutip_check.py