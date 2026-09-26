# Rebuild roadmap

Each step is one reviewable commit (or PR): code, tests, and an explainer in
`docs/`. A step is done only when CI is green and you can explain the explainer.

| Step | Scope | Key outputs | Status |
|------|-------|-------------|--------|
| 0 | Scaffold | `Project.toml`, CI (Julia 1.10 LTS + latest), Aqua, layout, `.gitignore` | ✅ |
| 1 | Noise models | `NoiseEnsemble`, `OUNoiseModel` (exact AR(1)), `WhiteNoiseModel`, `FractionalNoiseModel` (FIR + AR), stationary theory | ✅ code + tests (awaiting first CI run) |
| 2 | Estimators | FFT autocovariance with explicit `mean = :known / :trajectory / :ensemble`, periodogram, *corrected* Wiener–Khinchin, OLS and log-binned slope fits | ⏳ |
| 3 | Exact finite-record expectations + validation | E[ACF estimator] and E[periodogram] for every model and α, `validate_noise_model` returning a report struct with pass/fail from standard errors | ⏳ |
| 4 | Quantum layer | Pauli algebra, U₀ frame, Liouvillian, averaged propagator Λ(t), states, gate fidelity, with analytic tests (zero noise → identity, static-noise closed form, trace preservation, complete positivity) | ⏳ |
| 5 | IO + plotting | CSV round trip (fixes B4), Plots as a package **extension** so the core does not depend on Plots | ⏳ |
| 6 | Notebook | `notebooks/Noise_Validation.ipynb` rebuilt section by section on the tested API, with every claim computed rather than typed; own `Project.toml` | ⏳ |
| 7 | Release | README polish, docs index, `v0.2.0` tag, archived figures/tables | ⏳ |

See `AUDIT.md` for the list of problems each step fixes.
