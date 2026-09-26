# Rebuild roadmap

Each step is one reviewable commit (or PR): code, tests, and an explainer in
`docs/`. A step is done only when CI is green and you can explain the explainer.

| Step | Scope | Key outputs | Status |
|------|-------|-------------|--------|
| 0 | Scaffold | `Project.toml`, CI (Julia 1.10 LTS + latest), Aqua, layout, `.gitignore` | ✅ done |
| 1 | Noise models | `NoiseEnsemble`, `OUNoiseModel` (exact AR(1)), `WhiteNoiseModel`, `FractionalNoiseModel` (FIR + AR), stationary theory | ✅ done |
| 1b | Band-limited 1/f | `BandLimitedOneOverFNoiseModel` (Mathematica kernel), FFT synthesis with exact variance | ✅ done |
| 2 | Estimators | FFT autocovariance with explicit `demean = :known / :ensemble / :trajectory`, periodogram, *corrected* Wiener–Khinchin, OLS, t quantiles, power-law and exponential fits | ✅ done |
| 3 | Exact finite-record expectations + validation | E[ACF estimator] and E[periodogram] for every model and α, `validate_noise_model` returning a report struct with pass/fail from standard errors | ⏳ |
| 4 | Quantum layer | Pauli algebra, U₀ frame, Liouvillian, averaged propagator Λ(t), states, gate fidelity, with analytic tests (zero noise → identity, static-noise closed form, trace preservation, complete positivity) | ⏳ |
| 5 | IO + plotting | CSV save/load round trip (the old `load_propagator` could not read `save_propagator` files), Plots as a package **extension** so the core does not depend on Plots | ⏳ |
| 6 | Notebook | `notebooks/Noise_Validation.ipynb` rebuilt section by section on the tested API, with every claim computed rather than typed; own `Project.toml` | ⏳ |
| 7 | Release | README polish, docs index, `v0.2.0` tag, archived figures/tables | ⏳ |

Steps 0–2: 286 tests pass locally (Julia 1.13).

## Decisions (2026-09-26)

1. **One package, two layers.** The noise library exists to drive the qubit
   simulation, so the noise models and the quantum layer (Physics, Propagator,
   gate fidelity) live together in this repo.
2. **Keep `BandLimitedOneOverFNoiseModel`**, which matches the Mathematica kernel
   Ci(100|τ|) − Ci(|τ|). It is ported in Step 1b and returns a `NoiseEnsemble`.
   `PowerLawNoiseModel` (superseded by `FractionalNoiseModel`) and the unfinished
   `TelegraphNoiseModel` stub are dropped.
3. OU generator: exact AR(1) update, i.e. Wikipedia's *Formal solution* applied per step
   (see `docs/01_noise_models.md` §2). Confirmed.
4. Theory follows the user's sources: Kasdin (1995), Dutta & Horn (1981),
   Ruseckas & Kaulakys (2010), Wikipedia OU. See `docs/REFERENCES.md`.