# Step 5: Saving, loading and plotting

*How results leave the program (CSV), how a teammate's results come in, and how
everything is plotted without making Plots a hard dependency.*

## 1. One CSV layout for everything

```
# QuantumNoiseSimulator CSV v1        ← format line (checked on load)
# kind = delta_sweep                   ← what the file holds (checked on load)
# tau_c = 1.0                          ← any metadata, one "key = value" per line
delta,sigma,eps,eps_se,…               ← ONE header row of column names
0.02,0.1929,0.00331,0.000102,…         ← one row per record
```

The files open in Excel, pandas, Mathematica or a text editor. Two design rules make
the round trip safe, and the tests check both:

1. **Exact numbers.** Values are written with Julia's shortest exact representation,
   so loading a file gives back *bit-for-bit* the same Float64. `test_io.jl`
   compares with `==` (not `≈`) and includes −0.0, NaN, ±Inf, 5e-324 and 1e-300.
2. **Columns are found by name, never by position.** This fixes the old bug:
   `save_propagator` wrote `time, Λ…`, but `load_propagator` read column 1 as the
   parameter and column 2 as time, so every value came from the wrong column. The
   test scrambles the column order and requires an identical result.

A loader also checks `kind`, so loading a noise file as a channel is an error rather
than garbage.

| Function | `kind` | Columns | Use |
|----------|--------|---------|-----|
| `save_ensemble` / `load_ensemble` | `ensemble` | `t, xi_1 … xi_M` | archive or share the *exact* noise that drove a run |
| `save_channel` / `load_channel` | `channel` | `time, S11_re, S11_im, …, S44_im` | the averaged map ⟨S(t)⟩ (Eq. 3.21), from `simulate_gate` or `tcl2_evolution` |
| `save_delta_sweep` / `load_delta_sweep` | `delta_sweep` | `delta, sigma, eps, eps_se`, then per state `xp_rho11_re`, …, `xp_lambda1_re`, …, `xp_lambda1_se`, … | the reference curves of Figs. 3.2–3.22 |
| `read_table` | any | any | inspect a file: `(meta, cols, data)` |

`save_channel(path, r)` also records θ, t_g, the solver, the number of trajectories
and dt, so a file says how it was made.

## 2. The hand-in format: how a teammate's results come in

`load_delta_sweep` requires only `delta` and `eps`. **Every other column is optional**,
and a missing column loads as NaN. So a master-equation result can be handed in as

```
# QuantumNoiseSimulator CSV v1
# kind = delta_sweep
# source = SMNE 4th order, t_g = τ_c, time unit τ_c
delta,eps,xp_lambda1_re,xp_lambda1_im,xp_lambda2_re,xp_lambda2_im
0.1,0.0165,0.97,0.0,0.03,0.0
0.2,0.031,0.95,-0.001,0.05,0.001
```

It loads into exactly the same structure as `delta_sweep`, so it plots against
brute force with one line:

```julia
bf    = delta_sweep(g, τc, 0.0:0.05:0.5, 4000; rng = Xoshiro(1))
smne4 = load_delta_sweep("smne4.csv")
plot_error_vs_delta(bf; compare = ["SMNE 4th" => smne4, "2nd order" => tcl2_delta_sweep(g, τc, 0.0:0.05:0.5)])
plot_eigenvalues_vs_delta(bf, :xp; compare = ["SMNE 4th" => smne4])
```

Column names for the states: `zp, zm, xp, xm, yp, ym` (thesis §3.2.3), then
`_rhoab_re/_im` for element (a, b) of ⟨ρ(t_g)⟩ and `_lambda1/2_re/_im` for the eigenvalues.
They must use the same time unit and the same ρ convention (frame of Eq. 3.12, i.e.
the full ρ(t_g), not the interaction-picture one).

## 3. Plotting as a package extension

Plots is large and slow to load, and the library does not need it to compute anything.
So Plots is a **weak dependency** (`[weakdeps]` in `Project.toml`), and the plotting
code lives in `ext/QuantumNoiseSimulatorPlotsExt.jl`:

- `using QuantumNoiseSimulator` alone does not load Plots;
- `using QuantumNoiseSimulator, Plots` makes Julia load the extension automatically;
- calling a plot function without Plots gives a clear error ("run `using Plots` first").

| Function | Shows | Thesis figure |
|----------|-------|---------------|
| `plot_noise_traces(ens)` | a few raw ξ(t) | — |
| `plot_validation(report)` | measured C(τ) ± 2 SE and periodogram vs their exact expectations; PASS/FAIL | noise validation |
| `plot_state_evolution(r, ρ0)` | ⟨σx⟩, ⟨σy⟩, ⟨σz⟩, Tr ρ² vs t/t_g | — |
| `plot_timestep_convergence(rows)` | paired Δε ± 2 SE vs dt, with the ±1 SE band of ε | convergence check |
| `plot_error_vs_delta(rows; compare)` | ⟨ε⟩ ± 2 SE vs δ, log scale | Figs. 3.20–3.22 |
| `plot_eigenvalues_vs_delta(rows, :xp; compare)` | Re λ₁, Re λ₂ (± 2 SE), Im λ₁, Im λ₂ vs δ | Figs. 3.2–3.19 |

Every function *returns* the plot and has no side effects. Show it with `display(p)`,
save it with `savefig(p, "fig.png")`, and combine plots with `plot(p1, p2; layout = (2, 1))`.
Error bands are ± 2 SE (about 95%). Brute force is solid; comparisons are dashed.

## 4. Tests

| Test | What it proves |
|------|----------------|
| ensemble, channel and sweep files reload with `==` / `isequal` | the round trip is exact |
| awkward floats (−0.0, NaN, Inf, subnormal, 1e-300) survive | no precision or formatting loss |
| scrambled column order gives the same result | columns are read by name (old bug fixed) |
| TCL2 rows save and reload with NaN SEs | master-equation results fit the same format |
| a minimal `delta,eps,…` file loads | the hand-in format works |
| wrong kind, non-package file, missing `eps`, ragged row → `ArgumentError` | bad files fail loudly |
| a plot function before `using Plots` → clear error | the extension boundary works |
| after `using Plots`, the extension is loaded and all six plots build; one renders to PNG | the plotting code runs |