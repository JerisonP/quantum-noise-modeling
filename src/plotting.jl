# ─────────────────────────────────────────────────────────────────────────────
# Plotting entry points. The implementations live in ext/QuantumNoiseSimulatorPlotsExt.jl
# ─────────────────────────────────────────────────────────────────────────────
#
# Plots is a *weak dependency*: installing this package does not install Plots, and
# `using QuantumNoiseSimulator` does not load it. When a user also runs `using Plots`,
# Julia automatically loads the extension, which adds the real methods below.
# Without Plots, each function raises a clear error instead of a MethodError.

_needs_plots(name) = error("$name needs Plots: run `using Plots` first. (Plotting is a package " *
                           "extension, so the core library does not depend on Plots.)")

"""
    plot_noise_traces(ens::NoiseEnsemble; n_show=5, title) -> Plot

A few raw trajectories ξ(t), as a first look at the noise.
"""
plot_noise_traces(args...; kwargs...) = _needs_plots("plot_noise_traces")

"""
    plot_validation(report::ValidationReport) -> Plot

What `validate_noise_model` actually tested. Left: the measured autocovariance ± 2 SE
against its exact finite-record expectation. Right: the measured periodogram against
its exact expectation (log–log). The title shows PASS/FAIL.
"""
plot_validation(args...; kwargs...) = _needs_plots("plot_validation")

"""
    plot_state_evolution(r::GateResult, ρ0; title) -> Plot

⟨σx⟩, ⟨σy⟩, ⟨σz⟩ and the purity Tr ρ² of the noise-averaged state versus t/t_g.
"""
plot_state_evolution(args...; kwargs...) = _needs_plots("plot_state_evolution")

"""
    plot_timestep_convergence(rows) -> Plot

The output of `timestep_convergence`: the paired change Δε ± 2 SE versus dt, drawn
against a band of ±1 Monte-Carlo SE of ε. dt is converged when the points near the
smallest dt sit well inside the band.
"""
plot_timestep_convergence(args...; kwargs...) = _needs_plots("plot_timestep_convergence")

"""
    plot_error_vs_delta(rows; compare = ("label" => rows2, …), title) -> Plot

⟨ε⟩_ξ ± 2 SE versus δ on a log scale (thesis Figs. 3.20–3.22), from `delta_sweep`.
`compare` overlays other results with the same row shape, such as `tcl2_delta_sweep`
or a teammate's file read with `load_delta_sweep`. It takes one `label => rows` pair
or a vector of them.
"""
plot_error_vs_delta(args...; kwargs...) = _needs_plots("plot_error_vs_delta")

"""
    plot_eigenvalues_vs_delta(rows, state::Symbol; compare = ("label" => rows2, …)) -> Plot

The four panels of thesis Figs. 3.2–3.19 for one state (`:zp, :zm, :xp, :xm, :yp, :ym`):
Re λ₁, Re λ₂ (brute force ± 2 SE), Im λ₁, Im λ₂, each versus δ, with optional
comparison curves.
"""
plot_eigenvalues_vs_delta(args...; kwargs...) = _needs_plots("plot_eigenvalues_vs_delta")