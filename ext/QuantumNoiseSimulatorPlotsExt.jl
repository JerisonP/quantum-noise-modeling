# ─────────────────────────────────────────────────────────────────────────────
# Plots extension: loaded automatically when the user runs `using Plots`
# ─────────────────────────────────────────────────────────────────────────────
#
# Conventions (kept from the old Plotting.jl):
#   • every function RETURNS a plot and has no side effects. `display(p)` or
#     `savefig(p, "file.png")` is up to the caller;
#   • tab10 colours, line width 2 for primary series, dashed black for theory,
#     and dashed colours for comparison curves;
#   • error bands are ± 2 SE (about 95%).
# Plots is used only through `Plots.` so that none of its exported names can
# clash with the package's own.

module QuantumNoiseSimulatorPlotsExt

import Plots
using QuantumNoiseSimulator
import QuantumNoiseSimulator: plot_noise_traces, plot_validation, plot_state_evolution,
                              plot_timestep_convergence, plot_error_vs_delta,
                              plot_eigenvalues_vs_delta

const _LW = 2.0
_color(i) = Plots.palette(:tab10)[mod1(i, 10)]
_base(; kwargs...) = Plots.plot(; grid=true, titlefontsize=11, guidefontsize=10, tickfontsize=8, kwargs...)
_pairs(c) = c === nothing ? Pair[] : c isa Pair ? [c] : collect(c)

function plot_noise_traces(ens::NoiseEnsemble; n_show::Integer=5, title="Noise trajectories")
    p = _base(xlabel="t", ylabel="ξ(t)", title=title, size=(900, 400))
    for j in 1:min(n_show, ntrajectories(ens))
        Plots.plot!(p, times(ens), ens[j]; lw=1.0, alpha=0.8, color=_color(j),
                    label=(j == 1 ? "trajectories" : false))
    end
    Plots.hline!(p, [0.0]; color=:black, lw=0.6, ls=:dot, label=false)
    return p
end

function plot_validation(r::ValidationReport)
    a, s = r.acf, r.psd
    p1 = _base(xlabel="lag τ", ylabel="C(τ)", title="Autocovariance")
    Plots.plot!(p1, a.lags, a.C; ribbon=2 .* a.C_err, fillalpha=0.25, lw=_LW,
                color=_color(1), label="measured ± 2 SE")
    Plots.plot!(p1, a.lags, a.expected; lw=_LW, ls=:dash, color=:black, label="exact expectation")
    ok = (s.S .> 0) .& (s.expected .> 0)
    p2 = _base(xlabel="f", ylabel="S(f)", title="Periodogram", xscale=:log10, yscale=:log10)
    Plots.plot!(p2, s.freqs[ok], s.S[ok]; lw=1.0, color=_color(1), label="measured")
    Plots.plot!(p2, s.freqs[ok], s.expected[ok]; lw=_LW, ls=:dash, color=:black, label="exact expectation")
    status = passed(r) ? "PASS" : "FAIL"
    return Plots.plot(p1, p2; layout=(1, 2), size=(1000, 420),
                      plot_title="$(nameof(typeof(r.model))), M = $(r.M): $status",
                      plot_titlefontsize=12)
end

function plot_state_evolution(r::GateResult, ρ0::AbstractMatrix; title="Noise-averaged state")
    ev = evolve_state(r, ρ0)
    x = ev.times ./ r.gate.tg
    b = bloch_vector.(ev.ρ)
    p = _base(xlabel="t / t_g", ylabel="", title=title, size=(900, 420), legend=:outerright)
    for (k, name) in enumerate(("⟨σx⟩", "⟨σy⟩", "⟨σz⟩"))
        Plots.plot!(p, x, getindex.(b, k); lw=_LW, color=_color(k), label=name)
    end
    Plots.plot!(p, x, purity.(ev.ρ); lw=_LW, ls=:dash, color=:black, label="Tr ρ²")
    return p
end

function plot_timestep_convergence(rows::AbstractVector)
    dt = [r.dt for r in rows]
    p = _base(xlabel="dt", ylabel="Δε = ε(dt) − ε(finest dt)", xscale=:log10,
              title="Time-step convergence (paired, same trajectories)", size=(750, 450))
    εse = first(rows).ε_se
    Plots.hspan!(p, [-εse, εse]; color=:gray, alpha=0.2, label="± 1 Monte-Carlo SE of ε")
    Plots.hline!(p, [0.0]; color=:black, lw=0.6, label=false)
    Plots.plot!(p, dt, [r.Δε for r in rows]; yerror=2 .* [r.Δε_se for r in rows],
                marker=:circle, lw=_LW, color=_color(1), label="Δε ± 2 SE")
    return p
end

function plot_error_vs_delta(rows::AbstractVector; compare=nothing,
                             title="Average gate error vs δ")
    p = _base(xlabel="δ", ylabel="⟨ε⟩_ξ", yscale=:log10, title=title, size=(750, 500),
              legend=:bottomright)
    ok = [r.ε > 0 for r in rows]
    δ, ε, se = [r.δ for r in rows][ok], [r.ε for r in rows][ok], [r.ε_se for r in rows][ok]
    lower = min.(2 .* se, 0.99 .* ε)                   # keep the band above zero on a log axis
    Plots.plot!(p, δ, ε; ribbon=(lower, 2 .* se), fillalpha=0.25, lw=_LW, marker=:circle,
                markersize=3, color=_color(1), label="brute force ± 2 SE")
    for (k, (label, c)) in enumerate(_pairs(compare))
        okc = [x.ε > 0 for x in c]
        Plots.plot!(p, [x.δ for x in c][okc], [x.ε for x in c][okc]; lw=_LW, ls=:dash,
                    color=_color(k + 1), label=string(label))
    end
    return p
end

function plot_eigenvalues_vs_delta(rows::AbstractVector, state::Symbol; compare=nothing)
    state in keys(first(rows).states) ||
        throw(ArgumentError("state must be one of $(keys(first(rows).states)), got :$state"))
    cmp = _pairs(compare)
    δ = [r.δ for r in rows]
    st(r) = getproperty(r.states, state)
    panels = Any[]
    for (f, fname) in ((real, "Re"), (imag, "Im")), i in 1:2      # (a) Re λ₁ (b) Re λ₂ (c) Im λ₁ (d) Im λ₂
        sp = _base(xlabel="δ", ylabel="$fname λ$i", legend=:best)
        y = [f(st(r).λ[i]) for r in rows]
        if f === real
            Plots.plot!(sp, δ, y; ribbon=2 .* [st(r).λ_se[i] for r in rows], fillalpha=0.25,
                        lw=_LW, color=_color(1), label="brute force ± 2 SE")
        else
            Plots.plot!(sp, δ, y; lw=_LW, color=_color(1), label="brute force")
        end
        for (k, (label, c)) in enumerate(cmp)
            Plots.plot!(sp, [x.δ for x in c], [f(st(x).λ[i]) for x in c]; lw=_LW, ls=:dash,
                        color=_color(k + 1), label=string(label))
        end
        push!(panels, sp)
    end
    return Plots.plot(panels...; layout=(2, 2), size=(1000, 700),
                      plot_title="Eigenvalues of ⟨ρ_$state(t_g)⟩", plot_titlefontsize=12)
end

end # module