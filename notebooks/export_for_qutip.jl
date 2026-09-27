# Export noise paths + the package's per-trajectory results for an independent QuTiP check.
# Run from notebooks/:   julia --project=. export_for_qutip.jl
using QuantumNoiseSimulator, Random

θ, tg, τc = π / 2, 1.0, 1.0
σ = 0.965                                   # δ = 0.5, the strongest noise in Fig. 3.20
out = mkpath(joinpath(@__DIR__, "results"))

ens = generate_ensemble(OUNoiseModel(σ, τc), (0.0, tg), tg / 400, 300; rng=Xoshiro(11))
save_ensemble(joinpath(out, "qutip_noise.csv"), ens)

r = simulate_gate(CosineGate(θ, tg), ens)
ε = trajectory_errors(r)
open(joinpath(out, "package_out.csv"), "w") do io
    println(io, "# QuantumNoiseSimulator CSV v1\n# kind = check")
    println(io, "eps,zp_rho11,zp_rho22,zp_re_rho12,zp_im_rho12,xp_rho11,xp_rho22,xp_re_rho12,xp_im_rho12")
    for (e, S) in zip(ε, r.S_final)
        z = apply_channel(S, cardinal_states().zp)
        x = apply_channel(S, cardinal_states().xp)
        println(io, join(string.([e, real(z[1, 1]), real(z[2, 2]), real(z[1, 2]), imag(z[1, 2]),
                                  real(x[1, 1]), real(x[2, 2]), real(x[1, 2]), imag(x[1, 2])]), ","))
    end
end
println("wrote results/qutip_noise.csv and results/package_out.csv (300 trajectories, ε = ", sum(ε) / length(ε), ")")