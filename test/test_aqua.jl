# Automated package-hygiene checks: method ambiguities, unbound type
# parameters, undefined exports, stale/missing dependencies, [compat] bounds,
# and type piracy. See https://github.com/JuliaTesting/Aqua.jl
using Aqua
Aqua.test_all(QuantumNoiseSimulator)
