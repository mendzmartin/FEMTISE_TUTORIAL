# =============================================================================
# Tutorial: eigenvalues over two potential parameters (tilted double well)
#
#     H = -½ d²/dx² + λ (x² − 1)² + μ x
#
# λ is the barrier height and μ the tilt. In the symmetric well (μ = 0) the
# lowest states form tunnelling doublets whose splitting ε₁ − ε₀ decreases
# exponentially with the barrier; a tilt localizes the states in one well and the
# gap grows to ≈ 2μ.
#
# The input file sweeps both parameters (FEMTISE ≥ v0.2.0):
#     analysis_param   = 1 0.75 10.0 0.25    λ = 1.0, 1.25, …, 10.0
#     analysis_param_2 = 2 -0.05 1.0 0.05    μ = 0.0, 0.05, …, 1.0
# (the driver computes λᵢ + kΔλ for k ≥ 1, hence the starting values).
#
# Usage (from the repository root):
#     julia --project=. src/double_well_parameter_map/run_tutorial.jl [--rerun]
# =============================================================================

using FEMTISE, Gridap
include(joinpath(@__DIR__, "..", "..", "tools", "TutorialSimulation.jl"))
include(joinpath(@__DIR__, "..", "..", "tools", "TutorialFigures.jl"))
using .TutorialSimulation, .TutorialFigures
using CairoMakie       # loads the plotting functions of FEMTISE

cd(@__DIR__)
set_publication_theme!()

# ---- 1. Two-parameter sweep (777 eigenvalue problems) -----------------------------
_, sweep = simulate("input"; rerun = "--rerun" in ARGS)
λ, μ, E = sweep.λvector, sweep.μvector, sweep.ϵ_tensor          # E[n, i, j]
splitting = E[2, :, :] .- E[1, :, :]
println("splitting at μ = 0: from ", splitting[1, 1], " (λ = ", λ[1], ") to ",
        splitting[end, 1], " (λ = ", λ[end], ")")

# ---- 2. Figures (FEMTISE plotting functions) -------------------------------------------------
fig = plot_energy_map(λ, μ, E; gap = (1, 2), colorscale = log10, contours = 8,
                      labels = (axis_label("Barrier height", "\\lambda", "a.u."),
                                axis_label("Tilt", "\\mu", "a.u.")))
save_figure(fig, "figures", "double_well_splitting_map")

fig = plot_energy_curves(λ, E[:, :, 1]; parameter_label = axis_label("Barrier height", "\\lambda", "a.u."),
                         energy_label = axis_label("Eigenvalue", "\\varepsilon_n", "a.u."))
save_figure(fig, "figures", "double_well_levels_symmetric")
println("Figures written to ", joinpath(@__DIR__, "figures"))
