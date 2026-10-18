# =============================================================================
# Tutorial: two electrons and a nucleus (Z = 2) in a quasi-one-dimensional trap,
# a one-dimensional model of the helium atom.
#
#     H = -½ (∂²/∂z₁² + ∂²/∂z₂²) + λ [ v(z₁ - z₂; 1) - v(z₁; Z) - v(z₂; Z) ]
#
# The two electron coordinates (z₁, z₂) are the two dimensions of the mesh, so the
# correlated two-electron problem is solved exactly within the model. v is the
# effective interaction obtained by integrating out a harmonic transverse trap
# (see adhoc_potential_function.jl); λ ("switch factor") scales all interactions.
#
# 1. Runs three calculations with the FEMTISE input-file driver (results in ./output):
#      input_01.dat  without electron-electron repulsion, λ = 1
#      input_02.dat  full interaction, λ = 1
#      input_03.dat  full interaction, sweep of λ from 0 to 1
# 2. Exports the figures to ./figures (PDF and PNG).
#
# Usage (from the repository root):
#     julia --project=. src/coulomb_interaction_2d_potential/run_tutorial.jl [--rerun]
# =============================================================================

using FEMTISE, Gridap
include(joinpath(@__DIR__, "..", "..", "tools", "TutorialSimulation.jl"))
include(joinpath(@__DIR__, "..", "..", "tools", "TutorialFigures.jl"))
using .TutorialSimulation, .TutorialFigures

cd(@__DIR__)
set_publication_theme!()
rerun = "--rerun" in ARGS

# ---- 1. Simulations -------------------------------------------------------------
_, no_repulsion = simulate("input_01"; rerun = rerun)
_, full = simulate("input_02"; rerun = rerun)
_, sweep = simulate("input_03"; rerun = rerun)

ε_free, order_free = sorted_spectrum(no_repulsion)
ε_full, order_full = sorted_spectrum(full)
println("Ground state without / with repulsion: ", ε_free[1], " / ", ε_full[1])

# ---- 2a. Eigenvalues as functions of the interaction strength ---------------------------
λ = real.(sweep.λvector)
levels = sort(real.(sweep.ϵ_matrix), dims = 1)       # ascending order at every λ
fig = plot_parameter_sweep(λ, levels; states = 1:10, parameter_name = "Interaction strength",
                           parameter_symbol = "\\lambda", width = 330, height = 230)
save_figure(fig, "figures", "helium_model_sweep")

# ---- 2b. Two-electron densities ------------------------------------------------------------
z = collect(range(-4.0, 4.0, length = 161))
density_map(results, order, k) = abs2.(sample_on_grid(results.ϕ[order[k]], z, z))
densities = [density_map(no_repulsion, order_free, 1), density_map(no_repulsion, order_free, 2),
             density_map(full, order_full, 1), density_map(full, order_full, 2)]
titles = [math_label("\\mathrm{no\\ repulsion},\\ \\varepsilon_0=$(round(ε_free[1], digits = 4))"),
          math_label("\\mathrm{no\\ repulsion},\\ \\varepsilon_1=$(round(ε_free[2], digits = 4))"),
          math_label("\\mathrm{full},\\ \\varepsilon_0=$(round(ε_full[1], digits = 4))"),
          math_label("\\mathrm{full},\\ \\varepsilon_1=$(round(ε_full[2], digits = 4))")]
fig = plot_density_maps(z, z, densities, titles; columns = 2, xname = "Electron 1",
                        xsymbol = "z_1", yname = "Electron 2", ysymbol = "z_2", diagonal = true)
save_figure(fig, "figures", "helium_model_densities")
println("Figures written to ", joinpath(@__DIR__, "figures"))
