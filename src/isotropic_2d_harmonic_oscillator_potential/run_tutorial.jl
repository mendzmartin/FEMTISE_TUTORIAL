# =============================================================================
# Tutorial: isotropic two-dimensional harmonic oscillator
#
#     H = -½ (∂²ₓ + ∂²ᵧ) + ½ ω² [(x - x₀)² + (y - y₀)²]
#     exact levels ε = ω (n + 1),  n = nₓ + nᵧ,  degeneracy n + 1
#
# 1. Runs the FEMTISE input-file driver with `input.dat` (results in ./output).
# 2. Compares the spectrum with the exact, degenerate levels.
# 3. Exports the figures to ./figures (PDF and PNG).
#
# Usage (from the repository root):
#     julia --project=. src/isotropic_2d_harmonic_oscillator_potential/run_tutorial.jl [--rerun]
# =============================================================================

using FEMTISE, Gridap
include(joinpath(@__DIR__, "..", "..", "tools", "TutorialSimulation.jl"))
include(joinpath(@__DIR__, "..", "..", "tools", "TutorialFigures.jl"))
using .TutorialSimulation, .TutorialFigures

cd(@__DIR__)
set_publication_theme!()

# ---- 1. Simulation -------------------------------------------------------------
id, results = simulate("input"; rerun = "--rerun" in ARGS)
ω = id.params.params_potential[1]
ε, order = sorted_spectrum(results)
ϕ = results.ϕ[order]

# ---- 2. Comparison with the exact spectrum ----------------------------------------
exact = sort([ω * (nx + ny + 1) for nx in 0:10 for ny in 0:10])[1:length(ε)]
println("Largest relative error of the ", length(ε), " levels: ",
        maximum(abs.((ε .- exact) ./ exact)))

# ---- 3. Figures -----------------------------------------------------------------------
fig = plot_spectrum(ε; exact = exact, index_offset = 1, index_name = "State index",
                    index_symbol = "k")
save_figure(fig, "figures", "qho2d_spectrum")

# Densities of the six lowest states. Within a degenerate level any combination of
# eigenvectors is also an eigenvector, so the orientation of the nodal lines of
# degenerate states is arbitrary.
grid = collect(range(-4.0, 4.0, length = 161))
states = 1:6
densities = [abs2.(sample_on_grid(ϕ[k], grid, grid)) for k in states]
titles = [math_label("k=$k,\\ \\varepsilon_k=$(round(ε[k], digits = 4))") for k in states]
fig = plot_density_maps(grid, grid, densities, titles; columns = 3, xsymbol = "x", ysymbol = "y")
save_figure(fig, "figures", "qho2d_densities")
println("Figures written to ", joinpath(@__DIR__, "figures"))
