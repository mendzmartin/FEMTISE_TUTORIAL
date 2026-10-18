# =============================================================================
# Tutorial: isotropic one-dimensional harmonic oscillator
#
#     H = -½ d²/dx² + ½ ω² (x - x₀)²,     exact levels εₙ = ω (n + ½)
#
# 1. Runs the FEMTISE input-file driver with `input.dat` and the potential defined
#    in `adhoc_potential_function.jl` (results in ./output).
# 2. Compares the spectrum with the exact levels.
# 3. Exports the figures to ./figures (PDF and PNG).
#
# Usage (from the repository root):
#     julia --project=. src/isotropic_1d_harmonic_oscillator_potential/run_tutorial.jl [--rerun]
# =============================================================================

using FEMTISE, Gridap
include(joinpath(@__DIR__, "..", "..", "tools", "TutorialSimulation.jl"))
include(joinpath(@__DIR__, "..", "..", "tools", "TutorialFigures.jl"))
using .TutorialSimulation, .TutorialFigures

cd(@__DIR__)
set_publication_theme!()

# ---- 1. Simulation -------------------------------------------------------------
id, results = simulate("input"; rerun = "--rerun" in ARGS)
ω, x₀ = id.params.params_potential
ε, order = sorted_spectrum(results)
ϕ = results.ϕ[order]

# ---- 2. Comparison with the exact spectrum ----------------------------------------
exact = [ω * (n + 0.5) for n in 0:length(ε)-1]
relative_error = (ε .- exact) ./ exact
println("Largest relative error of the ", length(ε), " levels: ", maximum(abs.(relative_error)))

# ---- 3. Figures -----------------------------------------------------------------------
fig = plot_spectrum(ε; exact = exact)
save_figure(fig, "figures", "qho1d_spectrum")

x = collect(range(-6.0, 6.0, length = 1201))
densities = [abs2.(sample_on_grid(ψ, x)) for ψ in ϕ]
potential = @. 0.5 * ω^2 * (x - x₀)^2
fig = plot_densities_1d(x, densities, ε; potential = potential, states = 1:4,
                        scale = 1.5, xlimits = (-6, 6), ylimits = (0, 4.6))
save_figure(fig, "figures", "qho1d_densities")
println("Figures written to ", joinpath(@__DIR__, "figures"))
