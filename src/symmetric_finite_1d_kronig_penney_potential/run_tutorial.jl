# =============================================================================
# Tutorial: symmetric finite one-dimensional Kronig–Penney chain
#
#     H = -½ d²/dx² + V(x),   V = V₀ < 0 inside each of N (odd) wells of width b
#                              whose centres are separated by a; V = 0 elsewhere
#
# The bound levels of the N wells form a band of N states below the continuum
# threshold ε = 0 (a finite piece of the band structure of a periodic lattice).
#
# 1. Runs the FEMTISE input-file driver with `input.dat` (results in ./output).
# 2. Counts the bound states.
# 3. Exports the figures to ./figures (PDF and PNG).
#
# Usage (from the repository root):
#     julia --project=. src/symmetric_finite_1d_kronig_penney_potential/run_tutorial.jl [--rerun]
# =============================================================================

using FEMTISE, Gridap
include(joinpath(@__DIR__, "..", "..", "tools", "TutorialSimulation.jl"))
include(joinpath(@__DIR__, "..", "..", "tools", "TutorialFigures.jl"))
using .TutorialSimulation, .TutorialFigures
using CairoMakie

cd(@__DIR__)
include("adhoc_potential_function.jl")     # potential, also used to draw V(x)
set_publication_theme!()

# ---- 1. Simulation -------------------------------------------------------------
id, results = simulate("input"; rerun = "--rerun" in ARGS)
num_wells, period, width, depth = id.params.params_potential
ε, order = sorted_spectrum(results)
ϕ = results.ϕ[order]

# ---- 2. Bound states ----------------------------------------------------------------
bound = count(<(0.0), ε)
println("Bound states: ", bound, " (number of wells: ", num_wells, ")")

# ---- 3. Figures -----------------------------------------------------------------------
# (a) all computed levels and the continuum threshold; (b) zoom on the bound band.
n = 0:length(ε)-1
fig = Figure(size = (480, 215))
numeric = series_style(1)
ax1 = Axis(fig[1, 1], xlabel = axis_label("Quantum number", "n"),
           ylabel = axis_label("Eigenvalue", "\\varepsilon_n", "a.u."))
threshold = hlines!(ax1, [0.0]; color = :gray40, linestyle = :dash, linewidth = 0.8)
levels = scatter!(ax1, n, ε; color = numeric.color, marker = numeric.marker)
ax2 = Axis(fig[1, 2], xlabel = axis_label("Quantum number", "n"),
           ylabel = axis_label("Eigenvalue", "\\varepsilon_n", "a.u."),
           title = math_label("\\mathrm{bound\\ band},\\ $bound\\ \\mathrm{states}"))
scatterlines!(ax2, n[1:bound], ε[1:bound]; color = numeric.color, marker = numeric.marker,
              linestyle = numeric.linestyle)
Legend(fig[2, 1:2], [levels, threshold],
       [math_label("\\mathrm{FEMTISE\\ eigenvalues}"),
        math_label("\\mathrm{continuum\\ threshold},\\ \\varepsilon=0")];
       orientation = :horizontal, tellwidth = false)
TutorialFigures.panel_label!(fig[1, 1, TopLeft()], "a")
TutorialFigures.panel_label!(fig[1, 2, TopLeft()], "b")
save_figure(fig, "figures", "kronig_penney_spectrum")

L = id.params.L
x = collect(range(-L / 2, L / 2, length = 13801))     # 0.02 a.u. spacing
potential = [kronig_penney_potential((xi,), (num_wells, period, width, depth)) for xi in x]
states = (1, 2, bound)                               # bottom, next and top of the band
# evaluate only the states that are drawn (point evaluation of FE functions is costly)
densities = [n in states ? abs2.(sample_on_grid(ϕ[n], x)) : Float64[] for n in eachindex(ϕ)]
fig = plot_densities_stacked(x, densities, ε; potential = potential, states = states,
                             xlimits = (-L / 2, L / 2))
save_figure(fig, "figures", "kronig_penney_densities")
println("Figures written to ", joinpath(@__DIR__, "figures"))
