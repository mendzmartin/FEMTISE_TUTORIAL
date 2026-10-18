# =============================================================================
# Tutorial: vibrational levels of a diatomic molecule from a tabulated potential
#
#     H = -1/(2μ) d²/dR² + V(R)
#
# V(R) is read from the data file `pec.dat` (here the Morse potential of H₂, written
# by make_potential_table.jl) and μ = 918.08 mₑ is the reduced mass of the two
# protons. Both are given in the optional settings of the input files
# (FEMTISE ≥ v0.2.0):
#     masses              = 918.08
#     potential_data_file = ./pec.dat
#     interpolation       = linear | cubic
#
# The exact Morse levels, εₙ = ω(n + ½) − ω²(n + ½)²/(4D) with ω = a √(2D/μ), measure
# the accuracy of the two interpolations of the table.
#
# Usage (from the repository root):
#     julia --project=. src/tabulated_morse_vibrations/run_tutorial.jl [--rerun]
# =============================================================================

using FEMTISE, Gridap
include(joinpath(@__DIR__, "..", "..", "tools", "TutorialSimulation.jl"))
include(joinpath(@__DIR__, "..", "..", "tools", "TutorialFigures.jl"))
using .TutorialSimulation, .TutorialFigures
using CairoMakie

cd(@__DIR__)
set_publication_theme!()
rerun = "--rerun" in ARGS
isfile("pec.dat") || include("make_potential_table.jl")

# ---- 1. Simulations: the same table with linear and cubic interpolation ------------
_, linear = simulate("input_linear"; rerun = rerun)
id, cubic = simulate("input_cubic"; rerun = rerun)
ε_linear, _ = sorted_spectrum(linear)
ε_cubic, order = sorted_spectrum(cubic)
ϕ = cubic.ϕ[order]

# ---- 2. Exact Morse levels -------------------------------------------------------------
D, a, Rₑ = 0.1745, 1.0282, 1.4011          # must match make_potential_table.jl
μ = id.options.masses[1]
ω = a * sqrt(2D / μ)
n = 0:length(ε_cubic)-1
exact = [ω * (k + 0.5) - ω^2 * (k + 0.5)^2 / (4D) for k in n]
println("ω = ", ω, " hartree; largest relative error: linear ",
        maximum(abs.((ε_linear .- exact) ./ exact)), ", cubic ",
        maximum(abs.((ε_cubic .- exact) ./ exact)))

# ---- 3. Figures -------------------------------------------------------------------------
# (a) potential, levels and vibrational densities
R = collect(range(0.6, 3.0, length = 1201))
V = tabulated_potential("pec.dat"; interpolation = :cubic)
densities = [abs2.(sample_on_grid(ψ, R)) for ψ in ϕ[1:4]]
fig = plot_densities_1d(R, densities, ε_cubic; potential = V.(R), states = 1:4,
                        scale = 0.01, xlimits = (0.6, 3.0), ylimits = (0.0, 0.1),
                        position_name = "Internuclear distance", position_symbol = "R",
                        position_unit = "bohr", width = 400, height = 240)
save_figure(fig, "figures", "morse_levels")

# (b) accuracy of the two interpolations of the table
fig = Figure(size = (400, 230))
ax = Axis(fig[1, 1], yscale = log10, xlabel = axis_label("Vibrational quantum number", "n"),
          ylabel = axis_label("Relative error", "|\\varepsilon_n^h-\\varepsilon_n|/\\varepsilon_n"))
plots = []
for (k, ε) in enumerate((ε_linear, ε_cubic))
    s = series_style(k)
    error = (ε .- exact) ./ exact
    # filled markers: ε_h above the exact level; open markers: below
    push!(plots, scatterlines!(ax, n, abs.(error); color = s.color, marker = s.marker,
                               linestyle = s.linestyle, strokecolor = s.color, strokewidth = 0.9,
                               markercolor = [e > 0 ? s.color : RGBAf(1, 1, 1, 1) for e in error]))
end
bounds = [MarkerElement(color = :black, marker = :circle, markersize = 6),
          MarkerElement(color = :white, strokecolor = :black, strokewidth = 0.9, marker = :circle,
                        markersize = 6)]
methods = [[LineElement(color = series_style(k).color, linestyle = series_style(k).linestyle),
            MarkerElement(color = series_style(k).color, marker = series_style(k).marker, markersize = 6)]
           for k in 1:2]
Legend(fig[1, 2], [methods, bounds],
       [[math_label("\\mathrm{linear\\ interpolation}"), math_label("\\mathrm{cubic\\ spline}")],
        [math_label("\\varepsilon_n^h>\\varepsilon_n"), math_label("\\varepsilon_n^h<\\varepsilon_n")]],
       [math_label("\\mathrm{Table},\\ \\Delta R = 0.02\\ \\mathrm{bohr}"),
        math_label("\\mathrm{Sign\\ of\\ the\\ error}")])
save_figure(fig, "figures", "morse_table_accuracy")
println("Figures written to ", joinpath(@__DIR__, "figures"))
