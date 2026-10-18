# =============================================================================
# Tutorial: autocorrelation function and spectrum of a wave packet
#
# A Gaussian wave packet is released in the anharmonic Morse potential
#     V(x) = D [1 − exp(−a x)]²,     D = 8, a = 0.5, mass 1   (atomic units)
# whose bound levels are exactly εₙ = ω(n + ½) − χ(n + ½)², ω = a√(2D) = 2,
# χ = ω²/(4D) = 0.125.
#
# From the eigenpairs computed with FEMTISE:
#   cₙ = ⟨φₙ|ψ₀⟩,   C(t) = Σ |cₙ|² e^{−iεₙt},   σ(ω) = Σ |cₙ|² δ(ω − εₙ).
# The survival probability |C(t)|² oscillates with the classical period ≈ 2π/ω and
# its pattern changes as the levels dephase (anharmonicity). All phases realign at
# t = 2π/χ; for this potential ω/χ = 16 is an even integer, so the packet already
# revives completely at t = π/χ.
#
# This tutorial uses the low-level interface (no input file).
#
# Usage (from the repository root):
#     julia --project=. src/autocorrelation_spectrum/run_tutorial.jl
# =============================================================================

using FEMTISE, Gridap
include(joinpath(@__DIR__, "..", "..", "tools", "TutorialFigures.jl"))
using .TutorialFigures
using CairoMakie

cd(@__DIR__)
set_publication_theme!()

# ---- 1. Eigenpairs of the Morse oscillator ---------------------------------------------
D, a = 8.0, 0.5
morse(x) = D * (1 - exp(-a * x[1]))^2
model = CartesianDiscreteModel((-4.0, 20.0), (1200,))
sol = solve_tise(model, morse; nev = 30, sigma = 0.0)

ω, χ = a * sqrt(2D), a^2 / 2           # χ = ω²/(4D)
bound = count(<(D), sol.energies)
exact = [ω * (n + 0.5) - χ * (n + 0.5)^2 for n in 0:bound-1]
println("bound states: ", bound, "; largest relative error: ",
        maximum(abs.((sol.energies[1:bound] .- exact) ./ exact)))

# ---- 2. Wave packet and its expansion coefficients ---------------------------------------
x₀ = 2.0                                   # displacement from the minimum
ψ₀(x) = (ω / π)^0.25 * exp(-ω * (x[1] - x₀)^2 / 2)
c = expansion_coefficients(ψ₀, sol.states, sol.dΩ)
weights = abs2.(c)
println("Σ|cₙ|² = ", sum(weights), " (bound states: ", sum(weights[1:bound]), ")")

# ---- 3. Autocorrelation and spectrum ------------------------------------------------------
t = collect(0.0:0.02:300.0)
C = autocorrelation(c, sol.energies, t)
revival_times = (π / χ, 2π / χ)

energy_grid = collect(range(0.0, 7.0, length = 3501))
γ = 0.05                                     # Gaussian window: resolution 0.05 a.u.
σ_from_C = fourier_transform(t, C, energy_grid; window = :gaussian, width = γ)
σ_direct = autocorrelation_spectrum(c, sol.energies, energy_grid; width = γ, lineshape = :gaussian)
println("max |σ(from C) − σ(direct)| = ", maximum(abs.(σ_from_C .- σ_direct)))

# ---- 4. Figure ------------------------------------------------------------------------------
fig = Figure(size = (500, 220))
ax1 = Axis(fig[1, 1], xlabel = axis_label("Time", "t", "a.u."),
           ylabel = axis_label("Survival probability", "|C(t)|^2"))
window = t .<= 65
lines!(ax1, t[window], abs2.(C[window]); color = series_style(1).color, linewidth = 0.8)
vlines!(ax1, collect(revival_times); color = :gray40, linestyle = :dash, linewidth = 0.8)
for (time, label) in zip(revival_times, ("\\pi/\\chi", "2\\pi/\\chi"))
    text!(ax1, time, 1.04; text = math_label(label), fontsize = 8, align = (:center, :bottom))
end
xlims!(ax1, 0, 65); ylims!(ax1, 0, 1.17)

ax2 = Axis(fig[1, 2], xlabel = axis_label("Energy", "\\omega", "a.u."),
           ylabel = axis_label("Spectrum", "\\sigma(\\omega)", "1/a.u."))
shown = (weights .> 1e-4) .& (sol.energies .< 7.0)
scale = maximum(σ_direct) / maximum(weights)          # sticks drawn on the same scale as σ
sticks = linesegments!(ax2, vec([Point2f(e, 0) for e in sol.energies[shown]] .=>
                                [Point2f(e, scale * w) for (e, w) in zip(sol.energies[shown], weights[shown])]);
                       color = :gray50, linewidth = 1.2)
from_C = lines!(ax2, energy_grid, σ_from_C; color = series_style(1).color)
direct = lines!(ax2, energy_grid, σ_direct; color = series_style(2).color, linestyle = :dash)
xlims!(ax2, 0, 7); ylims!(ax2, 0, 1.1 * maximum(σ_direct))
Legend(fig[2, 1:2], [from_C, direct, sticks],
       [math_label("\\mathrm{Fourier\\ transform\\ of}\\ C(t)"),
        math_label("\\mathrm{Gaussian\\ lines\\ at}\\ \\varepsilon_n"),
        math_label("|c_n|^2\\ \\mathrm{(scaled)}")];
       orientation = :horizontal, tellwidth = false)
TutorialFigures.panel_label!(fig[1, 1, TopLeft()], "a")
TutorialFigures.panel_label!(fig[1, 2, TopLeft()], "b")
save_figure(fig, "figures", "autocorrelation_spectrum")
println("Figures written to ", joinpath(@__DIR__, "figures"))
