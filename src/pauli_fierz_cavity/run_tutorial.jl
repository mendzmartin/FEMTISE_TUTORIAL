# =============================================================================
# Tutorial: nucleus + electron + cavity photon (Pauli–Fierz Hamiltonian, 3 DOF)
#
# One vibrational coordinate R (nucleus, mass M), one electronic coordinate z
# (mass 1) and one cavity mode with photon displacement coordinate q, in the
# length gauge (atomic units):
#
#     H = -1/(2M) ∂²_R - ½ ∂²_z - ½ ∂²_q + V(R, z, q)
#     V = V_m(R, z) + ½ ω_c² q² + g √(2ω_c) q d + (g²/ω_c) d²
#
# with V_m = ½ k_R R² + ½ k_z z² and dipole d = R − z. The last term is the dipole
# self-energy. This quadratic model is solved exactly by its normal modes, so the
# finite element levels can be checked level by level.
#
# When the cavity frequency ω_c crosses the electronic excitation ω_e = √k_z the
# electron and the photon hybridize into a lower and an upper polariton separated
# by the Rabi splitting Ω_Rabi ≈ 2 g |d₀₁| √(2ω_c)|q₀₁| = √2 g (here).
#
# This tutorial uses the low-level interface: solve_tise / sweep_eigenvalues on a
# three-dimensional mesh with one mass per coordinate.
#
# Usage (from the repository root):
#     julia --project=. src/pauli_fierz_cavity/run_tutorial.jl [--rerun]
# The sweeps (19 three-dimensional problems, a few minutes) are stored in
# output/pauli_fierz_sweeps.jld2 and reused unless --rerun is given.
# =============================================================================

using FEMTISE, Gridap, JLD2, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "tools", "TutorialFigures.jl"))
using .TutorialFigures
using CairoMakie

cd(@__DIR__)
set_publication_theme!()

# ---- 1. Model ---------------------------------------------------------------------------
M, k_R, k_z = 4.0, 0.81, 1.0             # nuclear frequency √(k_R/M) = 0.45, ω_e = 1
masses = (M, 1.0, 1.0)                   # coordinates (R, z, q)
matter(r) = 0.5k_R * r[1]^2 + 0.5k_z * r[2]^2
dipole(r) = r[1] - r[2]
cavity(ω_c, g) = pauli_fierz_potential(matter, dipole; frequency = ω_c, coupling = g)

"Exact normal-mode frequencies (ascending): eigenvalues of the mass-weighted Hessian."
function normal_modes(ω_c, g)
    shift = g * sqrt(2ω_c) / ω_c^2                  # V = ½ω_c²(q + shift·d)² + V_m
    u, e_q = [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]
    hessian = Diagonal([k_R, k_z, ω_c^2]) + ω_c^2 * shift^2 * u * u' +
              ω_c^2 * shift * (u * e_q' + e_q * u')
    weights = Diagonal(1 ./ sqrt.(collect(masses)))
    return sort(sqrt.(eigvals(Symmetric(weights * hessian * weights))))
end

# ---- 2. Finite element sweeps -------------------------------------------------------------
# Box (R, z, q) ∈ [−3.5, 3.5] × [−5, 5]², Lagrange P2 elements.
model = CartesianDiscreteModel((-3.5, 3.5, -5.0, 5.0, -5.0, 5.0), (14, 18, 18))
frequencies = collect(0.70:0.05:1.30)    # cavity frequency ω_c, coupling g = 0.05
g₀ = 0.05
couplings = collect(0.0:0.02:0.10)       # coupling g at resonance, ω_c = ω_e = 1
nev = 10

results = "output/pauli_fierz_sweeps.jld2"
if isfile(results) && !("--rerun" in ARGS)
    println("Using existing results in ", results, " (pass --rerun to recompute).")
    E_frequency, E_coupling = load(results, "E_frequency", "E_coupling")
else
    settings = (; masses = masses, nev = nev, sigma = 0.0)     # the lowest nev levels
    E_frequency = sweep_eigenvalues(model, (x, ω_c) -> cavity(ω_c, g₀)(x), frequencies; settings...)
    E_coupling = sweep_eigenvalues(model, (x, g) -> cavity(1.0, g)(x), couplings; settings...)
    mkpath("output")
    jldsave(results; E_frequency, E_coupling, frequencies, couplings)
end
excitations(E) = E .- E[1:1, :]          # εₙ − ε₀

# ---- 3. Accuracy -----------------------------------------------------------------------------
"Exact levels Σₖ nₖΩₖ (above the ground state), lowest `count`."
exact_levels(Ω, count) = sort(vec([dot(collect(n), Ω) for n in Iterators.product(0:6, 0:3, 0:3)]))[1:count]
# levels inside the window of the figure (εₙ − ε₀ < 1.5)
largest_error = maximum(maximum(abs(e - x) for (e, x) in zip(excitations(E_frequency)[:, i],
                                    exact_levels(normal_modes(ω_c, g₀), nev)) if x < 1.5)
                        for (i, ω_c) in enumerate(frequencies))
println("largest error of the excitation energies below 1.5 hartree: ", largest_error, " hartree")

# at resonance (ω_c = 1) the levels are 0, Ω_nuc, 2Ω_nuc, lower and upper polariton
rabi_splitting(Ω) = Ω[3] - Ω[2]
splitting_fem = [E[5] - E[4] for E in eachcol(E_coupling)]
splitting_exact = [rabi_splitting(normal_modes(1.0, g)) for g in couplings]
println("Rabi splitting, g = ", couplings, "\n  FEM   ", splitting_fem, "\n  exact ", splitting_exact)

# ---- 4. Figure ---------------------------------------------------------------------------------
fig = Figure(size = (500, 250))
window = (0.35, 1.5)                      # excitation energies shown in panel (a)
ax1 = Axis(fig[1, 1], xlabel = axis_label("Cavity frequency", "\\omega_c", "a.u."),
           ylabel = axis_label("Excitation energy", "\\varepsilon_n-\\varepsilon_0", "a.u."))
fine = collect(range(frequencies[1], frequencies[end], length = 241))
modes = [normal_modes(ω_c, g₀) for ω_c in fine]
# exact levels: every combination of quanta (nuclear, lower, upper); the two
# polariton fundamentals are highlighted, the rest drawn as thin gray lines
for n in Iterators.product(0:3, 0:1, 0:1)
    sum(n) == 0 && continue
    level = [dot(collect(n), Ω) for Ω in modes]
    minimum(level) < window[2] || continue
    n == (0, 1, 0) || n == (0, 0, 1) || lines!(ax1, fine, level; color = :gray60, linewidth = 0.7)
end
bare = lines!(ax1, fine, fine; color = :black, linestyle = :dot, linewidth = 0.9)
hlines!(ax1, [sqrt(k_z)]; color = :black, linestyle = :dot, linewidth = 0.9)
lower = lines!(ax1, fine, [Ω[2] for Ω in modes]; color = series_style(1).color, linewidth = 1.6)
upper = lines!(ax1, fine, [Ω[3] for Ω in modes]; color = series_style(2).color,
               linestyle = series_style(2).linestyle, linewidth = 1.6)
others = LineElement(color = :gray60, linewidth = 0.7)
shown = [(ω_c, e) for (i, ω_c) in enumerate(frequencies) for e in excitations(E_frequency)[2:end, i]
         if e < window[2] - 0.01]
fem = scatter!(ax1, first.(shown), last.(shown); color = :white, strokecolor = :black,
               strokewidth = 0.8, marker = :circle, markersize = 5)
xlims!(ax1, frequencies[1], frequencies[end]); ylims!(ax1, window...)
Legend(fig[2, 1], [lower, upper, others, bare, fem],
       [math_label("\\mathrm{lower\\ polariton}"), math_label("\\mathrm{upper\\ polariton}"),
        math_label("\\mathrm{other\\ levels}"), math_label("\\mathrm{bare}\\ \\omega_c,\\ \\omega_e"),
        math_label("\\mathrm{FEM}")];
       orientation = :horizontal, nbanks = 2, tellwidth = false)

ax2 = Axis(fig[1, 2], xlabel = axis_label("Coupling strength", "g", "a.u."),
           ylabel = axis_label("Rabi splitting", "\\Omega_\\mathrm{R}", "a.u."))
fine_g = collect(range(0, couplings[end], length = 101))
exact_line = lines!(ax2, fine_g, [rabi_splitting(normal_modes(1.0, g)) for g in fine_g];
                    color = :black, linewidth = 1.2)
estimate = lines!(ax2, fine_g, sqrt(2) .* fine_g; color = series_style(3).color,
                  linestyle = series_style(3).linestyle, linewidth = 1.2)
fem_g = scatter!(ax2, couplings, splitting_fem; color = :white, strokecolor = :black,
                 strokewidth = 0.8, marker = :circle, markersize = 5)
xlims!(ax2, -0.003, couplings[end] + 0.003); ylims!(ax2, -0.004, 1.06 * maximum(splitting_exact))
Legend(fig[2, 2], [exact_line, estimate, fem_g],
       [math_label("\\mathrm{normal\\ modes}"), math_label("\\sqrt{2}\\,g"), math_label("\\mathrm{FEM}")];
       orientation = :horizontal, nbanks = 2, tellwidth = false)
TutorialFigures.panel_label!(fig[1, 1, TopLeft()], "a")
TutorialFigures.panel_label!(fig[1, 2, TopLeft()], "b")
colsize!(fig.layout, 1, Relative(0.58))
save_figure(fig, "figures", "pauli_fierz_polaritons")
println("Figures written to ", joinpath(@__DIR__, "figures"))
