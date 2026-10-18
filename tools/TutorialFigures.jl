"""
    TutorialFigures

Publication-quality figures for the FEMTISE.jl tutorials, built with CairoMakie.

Every figure follows the same rules:

* **LaTeX typography** for labels, legends and tick numbers.
* **Axis labels** give a descriptive name, the symbol and the unit:
  "Energy, ε (a.u.)" — see [`axis_label`](@ref).
* **Black-and-white and colour-blind safe**: each data series is identified by
  colour (Okabe–Ito palette) *and* line style *and* marker shape, so it remains
  distinguishable in grey-scale print and for readers with colour-vision deficiency.
* **Density maps** use `:inferno`, a perceptually uniform colour map whose lightness
  increases monotonically (it prints as a clean grey ramp), always with a colour bar.
* **No text on top of data**: legends are placed outside the data region.
* **Compact**: axis ranges are fitted to the data.

Figures are exported with [`save_figure`](@ref) as vector PDF (for articles) and
300-dpi PNG (for web pages).
"""
module TutorialFigures

using CairoMakie

const LaTeXString = CairoMakie.Makie.LaTeXStrings.LaTeXString

export set_publication_theme!, axis_label, math_label, series_style, save_figure,
       plot_spectrum, plot_densities_1d, plot_densities_stacked, plot_density_maps,
       plot_parameter_sweep

# -----------------------------------------------------------------------------
# Style
# -----------------------------------------------------------------------------

"Okabe–Ito colour-blind-safe palette."
const PALETTE = [colorant"#0072B2", colorant"#D55E00", colorant"#009E73",
                 colorant"#CC79A7", colorant"#E69F00", colorant"#56B4E9",
                 colorant"#000000"]
const LINE_STYLES = [:solid, :dash, :dot, :dashdot, (:dash, :dense), (:dot, :loose), :solid]
const MARKERS = [:circle, :rect, :diamond, :utriangle, :dtriangle, :pentagon, :xcross]

"Colour map used for every density plot."
const DENSITY_COLORMAP = :inferno

"""
    series_style(k)

Colour, line style and marker of the `k`-th data series. Using all three together
keeps the series distinguishable without colour.
"""
series_style(k::Integer) = (color = PALETTE[mod1(k, 7)],
                            linestyle = LINE_STYLES[mod1(k, 7)],
                            marker = MARKERS[mod1(k, 7)])

"""
    set_publication_theme!()

Activate the common theme: LaTeX fonts, 8-pt text, thin axes and frames.
Figure sizes are given in points (1 unit = 1 pt).
"""
function set_publication_theme!()
    set_theme!(merge(Theme(
        fontsize = 9,
        figure_padding = (4, 10, 4, 4),
        Axis = (; xgridvisible = false, ygridvisible = false,
                  xtickalign = 1, ytickalign = 1, xticksize = 3, yticksize = 3,
                  spinewidth = 0.7, xtickwidth = 0.7, ytickwidth = 0.7, titlesize = 9),
        Lines = (; linewidth = 1.2),
        Scatter = (; markersize = 5, strokewidth = 0.7),
        Legend = (; framevisible = false, labelsize = 8, patchsize = (18, 8), rowgap = 0),
        Colorbar = (; ticklabelsize = 8, labelsize = 9, size = 7, spinewidth = 0.6),
    ), theme_latexfonts()))
end

"""
    axis_label(name, symbol = "", unit = "")

LaTeX axis label "Name, symbol (unit)". `symbol` is written in LaTeX math syntax,
`name` and `unit` in plain text.

```julia
axis_label("Energy", raw"\\varepsilon_n", "a.u.")   # Energy, εₙ (a.u.)
```
"""
function axis_label(name::AbstractString, symbol::AbstractString = "", unit::AbstractString = "")
    text = "\\mathrm{" * replace(name, " " => "\\ ") * "}"
    isempty(symbol) || (text *= ",\\ " * symbol)
    isempty(unit) || (text *= "\\ \\mathrm{(" * replace(unit, " " => "\\ ") * ")}")
    return LaTeXString("\$" * text * "\$")
end

"""
    math_label(text)

LaTeX math-mode string for titles and legend entries, e.g.
`math_label(raw"\\varepsilon_0 = 0.5")`.
"""
math_label(text::AbstractString) = LaTeXString("\$" * text * "\$")

"""
    log_axis_ticks(values)

Tick positions and LaTeX labels for a logarithmic axis that covers `values`, and
the matching axis limits. Ticks are placed at powers of ten; when the data span less
than three decades, ticks at 2× and 5× are added so that the axis is not empty.
"""
function log_axis_ticks(values)
    vmin, vmax = extrema(filter(>(0), values))
    lower, upper = vmin / 1.6, vmax * 1.6
    k_min, k_max = floor(Int, log10(lower)), ceil(Int, log10(upper))
    mantissas = k_max - k_min <= 3 ? (1, 2, 5) : (1,)
    ticks = [m * 10.0^k for k in k_min:k_max for m in mantissas if lower <= m * 10.0^k <= upper]
    label(t) = (k = floor(Int, log10(t) + 1e-9); m = round(Int, t / 10.0^k);
                m == 1 ? math_label("10^{$k}") : math_label("$m\\times10^{$k}"))
    return (ticks, label.(ticks)), (lower, upper)
end

"Panel label (a), (b), … placed above the top-left corner of a grid position."
panel_label!(position, letter) =
    Label(position, LaTeXString("\$\\mathrm{($letter)}\$"); fontsize = 10,
          halign = :right, valign = :bottom, padding = (0, 6, 4, 0), tellwidth = false)

"""
    save_figure(fig, folder, name; png_dpi = 300)

Write `fig` to `folder/name.pdf` (vector) and `folder/name.png` (raster at
`png_dpi`). Returns the two paths.
"""
function save_figure(fig, folder::AbstractString, name::AbstractString; png_dpi::Real = 300)
    mkpath(folder)
    pdf = joinpath(folder, name * ".pdf")
    png = joinpath(folder, name * ".png")
    save(pdf, fig; pt_per_unit = 1)
    save(png, fig; px_per_unit = png_dpi / 72)
    return pdf, png
end

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------

"""
    plot_spectrum(energies; exact = nothing, thresholds = (), index_offset = 0,
                  index_name = "Quantum number", index_symbol = "n",
                  width = 480, height = 190)

Eigenvalues versus their index. If the exact eigenvalues are given, they are drawn
as a reference and a second panel shows the relative error
``(\\varepsilon_n^h - \\varepsilon_n)/\\varepsilon_n`` (filled markers: upper bounds,
open markers: lower bounds). `index_offset` is the index of the first eigenvalue
(0 for quantum numbers n = 0, 1, …; 1 for state indices k = 1, 2, …).
`thresholds` is a collection of `(energy, label)` pairs drawn as horizontal
reference lines (e.g. the continuum threshold).
"""
function plot_spectrum(energies::AbstractVector{<:Real}; exact = nothing, thresholds = (),
                       index_offset::Integer = 0,
                       index_name = "Quantum number", index_symbol = "n",
                       width = 480, height = 190)
    n = index_offset .+ (0:length(energies)-1)
    fig = Figure(size = isnothing(exact) ? (width ÷ 2 + 40, height) : (width, height))
    numeric, reference = series_style(1), series_style(7)

    ax1 = Axis(fig[1, 1], xlabel = axis_label(index_name, index_symbol),
               ylabel = axis_label("Eigenvalue", "\\varepsilon_{$index_symbol}", "a.u."))
    plots, labels = Any[], Any[]
    if !isnothing(exact)
        push!(plots, stairs!(ax1, n, exact; step = :center, color = reference.color, linewidth = 0.8))
        push!(labels, L"\mathrm{exact}")
    end
    for (k, (level, label)) in enumerate(thresholds)
        push!(plots, hlines!(ax1, [level]; color = :gray40, linewidth = 0.8,
                             linestyle = LINE_STYLES[mod1(k + 1, 7)]))
        push!(labels, label)
    end
    push!(plots, scatter!(ax1, n, energies; color = numeric.color, marker = numeric.marker))
    push!(labels, L"\mathrm{FEMTISE}")

    if isnothing(exact)
        Legend(fig[1, 2], plots, labels)
        return fig
    end
    relative_error = (energies .- exact) ./ abs.(exact)
    magnitudes = max.(abs.(relative_error), 1e-16)
    yticks, ylimits = log_axis_ticks(magnitudes)
    ax2 = Axis(fig[1, 2], yscale = log10, yticks = yticks, xlabel = axis_label(index_name, index_symbol),
               ylabel = axis_label("Relative error",
                                   "|\\varepsilon_{$index_symbol}^h-\\varepsilon_{$index_symbol}|/\\varepsilon_{$index_symbol}"))
    upper = relative_error .> 0
    ylims!(ax2, ylimits...)
    scatterlines!(ax2, n, magnitudes; color = numeric.color,
                  marker = numeric.marker, linestyle = numeric.linestyle,
                  markercolor = [u ? numeric.color : RGBAf(1, 1, 1, 1) for u in upper],
                  strokecolor = numeric.color)
    Legend(fig[1, 3], [plots...,
                       MarkerElement(color = numeric.color, marker = :circle, markersize = 6,
                                     strokecolor = numeric.color),
                       MarkerElement(color = :white, marker = :circle, markersize = 6,
                                     strokecolor = numeric.color, strokewidth = 0.9)],
           [labels..., L"\mathrm{upper\ bound}", L"\mathrm{lower\ bound}"])
    panel_label!(fig[1, 1, TopLeft()], "a")
    panel_label!(fig[1, 2, TopLeft()], "b")
    return fig
end

"""
    plot_densities_1d(x, densities, energies; potential = nothing, states = 1:3,
                      scale = 1.0, xlimits = extrema(x), ylimits = nothing,
                      energy_unit = "a.u.", position_unit = "a.u.",
                      position_name = "Position", position_symbol = "x",
                      width = 480, height = 220)

Probability densities ``\\rho_n(x)`` drawn at the height of their energies,
``\\varepsilon_n + c\\,\\rho_n(x)`` (`scale` = ``c``), on top of the potential
``V(x)`` (gray). `densities` is a vector of vectors sampled at `x`.
"""
function plot_densities_1d(x, densities, energies; potential = nothing, states = 1:3,
                           scale::Real = 1.0, xlimits = extrema(x), ylimits = nothing,
                           energy_unit = "a.u.", position_unit = "a.u.",
                           position_name = "Position", position_symbol = "x",
                           width = 480, height = 220)
    fig = Figure(size = (width, height))
    q = position_symbol
    ax = Axis(fig[1, 1], xlabel = axis_label(position_name, q, position_unit),
              ylabel = axis_label("Energy", "V($q),\\ \\varepsilon_n+c\\,\\rho_n($q)", energy_unit))
    plots, labels = Any[], Any[]
    if !isnothing(potential)
        # thick light-gray line: distinguishable from the states also in grey scale
        push!(plots, lines!(ax, x, potential; color = :gray75, linewidth = 2.6))
        push!(labels, math_label("V($q)"))
    end
    for (k, n) in enumerate(states)
        s = series_style(k)
        hlines!(ax, [energies[n]]; color = (s.color, 0.6), linewidth = 0.5, linestyle = :dot)
        push!(plots, lines!(ax, x, energies[n] .+ scale .* densities[n];
                            color = s.color, linestyle = s.linestyle))
        push!(labels, LaTeXString("\$n=$(n - 1)\$"))
    end
    xlims!(ax, xlimits...)
    isnothing(ylimits) || ylims!(ax, ylimits...)
    Legend(fig[1, 2], plots, labels)
    return fig
end

"""
    plot_densities_stacked(x, densities, energies; potential, states = 1:3,
                           xlimits = extrema(x), position_unit = "a.u.",
                           energy_unit = "a.u.", width = 480, panel_height = 55)

Potential ``V(x)`` in the top panel and, below it, one panel per state with its
probability density ``\rho_n(x)``; all panels share the position axis. Use this
instead of [`plot_densities_1d`](@ref) when the levels are too close in energy to be
drawn at their own heights (e.g. the narrow bands of a lattice). The energies of
the states are marked on the potential panel and written in each density panel.
"""
function plot_densities_stacked(x, densities, energies; potential, states = 1:3,
                                xlimits = extrema(x), position_unit = "a.u.",
                                energy_unit = "a.u.", width = 480, panel_height = 55)
    rows = length(states) + 1
    fig = Figure(size = (width, (rows + 1) * panel_height + 70))
    top = Axis(fig[1, 1], ylabel = axis_label("Potential", "V", energy_unit))
    lines!(top, x, potential; color = :gray40, linewidth = 0.8)
    for (k, n) in enumerate(states)
        s = series_style(k)
        hlines!(top, [energies[n]]; color = s.color, linestyle = s.linestyle, linewidth = 0.8)
    end
    hidexdecorations!(top, ticks = false)
    axes = [top]
    for (k, n) in enumerate(states)
        s = series_style(k)
        ax = Axis(fig[k + 1, 1], ylabel = math_label("\\rho_{$(n - 1)}"),
                  yticks = WilkinsonTicks(2))
        lines!(ax, x, densities[n]; color = s.color, linestyle = s.linestyle)
        text!(ax, 0.99, 0.95; space = :relative, align = (:right, :top), fontsize = 8,
              text = math_label("n=$(n - 1),\\ \\varepsilon_n=$(round(energies[n], sigdigits = 6))\\ \\mathrm{$(energy_unit)}"))
        ylims!(ax, 0, 1.35 * maximum(densities[n]))
        k < length(states) && hidexdecorations!(ax, ticks = false)
        push!(axes, ax)
    end
    axes[end].xlabel = axis_label("Position", "x", position_unit)
    linkxaxes!(axes...)
    xlims!(axes[end], xlimits...)
    rowsize!(fig.layout, 1, Relative(0.3))
    rowgap!(fig.layout, 3)
    Label(fig[2:rows, 0], axis_label("Probability density", "", "a.u."); rotation = π / 2,
          tellheight = false)
    return fig
end

"""
    plot_density_maps(x, y, densities, titles; columns = length(densities),
                      xname = "Coordinate", xsymbol = "x", yname = "Coordinate",
                      ysymbol = "y", unit = "a.u.", limits = nothing,
                      diagonal = false, panel_size = 120)

Two-dimensional densities as heat maps (`:inferno`), each with its own colour bar
and title. `densities` is a vector of matrices `ρ[i, j] = ρ(x[i], y[j])`. With
`diagonal = true` the line ``x = y`` is drawn (useful for two identical particles).
"""
function plot_density_maps(x, y, densities, titles; columns::Integer = length(densities),
                           xname = "Coordinate", xsymbol = "x", yname = "Coordinate",
                           ysymbol = "y", unit = "a.u.", limits = nothing,
                           diagonal::Bool = false, panel_size::Real = 120)
    rows = cld(length(densities), columns)
    fig = Figure(size = (columns * (panel_size + 115), rows * (panel_size + 70)))
    letters = 'a':'z'
    for (k, (ρ, title)) in enumerate(zip(densities, titles))
        row, col = fldmod1(k, columns)
        cell = fig[row, col] = GridLayout()
        ax = Axis(cell[1, 1]; title = title, width = panel_size, height = panel_size,
                  xlabel = axis_label(xname, xsymbol, unit),
                  ylabel = axis_label(yname, ysymbol, unit))
        hm = heatmap!(ax, x, y, ρ; colormap = DENSITY_COLORMAP)
        diagonal && lines!(ax, collect(extrema(x)), collect(extrema(x));
                           color = :gray60, linestyle = :dash, linewidth = 0.7)
        isnothing(limits) || limits!(ax, limits...)
        Colorbar(cell[1, 2], hm; height = panel_size,
                 label = axis_label("Density", "\\rho", "a.u."))
        colgap!(cell, 4)
        panel_label!(cell[1, 1, TopLeft()], letters[k])
    end
    return fig
end

"""
    plot_parameter_sweep(λ, energies; states = 1:size(energies, 1),
                         parameter_name = "Parameter", parameter_symbol = raw"\\lambda",
                         parameter_unit = "", width = 330, height = 250)

Eigenvalues as functions of a potential parameter, from a FEMTISE parameter sweep
(`energies[n, j]` is the `n`-th lowest eigenvalue for `λ[j]`). The ground state is
drawn as a thick coloured line and the excited states as thin gray lines, which
keeps the figure readable in black and white. Since the levels are sorted at every
`λ`, crossings between levels appear as sharp corners.
"""
function plot_parameter_sweep(λ, energies; states = 1:size(energies, 1),
                              parameter_name = "Parameter", parameter_symbol = "\\lambda",
                              parameter_unit = "", width = 330, height = 250)
    fig = Figure(size = (width, height))
    ax = Axis(fig[1, 1], xlabel = axis_label(parameter_name, parameter_symbol, parameter_unit),
              ylabel = axis_label("Eigenvalue", "\\varepsilon_n($parameter_symbol)", "a.u."))
    excited = nothing
    for n in states[2:end]
        excited = lines!(ax, λ, energies[n, :]; color = :gray45, linewidth = 0.7)
    end
    ground = lines!(ax, λ, energies[states[1], :]; color = series_style(1).color, linewidth = 1.8)
    xlims!(ax, extrema(λ)...)
    Legend(fig[2, 1], [ground, excited],
           [math_label("\\mathrm{ground\\ state}"),
            math_label("\\mathrm{excited\\ states}\\ (n=1,\\ldots,$(length(states) - 1))")];
           orientation = :horizontal, tellwidth = false)
    return fig
end

end # module
