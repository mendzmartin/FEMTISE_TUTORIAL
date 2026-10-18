# FEMTISE_TUTORIAL

Worked examples for [FEMTISE.jl](https://github.com/mendzmartin/FEMTISE.jl), a Julia
package that solves the time-independent Schrödinger equation with the finite
element method. Each tutorial runs a complete calculation with the FEMTISE
input-file driver, checks the result against exact or physical references, and
produces publication-quality figures.

| Tutorial | Physics | What it shows |
|---|---|---|
| [1D harmonic oscillator](docs/isotropic_1d_harmonic_oscillator_potential.md) | $\varepsilon_n=\omega(n+1/2)$ | input file, comparison with exact levels, upper bounds |
| [2D harmonic oscillator](docs/isotropic_2d_harmonic_oscillator_potential.md) | degenerate levels $\omega(n+1)$ | 2D meshes, degeneracies, density maps |
| [Finite Kronig–Penney chain](docs/symmetric_finite_1d_kronig_penney_potential.md) | band formation | discontinuous potentials, integer parameters |
| [Helium-like two-electron model](docs/coulomb_interaction_2d_potential.md) | electron correlation | two particles on a 2D mesh, parameter sweeps |
| [Vibrations from a tabulated potential](docs/tabulated_morse_vibrations.md) | H₂ Morse levels | potential data files, nuclear masses, interpolation accuracy |
| [Autocorrelation and spectrum](docs/autocorrelation_spectrum.md) | wave-packet revivals, Franck–Condon spectrum | $C(t)$, its Fourier transform $\sigma(\omega)$ |
| [Tilted double well](docs/double_well_parameter_map.md) | tunnelling splitting | two-parameter sweeps, energy maps |
| [Nucleus + electron + cavity photon](docs/pauli_fierz_cavity.md) | polaritons, Rabi splitting | Pauli–Fierz Hamiltonian, three coordinates with different masses |

The tutorials are written for **FEMTISE.jl v0.2.0** (the first four also run with v0.1.1).

## Quick start

```bash
git clone https://github.com/mendzmartin/FEMTISE_TUTORIAL.git
cd FEMTISE_TUTORIAL
julia --project=. setup.jl                    # once: installs FEMTISE v0.2.0 and CairoMakie
julia --project=. make_figures.jl             # runs every tutorial and regenerates all figures
```

or run a single tutorial:

```bash
julia --project=. src/isotropic_1d_harmonic_oscillator_potential/run_tutorial.jl
```

Results are stored in `src/<tutorial>/output/` (not versioned). A tutorial reuses
them when they exist; pass `--rerun` to recompute.

## Repository structure

```
FEMTISE_TUTORIAL/
├── Project.toml          # environment shared by all tutorials (FEMTISE v0.2.0, CairoMakie)
├── setup.jl              # installs the environment (also on Julia 1.10)
├── make_figures.jl       # runs all tutorials
├── tools/
│   ├── TutorialSimulation.jl   # simulate(), sorted_spectrum(), sample_on_grid()
│   └── TutorialFigures.jl      # plotting functions and figure style
├── src/<tutorial>/
│   ├── adhoc_potential_function.jl   # potential V(x; parameters)
│   ├── input*.dat                    # FEMTISE input files
│   ├── run_tutorial.jl               # simulation + checks + figures
│   ├── figures/                      # generated figures (PDF and PNG)
│   └── *.ipynb, Project.toml         # earlier interactive notebooks (optional)
└── docs/                 # one page per tutorial
```

## Figures

All figures are produced by the functions of `tools/TutorialFigures.jl` from the
results of the simulations. They follow the same rules:

- LaTeX typography; every axis shows a name, a symbol and the unit.
- Every data series is identified by colour (colour-blind-safe Okabe–Ito palette),
  line style *and* marker, so the figures can be printed in black and white.
- Density maps use the perceptually uniform `inferno` colour map, whose lightness
  increases monotonically, and always have a colour bar.
- Legends and labels never overlap data; axis ranges are fitted to the data.
- Each figure is exported as vector PDF and 300-dpi PNG.

The functions can be reused for your own FEMTISE calculations:

```julia
include("tools/TutorialFigures.jl"); using .TutorialFigures
set_publication_theme!()
fig = plot_spectrum(ε; exact = exact_levels)
save_figure(fig, "figures", "my_spectrum")
```

## Documentation and help

- FEMTISE.jl documentation: https://mendzmartin.github.io/FEMTISE.jl/stable/
- Questions and bug reports: https://github.com/mendzmartin/FEMTISE.jl/issues
