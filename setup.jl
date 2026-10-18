# Prepare the environment of the tutorials (run once):
#
#     julia --project=. setup.jl
#
# Installs FEMTISE.jl v0.1.1 from GitHub (the version the tutorials were written
# for) and the plotting packages. Julia ≥ 1.11 reads the same information from the
# [sources] section of Project.toml; this script also covers Julia 1.10.

using Pkg
Pkg.activate(@__DIR__)
Pkg.add(url = "https://github.com/mendzmartin/FEMTISE.jl", rev = "v0.1.1")
Pkg.instantiate()
Pkg.precompile()
println("Environment ready. Run a tutorial with, e.g.,")
println("    julia --project=. src/isotropic_1d_harmonic_oscillator_potential/run_tutorial.jl")
