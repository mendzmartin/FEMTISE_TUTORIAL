# Run every tutorial and regenerate all figures.
#
#     julia --project=. make_figures.jl            # reuse existing results in */output
#     julia --project=. make_figures.jl --rerun    # recompute everything
#
# Each tutorial runs in its own Julia process, exactly as when it is run alone.

const TUTORIALS = [
    "isotropic_1d_harmonic_oscillator_potential",
    "isotropic_2d_harmonic_oscillator_potential",
    "symmetric_finite_1d_kronig_penney_potential",
    "coulomb_interaction_2d_potential",
]

for tutorial in TUTORIALS
    script = joinpath(@__DIR__, "src", tutorial, "run_tutorial.jl")
    println("\n=== ", tutorial)
    run(`$(Base.julia_cmd()) --project=$(@__DIR__) $script $ARGS`)
end
