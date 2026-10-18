"""
    TutorialSimulation

Helpers shared by the tutorial scripts: run the FEMTISE input-file driver only when
its results are missing, and sample finite element functions on regular grids.
"""
module TutorialSimulation

using FEMTISE
using Gridap

export simulate, sample_on_grid, sorted_spectrum

"""
    simulate(input_name; rerun = false)

Run the FEMTISE input-file driver for `input_name.dat` (path relative to the current
directory) and return `(id, results)` as given by `collect_result_data`.

The calculation is skipped if its output file already exists, unless `rerun = true`.
Folders named in `full_path_name` are created if needed.
"""
function simulate(input_name::AbstractString; rerun::Bool = false)
    settings = input_data(input_name)
    output_file = settings.full_path_name * "_eigen_data.jld2"
    if rerun || !isfile(output_file)
        mkpath(dirname(settings.full_path_name))
        simulation = set_type_potential(input_name)   # loads the potential file
        # The potential function was defined after this function was compiled, so
        # the driver is called in the latest "world age" (see the FEMTISE manual).
        Base.invokelatest(run_default_eigen_problem, simulation)
    else
        println("Using existing results in ", output_file, " (pass rerun = true to recompute).")
    end
    return collect_result_data(true, input_name)
end

"""
    sorted_spectrum(results)

Eigenvalues (real parts) in ascending order and the permutation that sorts them.
ARPACK returns the eigenpairs ordered by their distance to the shift `sigma`.
"""
function sorted_spectrum(results)
    ε = real.(results.ϵ)
    order = sortperm(ε)
    return ε[order], order
end

"""
    sample_on_grid(ψ, x)            -> Vector   (1D)
    sample_on_grid(ψ, x, y)         -> Matrix   (2D), M[i, j] = ψ(x[i], y[j])

Evaluate a finite element function on a regular grid.
"""
sample_on_grid(ψ, x::AbstractVector) = [ψ(Point(xi)) for xi in x]
sample_on_grid(ψ, x::AbstractVector, y::AbstractVector) = [ψ(Point(xi, yj)) for xi in x, yj in y]

end # module
