begin
    develop_package = false
    computer = "notebook"

    if computer == "notebook"
        path_repo = "/home/mendez/github_repositories/my_repositories/"
    elseif computer == "pcfamaf"
        path_repo = "/home/martin/github_repositories/my_repositories/"
    elseif computer == "ccad"
        path_repo = "/home/martinmendez/github_repositories/my_repositories/"
    end

    using Pkg
    Pkg.activate("../")

    develop_package ? Pkg.develop(path=path_repo*"FEMTISE.jl") : nothing

    Pkg.instantiate()
    # Pkg.status()
    # Pkg.precompile()
    # Pkg.build("GridapGmsh")
    
    # for update packages from old version to new version
    # Pkg.update("Gridap")
    # Pkg.update("JLD2")
    # Pkg.update("Revise")

    using Revise;
    using FEMTISE;

    run_default_eigen_problem(set_type_potential("/home/mendez/github_repositories/my_repositories/FEMTISE_TUTORIALS/030_experimentation/030_test/quantum_harmonic_oscillator_1D/020_input/input"))
end