
function charge_input_data(full_path_input_data::String)
    return input_data(full_path_input_data)
end

function build_input_data(full_path_name::String)
    attributes=readdlm(full_path_name*"_eigen_problem_attributes.dat",String);
    dimension=attributes[2,8];
    if dimension == "1D"
        params=Params1D(dimension,parse(Float64,attributes[4,8]),attributes[5,5],parse(Float64,attributes[7,7]),
            parse(Float64,attributes[3,6]),parse(Float64,attributes[6,10]),"",tuple(nothing))
    elseif dimension == "2D"
        params=Params2D(dimension,parse(Float64,attributes[4,9]),parse(Float64,attributes[4,10]),attributes[5,5],parse(Int64,attributes[7,10]),
            parse(Int64,attributes[8,10]),parse(Float64,attributes[3,6]),parse(Float64,attributes[6,10]),
            "",tuple(nothing))
    end
    adhoc_file_name = ""; analysis_param = false; different_masses = false;
    return InputData(full_path_name,adhoc_file_name,params,analysis_param,different_masses);
end

struct DefaultBinEigenProblem{T}
    ϵ::Vector{ComplexF64}
    ϕ::Vector{CellField}
    r::T
end

struct DefaultBinEigenProblemReducedDensity{T}
    ϵ::Vector{ComplexF64}
    ϕ::Vector{CellField}
    r::T
    rhoDOF1::Vector{Float64}
    rhoDOF2::Vector{Float64}
end

struct DefaultJLD2EigenProblem{T1,T2}
    ϵ::Vector{ComplexF64}
    ϕ::Vector{CellField}
    r::T1
    pts::T2
end

struct DefaultJLD2EigenProblemReducedDensity{T1,T2}
    ϵ::Vector{ComplexF64}
    ϕ::Vector{CellField}
    r::T1
    pts::T2
    rhoDOF1::Vector{Float64}
    rhoDOF2::Vector{Float64}
end

struct DefaultJLD2AllEigenProblem{T1,T2,T3,T4,T5}
    ϵ::Vector{ComplexF64}
    ϕ::Vector{CellField}
    r::T1
    pts::T2
    Ω::T3
    dΩ::Gridap.CellData.GenericMeasure
    Γ::T4
    dΓ::Gridap.CellData.GenericMeasure
    USpace::FESpace
    VSpace::FESpace
    model::T5
end

struct DefaultJLD2AllEigenProblemReducedDensity{T1,T2,T3,T4,T5}
    ϵ::Vector{ComplexF64}
    ϕ::Vector{CellField}
    r::T1
    pts::T2
    Ω::T3
    dΩ::Gridap.CellData.GenericMeasure
    Γ::T4
    dΓ::Gridap.CellData.GenericMeasure
    USpace::FESpace
    VSpace::FESpace
    model::T5
    rhoDOF1::Vector{Float64}
    rhoDOF2::Vector{Float64}
end

struct AnalysisParam
    ϵ_matrix::Matrix{ComplexF64}
    λvector::Vector{ComplexF64}
end

function charge_results(id::InputData)
    if id.params.dimension == "1D"
        eigen_vectors_output= id.full_path_name*"_eigen_vectors.bin"
        eigen_values_output= id.full_path_name*"_eigen_values.bin"
        coordinates_output = id.full_path_name*"_coordinates.bin"
        ϵ=read_bin(eigen_values_output;matrix_data=false)
        ϕ=read_bin(eigen_vectors_output;matrix_data=true,c_dim=id.params.nev)
        r=read_bin(coordinates_output;matrix_data=false,c_dim=1);
    elseif id.params.dimension == "2D"
        eigen_vectors_output= id.full_path_name*"_eigen_vectors.bin"
        eigen_values_output= id.full_path_name*"_eigen_values.bin"
        coordinates_output = id.full_path_name*"_coordinates.bin"
        ϵ=read_bin(eigen_values_output;matrix_data=false)
        ϕ=read_bin(eigen_vectors_output;matrix_data=true,c_dim=id.params.nev)
        r=read_bin(coordinates_output;matrix_data=false,c_dim=1);
    end
    return DefaultBinEigenProblem(ϵ,ϕ,r)
end


function interpolation_eigenstates!(ϕ::Vector{CellField},USpace::FESpace)
    for i in eachindex(ϕ)
        ϕi = Interpolable(ϕ[i])
        ϕ[i] = interpolate_everywhere(ϕi,USpace)
    end
    return ϕ
end

function triangulation_repair(model,grid_type::String)
    FullDirichlet_values,FullDirichlet_tags=make_boundary_conditions(grid_type,"FullDirichlet",ComplexF64);
    Ω,dΩ,Γ,dΓ=measures(model,3,FullDirichlet_tags);
    reff = ReferenceFE(lagrangian,Float64,2);
    VSpace,USpace=fe_spaces(model,reff;BC_data=(FullDirichlet_values,FullDirichlet_tags),BC_type="Dirichlet")
    return Ω,dΩ,Γ,dΓ,VSpace,USpace
end

function charge_results(id::InputData1D)
    if id.analysis_param == false
        if id.output_format_type == ("bin","eigen") # use DelimitedFiles package
            eigen_vectors_output= id.full_path_name*"_eigen_vectors.bin"
            eigen_values_output= id.full_path_name*"_eigen_values.bin"
            coordinates_output = id.full_path_name*"_coordinates.bin"
            ϵ=read_bin(eigen_values_output;matrix_data=false)
            ϕ=read_bin(eigen_vectors_output;matrix_data=true,c_dim=id.params.nev)
            r=read_bin(coordinates_output;matrix_data=false,c_dim=1);
            result = DefaultBinEigenProblem(ϵ,ϕ,r)
        elseif id.output_format_type == ("jld2","eigen") # use JLD2 package
            ϵ = load("$(id.full_path_name)_eigen_data.jld2", "ϵ")
            ϕ = load("$(id.full_path_name)_eigen_data.jld2", "ϕ")
            r = load("$(id.full_path_name)_eigen_data.jld2", "r")
            pts = load("$(id.full_path_name)_eigen_data.jld2", "pts")
            result = DefaultJLD2EigenProblem(ϵ,ϕ,r,pts)
        elseif id.output_format_type == ("jld2","all")
            ϵ = load("$(id.full_path_name)_eigen_data.jld2", "ϵ")
            ϕ = load("$(id.full_path_name)_eigen_data.jld2", "ϕ")
            r = load("$(id.full_path_name)_eigen_data.jld2", "r")
            pts = load("$(id.full_path_name)_eigen_data.jld2", "pts")
            model = load("$(id.full_path_name)_eigen_data.jld2", "model")
            Ω,dΩ,Γ,dΓ,VSpace,USpace = triangulation_repair(model,"simple_line")
            ϕ = interpolation_eigenstates!(ϕ,USpace)
            result = DefaultJLD2AllEigenProblem(ϵ,ϕ,r,pts,Ω,dΩ,Γ,dΓ,USpace,VSpace,model)
        end
    elseif id.analysis_param ≠ false
        if id.output_format_type == ("jld2","eigen") # use JLD2 package
            ϵ_matrix = load("$(id.full_path_name)_eigen_data.jld2", "ϵ_matrix")
            λvector = load("$(id.full_path_name)_eigen_data.jld2", "λvector")
        elseif id.output_format_type == ("bin","eigen") # use DelimitedFiles package
            eigen_values_output= id.full_path_name*"_eigen_values.bin";
            param_values_output = id.full_path_name*"_param_values.bin";
            λvector=read_bin(param_values_output;matrix_data=false);
            ϵ_matrix=read_bin(eigen_values_output;matrix_data=true,c_dim=length(λvector));
        end
        result = AnalysisParam(ϵ_matrix,λvector)
    end
    return result
end

function charge_results(id::InputData2D)
    if id.analysis_param == false
        if id.output_format_type == ("bin","eigen") # use DelimitedFiles package
            eigen_vectors_output= id.full_path_name*"_eigen_vectors.bin"
            eigen_values_output= id.full_path_name*"_eigen_values.bin"
            coordinates_output = id.full_path_name*"_coordinates.bin"
            ϵ=read_bin(eigen_values_output;matrix_data=false)
            ϕ=read_bin(eigen_vectors_output;matrix_data=true,c_dim=id.params.nev)
            r=read_bin(coordinates_output;matrix_data=false,c_dim=1);
            if id.reduced_density
                reduced_density_DOF1_output = id.full_path_name*"_reduced_density_DOF1.bin";
                reduced_density_DOF2_output = id.full_path_name*"_reduced_density_DOF2.bin";
                rhoDOF1=read_bin(reduced_density_DOF1_output;matrix_data=true,c_dim=id.params.nev);
                rhoDOF2=read_bin(reduced_density_DOF2_output;matrix_data=true,c_dim=id.params.nev);
            end
            id.reduced_density ? (result = DefaultBinEigenProblemReducedDensity(ϵ,ϕ,r,rhoDOF1,rhoDOF2)) : 
                (result = DefaultBinEigenProblem(ϵ,ϕ,r))
        elseif id.output_format_type == ("jld2","eigen") # use JLD2 package
            ϵ = load("$(id.full_path_name)_eigen_data.jld2", "ϵ")
            ϕ = load("$(id.full_path_name)_eigen_data.jld2", "ϕ")
            r = load("$(id.full_path_name)_eigen_data.jld2", "r")
            pts = load("$(id.full_path_name)_eigen_data.jld2", "pts")
            if id.reduced_density
                rhoDOF1 = load("$(id.full_path_name)_eigen_data.jld2", "rhoDOF1")
                rhoDOF2 = load("$(id.full_path_name)_eigen_data.jld2", "rhoDOF2")
            end
            id.reduced_density ? (result = DefaultJLD2EigenProblemReducedDensity(ϵ,ϕ,r,pts,rhoDOF1,rhoDOF2)) : 
                (result = DefaultJLD2EigenProblem(ϵ,ϕ,r,pts))
        elseif id.output_format_type == ("jld2","all")
            ϵ = load("$(id.full_path_name)_eigen_data.jld2", "ϵ")
            ϕ = load("$(id.full_path_name)_eigen_data.jld2", "ϕ")
            r = load("$(id.full_path_name)_eigen_data.jld2", "r")
            pts = load("$(id.full_path_name)_eigen_data.jld2", "pts")
            model = load("$(id.full_path_name)_eigen_data.jld2", "model")
            Ω,dΩ,Γ,dΓ,VSpace,USpace = triangulation_repair(model,"Cartesian2D")
            ϕ = interpolation_eigenstates!(ϕ,USpace)
            if id.reduced_density
                rhoDOF1 = load("$(id.full_path_name)_eigen_data.jld2", "rhoDOF1")
                rhoDOF2 = load("$(id.full_path_name)_eigen_data.jld2", "rhoDOF2")
            end
            id.reduced_density ? (result = DefaultJLD2AllEigenProblemReducedDensity(ϵ,ϕ,r,pts,Ω,dΩ,Γ,dΓ,USpace,VSpace,model,rhoDOF1,rhoDOF2)) : 
                (result = DefaultJLD2AllEigenProblem(ϵ,ϕ,r,pts,Ω,dΩ,Γ,dΓ,USpace,VSpace,model))
        end
    else
        if id.output_format_type == ("jld2","eigen") # use JLD2 package
            ϵ_matrix = load("$(id.full_path_name)_eigen_data.jld2", "ϵ_matrix")
            λvector = load("$(id.full_path_name)_eigen_data.jld2", "λvector")
        elseif id.output_format_type == ("bin","eigen") # use DelimitedFiles package
            eigen_values_output= id.full_path_name*"_eigen_values.bin";
            param_values_output = id.full_path_name*"_param_values.bin";
            λvector=read_bin(param_values_output;matrix_data=false);
            ϵ_matrix=read_bin(eigen_values_output;matrix_data=true,c_dim=length(λvector));
        end
        result = AnalysisParam(ϵ_matrix,λvector)
    end
    return result
end

function collect_result_data(switch_input_file::Bool,full_path_name::String)
    switch_input_file ? (id = charge_input_data(full_path_name)) : (id = build_input_data(full_path_name))
    return id, charge_results(id)
end


function plot_eigenvalues(id,results,range_to_show::StepRange{Int, Int})
    plotlyjs()
    if id.analysis_param ≠ false
        figure = plot()
        for i in range_to_show
            figure = scatter!(real.(results.λvector),real(results.ϵ_matrix[i,:]), label="eigenenergy (ϵ$(i)) [au]",legend=:top)
            figure = plot!(real.(results.λvector),real(results.ϵ_matrix[i,:]),label="")
        end
        figure = plot!(xlabel="parameter λ",ylabel="eigenenergy (ϵₙ(λ)) [au]",ticks = :native)
    end
    return figure
end

function plot_eigenvalues(id,results)
    plotlyjs()
    if id.analysis_param == false
        figure = scatter(real(results.ϵ),xlabel="n",ylabel="eigenenergy (ϵₙ) [au]",legend=:bottomright)
    end
    return figure
end

# Pkg.add("PlotlyJS")

function plot_eigenstates(id,results,range_to_show::StepRange{Int, Int})
    if id.params.dimension == "1D"
        plotlyjs();
        figure = plot()
        rho=zeros(Float64,length(results.r))
        Threads.@threads for i in range_to_show
            if ((typeof(id) <: InputData) || (typeof(id) <: InputData1D && id.output_format_type == ("bin","eigen")))
                rho=real.(conj.((results.ϕ)[:,i]).*((results.ϕ)[:,i]))
            elseif (typeof(id) <: InputData1D && id.output_format_type in [("jld2","eigen"),("jld2","all")])
                rho = real.(conj.((results.ϕ[i]).(results.pts)).*((results.ϕ[i]).(results.pts)))
            end
            figure = plot!(results.r,rho,lw=2,lc=:"black",label="")
            figure = scatter!(results.r,rho,label="ϵ$(i)",lw=0.1)
        end
        figure = plot!(xlabel="x coordinate [au]",ylabel="probability ρ(x)",ticks = :native)
        return figure
    else
        println("PLOT ERROR.")
        println("You can not display 2D plotting")
    end
end

function density(x,y,phi_n)
    density=zeros(length(y),length(x))
    for i in eachindex(y)
        for j in eachindex(x)
            index=(j-1)*length(y)+i
            density[i,j]=real.(conj.(phi_n[index]).*phi_n[index])
        end
    end
    return density
end

struct EvolutionWavePacket2D{T1,T2}
    ϕ::Vector{CellField}
    r::T1
    pts::T2
end

function plot_eigenstates(id,results,index_nev::Int;mapcolor::Symbol=:rainbow1)
    if id.analysis_param == false
        if id.params.dimension == "1D"
            plotlyjs();
            if ((typeof(id) <: InputData) || (typeof(id) <: InputData1D && id.output_format_type == ("bin","eigen")))
                rho = real.(conj.(results.ϕ[:,index_nev]).*(results.ϕ[:,index_nev]))
            elseif (typeof(id) <: InputData1D && id.output_format_type in [("jld2","eigen"),("jld2","all")])
                rho = real.(conj.(results.ϕ[index_nev].(results.pts)).*(results.ϕ[index_nev].(results.pts)))
            end
            figure = plot(results.r,rho,lw=2,lc=:"black",label="")
            figure = scatter!(results.r,rho,label="ϵ$(index_nev)",lw=0.1)
            figure = plot!(xlabel="x coordinate [au]",ylabel="probability ρ(x)",ticks = :native)
        elseif id.params.dimension == "2D"
            gr();
            if ((typeof(id) <: InputData) || (typeof(id) <: InputData2D && id.output_format_type == ("bin","eigen")))
                if id.params.nx==id.params.ny
                    rho = density(results.r[:,1],results.r[:,2],results.ϕ[:,index_nev])
                    figure1 = contour(results.r[:,1],results.r[:,2],rho,levels=10, color=mapcolor, fill=true,lw=0)
                else
                    rho = density(results.r[1:id.params.nx,1],results.r[1:id.params.ny,2],results.ϕ[:,index_nev])
                    figure1 = contour(results.r[1:id.params.nx,1],results.r[1:id.params.ny,2],rho, levels=10, color=mapcolor, fill=true,lw=0)
                end
            elseif (typeof(id) <: InputData2D && id.output_format_type in [("jld2","eigen"),("jld2","all")])
                rho = density(results.r[1],results.r[2],results.ϕ[index_nev].(results.pts))
                figure1 = contour(results.r[1],results.r[2],rho,levels=10, color=mapcolor, fill=true,lw=0)
            end
            figure1 = plot(figure1,title="Probability \$\\rho(x,y)\$ for \$\\epsilon_{$(index_nev)}\$",
                xlabel="\$x\$ coordinate [au]", ylabel="\$y\$ coordinate [au]")
    
            if id.reduced_density
                plotlyjs()
                if ((typeof(id) <: InputData) || (typeof(id) <: InputData1D && id.output_format_type == ("bin","eigen")))
                    if id.params.nx==id.params.ny
                        figure2=plot(results.r[:,1],results.rhoDOF1[:,index_nev],label="ρ(x)")
                        figure2=plot!(results.r[:,2],results.rhoDOF2[:,index_nev],label="ρ(y)")
                    else
                        figure2=plot(results.r[1:id.params.nx,1],results.rhoDOF1[:,index_nev],label="ρ(x)")
                        figure2=plot!(results.r[1:id.params.ny,2],results.rhoDOF2[:,index_nev],label="ρ(y)")
                    end
                elseif (typeof(id) <: InputData2D && id.output_format_type in [("jld2","eigen"),("jld2","all")])
                    figure2 = plot(results.r[1],results.rhoDOF1[:,index_nev],label="ρ(x)")
                    figure2 = plot!(results.r[2],results.rhoDOF2[:,index_nev],label="ρ(y)")
                end
                figure2=plot!(title="Reduced probability for ϵ$(index_nev)")
                figure2=plot!(xlabel="x or y coordinate [au]",ylabel="",ticks = :native)

                figure = plot(figure1,figure2,layout=2)
            else
                figure = figure1
            end
        end
        return figure
    else
        println("PLOT ERROR.")
        println("You can not plot eigenstate with activated analysis params.")
    end
end