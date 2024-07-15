# Pkg.add("PlotlyJS")
"""
    plot_eigenvalues(id,results,range_to_show;<keyword arguments>)

# Aim
- Plot the eigenvalues of the Hamiltonian operator.
The eigenvalues are plotted as a function of the parameter λ.
The eigenvalues are obtained from the diagonalization of the Hamiltonian operator.
The eigenvalues are plotted for the range of energy levels specified by the range_to_show variable.
The keyword arguments are used to set the title, xlabel, ylabel, and legend of the plot.

# Arguments
- `id`: InputData or InputData1D or InputData2D object.
- `results`: Results object.
- `range_to_show::StepRange{Int, Int}`: Range of energy levels to plot.
- `keyword arguments`:
    - `set_title::String`: Title of the plot.
    - `set_xlabel::String`: Label of the x-axis.
    - `set_ylabel::String`: Label of the y-axis.
    - `set_legend::Symbol`: Position of the legend.
"""
function plot_eigenvalues(id,results,range_to_show::StepRange{Int, Int};show_label=true)
    if id.analysis_param == false
        println("PLOT ERROR.")
        println("Check attributes, you are using the wrong function method. Analysis parameter is not activated.")
        figure = nothing
    else
        plotlyjs()
        figure = plot(xlabel="Parameter λ",ylabel="Eigen-energies (ϵn(λ) [au])",ticks = :native)
        for i in range_to_show
            if show_label
                figure = scatter!(real.(results.λvector),real(results.ϵ_matrix[i,:]), label="n=$(i)",legend=:top)
            else
                figure = scatter!(real.(results.λvector),real(results.ϵ_matrix[i,:]), label="")
            end
            figure = plot!(real.(results.λvector),real(results.ϵ_matrix[i,:]),label="")
        end
    end
    return figure
end

"""
    plot_eigenvalues(id,results;<keyword arguments>)

# Aim
- Plot the eigenvalues of the Hamiltonian operator.
The eigenvalues are obtained from the diagonalization of the Hamiltonian operator.
The keyword arguments are used to set the title, xlabel, ylabel, and legend of the plot.

# Arguments
- `id`: InputData or InputData1D or InputData2D object.
- `results`: Results object.
- `keyword arguments`:
    - `set_title::String`: Title of the plot.
    - `set_xlabel::String`: Label of the x-axis.
    - `set_ylabel::String`: Label of the y-axis.
    - `set_legend::Symbol`: Position of the legend.
"""
function plot_eigenvalues(id,results;
    set_title::String="",
    set_xlabel::String="Energy level (n)",set_ylabel::String="Eigen-energies (ϵn [au])",set_legend::Symbol=:bottomright)
    if id.analysis_param == false
        plotlyjs()
        figure = scatter(real(results.ϵ),title=set_title,xlabel=set_xlabel,ylabel=set_ylabel,legend=set_legend)
    else
        println("PLOT ERROR.")
        println("Check attributes, you are using the wrong function method. Analysis parameter is activated.")
        figure = nothing
    end
    return figure
end

"""
    plot_eigenstates(id,results,range_to_show;<keyword arguments>)

# Aim
- Plot the eigenstates of the Hamiltonian operator.
The eigenstates are obtained from the diagonalization of the Hamiltonian operator.
The eigenstates are plotted for the range of energy levels specified by the range_to_show variable.
The keyword arguments are used to set the title, xlabel, ylabel, and legend of the plot.

# Arguments
- `id`: InputData or InputData1D or InputData2D object.
- `results`: Results object.
- `range_to_show::StepRange{Int, Int}`: Range of energy levels to plot.
- `keyword arguments`:
    - `set_xlabel::String`: Label of the x-axis.
    - `set_ylabel::String`: Label of the y-axis.
"""
function plot_eigenstates(id,results,range_to_show::StepRange{Int, Int};
    set_xlabel::String="Coordinate (x [au])",set_ylabel::String="Probability density (ρ(x))")
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
            figure = plot!(results.r,rho,lw=2,lc=:black;label="")
            figure = scatter!(results.r,rho,label="n=$(i)",lw=0.1)
        end
        figure = plot!(xlabel=set_xlabel,ylabel=set_ylabel,ticks=:native)
    else
        println("PLOT ERROR.")
        println("Check attributes, you are using the wrong function method. 2D dimension problem is activated.")
        figure = nothing
    end
    return figure
end

"""
    density2D(x,y,phi_n)

# Aim
- Calculate the probability density of the eigenstates in 2D.

# Arguments
- `x::Array{Float64,1}`: Array of x-coordinates.
- `y::Array{Float64,1}`: Array of y-coordinates.
- `phi_n::Array{Complex{Float64},1}`: Array of eigenstates.
"""
function density2D(x,y,phi_n)
    density=zeros(length(y),length(x))
    for i in eachindex(y)
        for j in eachindex(x)
            index=(j-1)*length(y)+i
            density[i,j]=real.(conj.(phi_n[index]).*phi_n[index])
        end
    end
    return density
end

"""
    plot_eigenstates(id,results,index_nev;<keyword arguments>)

# Aim
- Plot the eigenstates of the Hamiltonian operator.
The eigenstates are obtained from the diagonalization of the Hamiltonian operator.
The eigenstates are plotted for the energy level specified by the index_nev variable.
The keyword arguments are used to set the color map of the plot (only for 2D plot).

# Arguments
- `id`: InputData or InputData1D or InputData2D object.
- `results`: Results object.
- `index_nev::Int`: Energy level to plot.
- `keyword arguments`:
    - `mapcolor::Symbol=:rainbow1`: Color map of the plot (only for 2D plot).
"""
function plot_eigenstates(id,results,index_nev::Int;mapcolor::Symbol=:rainbow1)
    if id.analysis_param == false
        if id.params.dimension == "1D"
            plotlyjs();
            if ((typeof(id) <: InputData) || (typeof(id) <: InputData1D && id.output_format_type == ("bin","eigen")))
                rho = real.(conj.(results.ϕ[:,index_nev]).*(results.ϕ[:,index_nev]))
            elseif (typeof(id) <: InputData1D && id.output_format_type in [("jld2","eigen"),("jld2","all")])
                rho = real.(conj.(results.ϕ[index_nev].(results.pts)).*(results.ϕ[index_nev].(results.pts)))
            end
            figure = plot(results.r,rho,lw=2,label="")
            figure = scatter!(results.r,rho,label="n=$(index_nev)",lw=0.1)
            figure = plot!(xlabel="Coordinate (x [au])",ylabel="Probability density (ρn(x))",ticks = :native)
        elseif id.params.dimension == "2D"
            gr();
            if ((typeof(id) <: InputData) || (typeof(id) <: InputData2D && id.output_format_type == ("bin","eigen")))
                if id.params.nx==id.params.ny
                    rho = density2D(results.r[:,1],results.r[:,2],results.ϕ[:,index_nev])
                    figure1 = contour(results.r[:,1],results.r[:,2],rho,levels=10, color=mapcolor, fill=true,lw=0)
                else
                    rho = density2D(results.r[1:id.params.nx,1],results.r[1:id.params.ny,2],results.ϕ[:,index_nev])
                    figure1 = contour(results.r[1:id.params.nx,1],results.r[1:id.params.ny,2],rho, levels=10, color=mapcolor, fill=true,lw=0)
                end
            elseif (typeof(id) <: InputData2D && id.output_format_type in [("jld2","eigen"),("jld2","all")])
                rho = density2D(results.r[1],results.r[2],results.ϕ[index_nev].(results.pts))
                figure1 = contour(results.r[1],results.r[2],rho,levels=10, color=mapcolor, fill=true,lw=0)
            end
            figure1 = plot(figure1,title="Probability density (ρ$(index_nev)(x))",xlabel="Coordinate (x [au])", ylabel="Coordinate (y [au])")
    
            if id.reduced_density
                plotlyjs()
                if ((typeof(id) <: InputData) || (typeof(id) <: InputData1D && id.output_format_type == ("bin","eigen")))
                    if id.params.nx==id.params.ny
                        figure2=plot(results.r[:,1],results.rhoDOF1[:,index_nev],label="ρn(x)")
                        figure2=plot!(results.r[:,2],results.rhoDOF2[:,index_nev],label="ρn(y)")
                    else
                        figure2=plot(results.r[1:id.params.nx,1],results.rhoDOF1[:,index_nev],label="ρn(x)")
                        figure2=plot!(results.r[1:id.params.ny,2],results.rhoDOF2[:,index_nev],label="ρn(y)")
                    end
                elseif (typeof(id) <: InputData2D && id.output_format_type in [("jld2","eigen"),("jld2","all")])
                    figure2 = plot(results.r[1],results.rhoDOF1[:,index_nev],label="ρn(x)")
                    figure2 = plot!(results.r[2],results.rhoDOF2[:,index_nev],label="ρn(y)")
                end
                figure2=plot!(title="Reduced probability density for level n=$(index_nev)")
                figure2=plot!(xlabel="Coordinate (x or y [au])",ylabel="",ticks = :native)

                figure = plot(figure1,figure2,layout=2)
            else
                figure = figure1
            end
        end
    else
        println("PLOT ERROR.")
        println("You can not plot eigenstate with activated analysis params.")
        figure = nothing
    end
    return figure
end