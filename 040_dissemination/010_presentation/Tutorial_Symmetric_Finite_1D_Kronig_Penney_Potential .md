# Tutorial To Simulate Symmetric Finite One Dimensional Kronig-Penney Potential

### Environment Activation

Activamos el environment particular para esta simulación, notemos que dentro de la función `activate` debemos colocar el path donde queremos localizar dicho environment.

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

#### Develop Package

En el caso de que sea necesario debemos agregar al environment el paquete FEMTISE con la ubicación local.

```julia
develop_package = false
path_repo="~/my_path_repo/"
develop_package ? Pkg.develop(path=path_repo*"FEMTISE.jl") : nothing
```

#### Adding Packages

Instalamos (si es necesario) y utilizamos paquetes específicos para la simulación.

```julia
install_pkg = false
if install_pkg
    Pkg.add("Revise")
    Pkg.add("Gridap")
    Pkg.add("Plots")
end
using Revise;
using FEMTISE;
using Gridap;
using Plots;
```

### Miscellaneous functions

Incluimos las funciones definidas dentro del paquete FEMTISE.

```julia
include(path_repo*"FEMTISE.jl/test/test_1d_kronig_penney/miscellaneous_functions.jl")
```

### Potential parameters

Definimos las propiedades del potencial

```julia
grid_size_length=276;
potential_depth=-0.5;
distance_between_wells=11;
well_width=1;
num_ions=23;
space_discretization=0.01;

unit_cell_potential=distance_between_wells+well_width;
```

#### Checking representation

Teniendo en cuenta el tamaño finito de la grilla de FE y las dimensiones de los pozos de potencial (ancho y separación) podremos chequear si la cantidad de sitios se pueden representar correctamente

```julia
quantity_check = num_ions*unit_cell_potential;
(grid_size_length ≥ quantity_check) ? println("The value of grid size length is ok (≥ $(quantity_check)).") : println("Increase grid size length, must be grid_size_length ≥ $(quantity_check).")
```

### Grid Building

Creamos la grilla de FE unidimensional

```julia
grid_type="simple_line";
params_model=("./","model1D",(-0.5*grid_size_length,0.5*grid_size_length),space_discretization);
model1D=make_model(grid_type,params_model);
rm(params_model[1]*params_model[2]*".msh")
```

El último paso podría omitirse si se quiere guardar la grilla en formato `.msh` para visualización externa.

#### Grid points

Construimos los vectores de puntos (puntos de evaluación de la grilla) y de coordenadas.

```julia
point_number=round(Int,abs(grid_size_length/space_discretization)+1)
space_coordinate,points=space_coord((-0.5*grid_size_length,0.5*grid_size_length),space_discretization,point_number-1;dimension="1D")
```

#### Plotting potential function

```julia
fig = plot(space_coordinate,symetric_kronig_penney(space_coordinate,num_ions,unit_cell_potential,well_width,potential_depth),label="")
fig = plot!(fig,xlabel="space coordinate (x [au])",ylabel="potential depth (v [au])")
display(fig)
```

### Boundary conditions

Definimos las condiciones de borde del sistema, en nuestro caso definimos condiciones de borde homogéneas en toda la frontera.

```julia
BC_type="FullDirichlet";
FullDirichlet_values,FullDirichlet_tags=make_boundary_conditions(grid_type,BC_type,ComplexF64);
```

### FE Domains

Construimos los dominios FE de la grilla: interiores y frontera. Además construimos los diferenciales de dichos dominios.

```julia
interior_FE_domain,differential_interior_FE_domain,boundary_FE_domain,differential_boundary_FE_domain = measures(model1D,3,FullDirichlet_tags)
```

### FE Reference

Definimos los polinomios de interpolación que se utilizarán y creamos los espacios Test y Trial asociados a las formulaciones débiles del problema.

```julia
reff = ReferenceFE(lagrangian,Float64,2)
TestSpace,TrialSpace = fe_spaces(model1D,reff,grid_type;BC_type=BC_type,TypeData=ComplexF64)
```

### Sturm-Liouville formulation

Definimos las funciones para utilizar la formulación de tipo Sturm-Liouville

```julia
p,q,r = kronig_penney_sturm_liouville((num_ions,unit_cell_potential,well_width,potential_depth))
```

### Eigen value problem

Resolvemos el problema de autovalores

```julia
eigen_energies,eigen_states = eigen_values_and_eigen_vectors(p,q,r,differential_interior_FE_domain,TrialSpace,TestSpace;
params=(10,1e-9,500,:none,potential_depth))
```

### Plotting results

```julia
fig=plot()
probability_densities=FEMTISE.density(eigen_states)
for i in 1:3#eachindex(ϕ)
    fig=plot!(space_coordinate,probability_densities[i].(points),
    label="probability density of energy=$(round(real(eigen_energies[i]),digits=4))",
    legend=:bottomleft)
end

fig=plot!(xlabel="space coordinate (x [au])",ylabel=" ")
fig=plot!(space_coordinate,0.1 .* symetric_kronig_penney(space_coordinate,num_ions,unit_cell_potential,well_width,potential_depth),
label="0.1*Kronig-Penney potential [au]")

display(fig)
```

We can save de figure using `savefig(fig,"010_example_kronig-penney.pdf")`.