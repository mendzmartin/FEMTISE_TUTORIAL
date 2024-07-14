# Tutorial To Simulate Isotropic One Dimensional Harmonic Oscillator

### Create simulation directory
First of all we need to create a specific directory to save this specific simulation results 

```bash
@prompt:~$ mkdir ~/my_directory_path/QHO1D
@prompt:~$ cd ~/my_directory_path/QHO1D
```

### Input file

We need to create an input file to simulate using default solver function inside FEMTISE package

```bash
@prompt:~/my_directory_path/QHO1D$ vi input.dat
```
Inside `input.dat` we need to write the following.

```text
full_path_name              = ~/my_directory_path/QHO1D/name_output_file
dom_type                    = s
nev                         = 10
dimension                   = 1D
sigma                       = 0.0
adhoc_file_name             = ~/my_path_repo/FEMTISE.jl/src/functions/default_running_miscellaneous_functions
potential_function_name     = qho_1d
params_potential_types      = f f
params_potential            = 1.0 0.0
analysis_param              = false
output_format_type          = jld2 all
## ONLY FOR 1D EIGENPROBLEMS
L                           = 100
Δx                          = 0.1
## ONLY FOR 2D EIGENPROBLEMS
Lx                          = 
Ly                          = 
nx                          = 
ny                          = 
different_masses            = 
reduced_density             = 
```
Note that `~/my_path_repo/` is the directory path where we can find FEMTISE.jl package.

### Run script

Create Julia code as
```bash
@prompt:~/my_directory_path/QHO1D$ vi run.jl
```
Inside `run.jl` we need to write the following.

```julia
begin
    using Pkg
    Pkg.activate("../")
    develop_package = true; develop_package ? Pkg.develop(path="~/my_path_repo/FEMTISE.jl") : nothing
    Pkg.instantiate()
    install_pkg = true; install_pkg ? Pkg.add("Revise") : nothing
    using Revise, FEMTISE;
    run_default_eigen_problem(set_type_potential("~/my_directory_path/QHO1D/input.dat"))
end
```

After this we can run the simulation using Julia compiler (for example: using multithread running with four threads)

```bash
@prompt:~/my_directory_path/QHO1D$ julia -t 4 run.jl 
```

### Analysis

After running we obtain an output data file in jld2 format called `name_output_file_eigen_data.jld2`.

Then using Jupyter Notebook (by intermediate Visual Studio Code) we can analyse output file so:

```bash
@prompt:~/my_directory_path/QHO1D$ code QHO1D.ipynb
```
Inside `QHO1D.ipynb` we need to write the following:

#### Environment

Activate Julia environment

```julia
using Pkg
Pkg.activate("../")
Pkg.instantiate()
```
Is necessary to mark FEMTISE package as developed mackage using specific path repository:

```julia
develop_package = true; develop_package ? Pkg.develop(path="~/my_path_repo/FEMTISE.jl") : nothing
```

Now we install package (if is nesseary) and use specific packages to analyse output data:

```julia
install_pkg = true
if install_pkg
    Pkg.add("Revise")
    Pkg.add("Gridap")
    Pkg.add("Plots")
    Pkg.add("JLD2")
    Pkg.add("DelimitedFiles")
    Pkg.rm("SpecialFunctions")
    Pkg.add("PlotlyJS")
end
using Revise, FEMTISE;
using Gridap, Gridap.CellData;
using Plots, JLD2, DelimitedFiles;
```

#### Include Julia codes

Include specific Julia codes from FEMTISE package (where is defined QHO1D potential) and specific code with essential function to plot and read output data.

```julia
include("~/my_path_repo/FEMTISE.jl/test/test_1d_kronig_penney/miscellaneous_functions.jl")
include("./post_proccesing_result_data.jl")
```
#### Read output data

All the information that we need to specify is where we find input file then using specific functions we can collect output data

```julia
path_input_file_name="~/my_directory_path/QHO1D/input.dat"
simulation_info, output_data = collect_result_data(true,path_input_file_name)
```
#### Plotting figures

Now we can plot eigenenergies:
```julia
fig1 = plot_eigenvalues(simulation_info, output_data)
display(fig1)
```

and eigenfunctions:
```julia
range_to_show=range(1,step=1,length=3)
fig2 = plot_eigenstates(simulation_info, output_data,range_to_show)
display(fig2)
```