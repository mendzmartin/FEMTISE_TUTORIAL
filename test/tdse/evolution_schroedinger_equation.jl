"""
    initial_coefficients(initial_wave_function,eigen_states,differential_interior_FE_domain)

# Aim
    - This function computes the coefficients of the linear combination of the eigenstates of the Hamiltonian operator.
    The coefficients are computed by the inner product of the initial wave function and the eigenstates.

# Arguments
    - `initial_wave_function::CellField`: Initial wave function.
    - `eigen_states::Vector{CellField}`: Eigenstates of the Hamiltonian operator.
    - `differential_interior_FE_domain::Gridap.CellData.GenericMeasure`: Measure of the interior of the finite element domain.

# Returns
    - `coefficients::Vector{ComplexF64}`: Coefficients of the linear combination of the eigenstates.
"""
function initial_coefficients(
    initial_wave_function::CellField,
    eigen_states::Vector{CellField},
    differential_interior_FE_domain::Gridap.CellData.GenericMeasure)
    coefficients=zeros(ComplexF64,length(eigen_states))
    Threads.@threads for i in eachindex(eigen_states)
        coefficients[i]=sum(∫(conj(eigen_states[i])*initial_wave_function)*differential_interior_FE_domain)
    end
    return coefficients;
end

"""
    evolution_schrodinger(initial_wave_function,eigen_states,eigen_energies,trial_space,differential_interior_FE_domain,time)

# Aim
    - This function computes the evolution of the wave function in time.
    The wave function is given by the linear combination of the eigenstates of the Hamiltonian operator.
    The coefficients of the linear combination are computed by the inner product of the initial wave function and the eigenstates.
    The time evolution is computed by the Schrödinger equation.  

# Arguments
    - `initial_wave_function::CellField`: Initial wave function.
    - `eigen_states::Vector{CellField}`: Eigenstates of the Hamiltonian operator.
    - `eigen_energies::Vector{ComplexF64}`: Eigenvalues of the Hamiltonian operator.
    - `trial_space::FESpace`: Finite element space where the wave function is defined.
    - `differential_interior_FE_domain::Gridap.CellData.GenericMeasure`: Measure of the interior of the finite element domain.
    - `time::Vector{Float64}`: Time vector.

# Returns
    - `wave_function::Vector{CellField}`: Time evolution of the wave function.
"""
function evolution_schrodinger(
    initial_wave_function::CellField,
    eigen_states::Vector{CellField},
    eigen_energies::Vector{ComplexF64},
    trial_space::FESpace,
    differential_interior_FE_domain::Gridap.CellData.GenericMeasure,
    time::Vector{Float64})
    Planck_constant=1.0
    coeffvec=initial_coefficients(initial_wave_function,eigen_states,differential_interior_FE_domain)
    wave_function=Vector{CellField}(undef,length(time))
    for i in eachindex(time)
        factor = coeffvec[1]
        wave_function[i]=interpolate_everywhere((factor*eigen_states[1]),trial_space)
        for j in 2:length(eigen_energies)
            factor = coeffvec[j]*exp(-im*(1.0/Planck_constant)*time[i]*eigen_energies[j])
            wave_function[i]=interpolate_everywhere((wave_function[i]+factor*eigen_states[j]),trial_space)
        end
        wave_function[i]=interpolate_everywhere(wave_function[i],trial_space)
        norm_switch=false
        norm_switch ? (wave_function[i]=interpolate_everywhere((wave_function[i]*(1.0/norm_l2(wave_function[i],differential_interior_FE_domain))),trial_space)) : nothing
    end
    return wave_function;
end

function evolution_schrodinger(
    initial_wave_function::CellField,
    eigen_states::Vector{CellField},
    eigen_energies::Vector{ComplexF64},
    trial_space::FESpace,
    differential_interior_FE_domain::Gridap.CellData.GenericMeasure,
    time::Float64)
    Planck_constant=1.0
    coeffvec=initial_coefficients(initial_wave_function,eigen_states,differential_interior_FE_domain)
    # wave_function=interpolate_everywhere((coeffvec[1]*eigen_states[1]),trial_space)
    wave_function=0.0*eigen_states[1]
    for j in 1:length(eigen_energies)
        factor = coeffvec[j]*exp(-im*(1.0/Planck_constant)*time*eigen_energies[j])
        wave_function=interpolate_everywhere((wave_function+factor*eigen_states[j]),trial_space)
    end
    return wave_function;
end

"""
    first_euler_bilineal_forms(c,p,q,Δtime,differential_interior_FE_domain)

# Aim
    - This function computes the first bilinear form of the Schrödinger equation.

# Arguments
    - `c::Float64`: Constant.
    - `p::CellField`: Potential.
    - `q::CellField`: Potential.
    - `Δtime::Float64`: Time step.
    - `differential_interior_FE_domain::Gridap.CellData.GenericMeasure`: Measure of the interior of the finite element domain.

# Returns
    - `a::Function`: Bilinear form.
"""
function first_euler_bilineal_forms(c,p,q,
    Δtime::Float64,differential_interior_FE_domain::Gridap.CellData.GenericMeasure)
    a(trial_solution,test_solution)=((c*Δtime)*∫(p*(∇(trial_solution)⋅∇(test_solution)))differential_interior_FE_domain
    +(2.0)*∫(trial_solution*test_solution)differential_interior_FE_domain
    -(c*Δtime)*∫(q*(trial_solution*test_solution))differential_interior_FE_domain)
    return a;
end

"""
    second_euler_bilineal_forms(c,p,q,initial_trial_solution,Δtime,differential_interior_FE_domain)

# Aim
    - This function computes the second bilinear form of the Schrödinger equation.

# Arguments
    - `c::Float64`: Constant.
    - `p::CellField`: Potential.
    - `q::CellField`: Potential.
    - `initial_trial_solution::CellField`: Initial trial solution.
    - `Δtime::Float64`: Time step.
    - `differential_interior_FE_domain::Gridap.CellData.GenericMeasure`: Measure of the interior of the finite element domain.  

# Returns
    - `b::Function`: Bilinear form.
"""
function second_euler_bilineal_forms(c,p,q,initial_trial_solution,
    Δtime::Float64,differential_interior_FE_domain::Gridap.CellData.GenericMeasure)
    b(test_solution)=((c*Δtime)*∫(p*(∇(initial_trial_solution)⋅∇(test_solution)))differential_interior_FE_domain
    +(2.0)*∫(initial_trial_solution*test_solution)differential_interior_FE_domain
    +(c*Δtime)*∫(q*(initial_trial_solution*test_solution))differential_interior_FE_domain)
    return b;
end

# we need to use Gridap.FESpaces package
"""
    evolution_schroedinger_by_euler_01(initial_wave_function,trial_space,test_space,differential_interior_FE_domain,time,params)

# Aim
    - This function computes the evolution of the wave function in time.
    The wave function is given by the linear combination of the eigenstates of the Hamiltonian operator.
    The time evolution is computed by the Schrödinger equation using the Euler method.

# Arguments
    - `initial_wave_function::CellField`: Initial wave function.
    - `trial_space::FESpace`: Finite element space where the trial solution is defined.
    - `test_space::FESpace`: Finite element space where the test solution is defined.
    - `differential_interior_FE_domain::Gridap.CellData.GenericMeasure`: Measure of the interior of the finite element domain.
    - `time::Vector{Float64}`: Time vector.
    - `params::Tuple`: Parameters of the Schrödinger equation.  

# Returns
    - `wave_function::Vector{CellField}`: Time evolution of the wave function.
"""
function evolution_schroedinger_by_euler_01(
    initial_wave_function::CellField,
    trial_space::FESpace,test_space::FESpace,
    differential_interior_FE_domain::Gridap.CellData.GenericMeasure,
    time::Vector{Float64},params::Tuple)
    c,p,q=params
    Δtime=abs(time[2]-time[1])
    a=first_euler_bilineal_forms(c,p,q,Δtime,differential_interior_FE_domain);
    assem=SparseMatrixAssembler(trial_space,test_space)
    test_solution=get_fe_basis(test_space)
    trial_solution=get_trial_fe_basis(trial_space)
    mat_contribs=a(trial_solution,test_solution)
    data=collect_cell_matrix(trial_space,test_space,mat_contribs)
    A=assemble_matrix(assem,data)
    wave_function=Vector{CellField}(undef,length(time))
    initial_psi=initial_wave_function
    wave_function[1] = initial_psi
    Threads.@threads for i in 2:length(time)
        b=second_euler_bilineal_forms(c,p,q,initial_psi,Δtime,differential_interior_FE_domain)
        vec_contribs=b(test_solution)
        data=collect_cell_vector(test_space,vec_contribs)
        B=assemble_vector(assem,data)
        x=A\B
        trial_solution = FEFunction(trial_space,x)
        norm_switch=false
        if norm_switch
            wave_function[i]=interpolate_everywhere((trial_solution*(1.0/norm_l2(trial_solution,differential_interior_FE_domain))),trial_space)
        else
            wave_function[i]=interpolate_everywhere(trial_solution,trial_space)
        end
        initial_psi=wave_function[i]
    end
    return wave_function;
end

"""
    evolution_schroedinger_by_euler_02(initial_wave_function,trial_space,test_space,differential_interior_FE_domain,time,params)

# Aim
    - This function computes the evolution of the wave function in time.
    The wave function is given by the linear combination of the eigenstates of the Hamiltonian operator.
    The time evolution is computed by the Schrödinger equation using the Euler method.

# Arguments
    - `initial_wave_function::CellField`: Initial wave function.
    - `trial_space::FESpace`: Finite element space where the trial solution is defined.
    - `test_space::FESpace`: Finite element space where the test solution is defined.
    - `differential_interior_FE_domain::Gridap.CellData.GenericMeasure`: Measure of the interior of the finite element domain.
    - `time::Vector{Float64}`: Time vector.
    - `params::Tuple`: Parameters of the Schrödinger equation.

# Returns
    - `wave_function::Vector{CellField}`: Time evolution of the wave function.
"""
function evolution_schroedinger_by_euler_02(
    initial_wave_function::CellField,
    trial_space::FESpace,test_space::FESpace,
    differential_interior_FE_domain::Gridap.CellData.GenericMeasure,
    time::Vector{Float64},params::Tuple)
    c,p,q=params
    Δtime=abs(time[2]-time[1])
    a=first_euler_bilineal_forms(c,p,q,Δtime,differential_interior_FE_domain)
    ls = LUSolver()
    solver = LinearFESolver(ls)
    wave_function=Vector{CellField}(undef,length(time))
    new_initial_wave_function=initial_wave_function
    wave_function[1] = new_initial_wave_function
    for i in 2:length(time)
        b=second_euler_bilineal_forms(c,p,q,new_initial_wave_function,Δtime,differential_interior_FE_domain)
        op = AffineFEOperator(a,b,trial_space,test_space)
        trial_solution = Gridap.solve(solver,op)
        norm_switch=false
        if norm_switch
            wave_function[i]=interpolate_everywhere((trial_solution*(1.0/norm_l2(trial_solution,differential_interior_FE_domain))),trial_space)
        else
            wave_function[i]=interpolate_everywhere(trial_solution,trial_space)
        end
        new_initial_wave_function=wave_function[i]
    end
    return wave_function;
end