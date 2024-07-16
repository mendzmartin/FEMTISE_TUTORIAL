distance(DOF1,DOF2)=abs(DOF2-DOF1)

effective_interaction(DOF) = exp(DOF^2)*erfc(DOF)
asymptotical_interaction(DOF) = 1.0/(sqrt(π)*DOF)

function avoided_divergence(DOF,b)
    DOF > 25.7 ? func=asymptotical_interaction(DOF+b) : func=effective_interaction(DOF+b)
    return func
end

function reduced_confinement_function(DOF1,DOF2,confinement_length,Yukawa_length,atomic_number)
    a=atomic_number*π*sqrt(π)*0.5*(1.0/confinement_length);
    b=confinement_length*0.5*(1.0/abs(Yukawa_length));
    coordinate=distance(DOF1,DOF2)*(1.0/confinement_length);
    return a*exp((-2.0*coordinate+b)*b)*avoided_divergence(coordinate,b)
end

export effective_Yukawa_potential_2e
"""
    effective_Yukawa_potential_2e(x,params)

# Aim:
    This function calculates the effective Yukawa potential between two electrons in a 1D system. 
    The potential is calculated using the reduced confinement function and the avoided divergence function.

# Arguments
    x::Array{Float64,1} : The position of the two electrons in 1D.
    params::Tuple : A tuple containing the parameters of the potential. 
        params[1]::Float64 : The confinement length of the potential.
        params[2]::Float64 : The Yukawa length of the potential.
        params[3]::Float64 : The position of the nuclei.
        params[4]::Float64 : The atomic number of the nuclei.
        params[5]::Float64 : The switch factor of the potential.

# Returns
    Float64 : The value of the effective Yukawa potential at the position x.
"""
function effective_Yukawa_potential_2e(x,params::Tuple)
    confinement_length,Yukawa_length,nuclei_coord,atomic_number,switch_factor=params;
    return (switch_factor*(reduced_confinement_function(x[1],x[2],confinement_length,Yukawa_length,1.0)
    -reduced_confinement_function(x[1],nuclei_coord,confinement_length,Yukawa_length,atomic_number)
    -reduced_confinement_function(nuclei_coord,x[2],confinement_length,Yukawa_length,atomic_number)))
end

export effective_Yukawa_potential_2e_without_repulsion
"""
    effective_Yukawa_potential_2e_without_repulsion(x,params)

# Aim:
    This function calculates the effective Yukawa potential between two electrons in a 1D system. 
    The potential is calculated using the reduced confinement function and the avoided divergence function.
    This function does not include the repulsion between the electrons.

# Arguments
    x::Array{Float64,1} : The position of the two electrons in 1D.
    params::Tuple : A tuple containing the parameters of the potential. 
        params[1]::Float64 : The confinement length of the potential.
        params[2]::Float64 : The Yukawa length of the potential.
        params[3]::Float64 : The position of the nuclei.
        params[4]::Float64 : The atomic number of the nuclei.
        params[5]::Float64 : The switch factor of the potential.

# Returns
    Float64 : The value of the effective Yukawa potential at the position x.
"""
function effective_Yukawa_potential_2e_without_repulsion(x,params::Tuple)
    confinement_length,Yukawa_length,nuclei_coord,atomic_number,switch_factor=params;
    return (switch_factor*(-reduced_confinement_function(x[1],nuclei_coord,confinement_length,Yukawa_length,atomic_number)
    -reduced_confinement_function(nuclei_coord,x[2],confinement_length,Yukawa_length,atomic_number)))
end
