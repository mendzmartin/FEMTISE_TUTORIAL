export qho_1d
"""
    qho_1d(x,params)

# Aim: 
    This function is a simple implementation of the 1D quantum harmonic oscillator potential. 
    It is used to test the simulation of the isotropic quantum harmonic oscillator in 1D.

# Arguments
    x::Array{Float64,1} : The position of the particle in 1D.
    params::Tuple : A tuple containing the parameters of the potential. 
        params[1]::Float64 : The frequency of the oscillator.
        params[2]::Float64 : The position of the minimum of the potential.

# Returns
    Float64 : The value of the potential at the position x.
"""
function qho_1d(x,params::Tuple)
    ω,x₁=params
    return 0.5*(ω*ω)*((x[1]-x₁)*(x[1]-x₁))
end
