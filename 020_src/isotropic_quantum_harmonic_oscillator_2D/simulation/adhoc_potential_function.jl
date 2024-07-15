export qho_2d
"""
    qho_2d(x,params)

# Aim:
    This function is a simple implementation of the 2D quantum harmonic oscillator potential. 
    It is used to test the simulation of the isotropic quantum harmonic oscillator in 2D.

# Arguments
    x::Array{Float64,1} : The position of the particle in 2D.
    params::Tuple : A tuple containing the parameters of the potential. 
        params[1]::Float64 : The frequency of the oscillator.
        params[2]::Float64 : The position of the minimum of the potential in the x direction.
        params[3]::Float64 : The position of the minimum of the potential in the y direction.
"""
function qho_2d(x,params::Tuple)
    ω,x₁,y₁=params
    return 0.5*(ω*ω)*((x[1]-x₁)*(x[1]-x₁)+(x[2]-y₁)*(x[2]-y₁))
end
