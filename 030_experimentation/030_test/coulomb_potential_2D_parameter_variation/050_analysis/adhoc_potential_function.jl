function qho_2D(x,params::Tuple)
    frecuency,qho1_center,qho2_center=params
    return 0.5*(frecuency*frecuency)*((x[1]-qho1_center)*(x[1]-qho1_center)+(x[2]-qho2_center)*(x[2]-qho2_center))
end

distance(coord_1,coord_2) = abs(coord_2-coord_1)

function electronic_coulomb_potential_interaction_2D(x,params::Tuple)
    switch_interaction,minimun_interparticle_distance=params
    (x[1] ≠ x[2]) ? interparticle_coord=distance(x[1],x[2]) : interparticle_coord=minimun_interparticle_distance
    return switch_interaction*(1.0/interparticle_coord)
end

export effective_electronic_potential_2D
function effective_electronic_potential_2D(x,params::Tuple)
    switch_interaction,minimun_interparticle_distance=params
    qho_2D_params=tuple(1.0,0.0,0.0)
    return qho_2D(x,qho_2D_params)+electronic_coulomb_potential_interaction_2D(x,params)
end


# +++++++++++++++++++++++++++++++++++++++++++++++
measure_length(DOF1,DOF2)=abs(DOF2-DOF1)

effective_interaction(absx) = exp(absx^2)*erfc(absx)
asymptotical_interaction(absx) = (1/sqrt(π))*(1/absx)

function avoided_divergence(absx,b)
    absx > 25.7 ? func=asymptotical_interaction(absx+b) : func=effective_interaction(absx+b)
    return func
end

function reduced_confinement_function(DOF1,DOF2,confinement_length,Yukawa_length,atomic_number)
    a=atomic_number*π*sqrt(π)*0.5*(1.0/confinement_length);
    b=confinement_length*0.5*(1.0/abs(Yukawa_length));
    coordinate=measure_length(DOF1,DOF2)*(1.0/confinement_length);
    return a*exp((-2.0*coordinate+b)*b)*avoided_divergence(coordinate,b)
end

export effective_Yukawa_potential_2e
function effective_Yukawa_potential_2e(x,params::Tuple)
    confinement_length,Yukawa_length,nuclei_coord,nuclei_mass,atomic_number=params;
    # return (reduced_confinement_function(x[1],x[2],confinement_length,Yukawa_length,1.0)
    # -reduced_confinement_function(x[1],nuclei_coord,confinement_length,Yukawa_length,atomic_number)
    # -reduced_confinement_function(nuclei_coord,x[2],confinement_length,Yukawa_length,atomic_number))
    return (-reduced_confinement_function(x[1],nuclei_coord,confinement_length,Yukawa_length,atomic_number)
    -reduced_confinement_function(nuclei_coord,x[2],confinement_length,Yukawa_length,atomic_number))
end