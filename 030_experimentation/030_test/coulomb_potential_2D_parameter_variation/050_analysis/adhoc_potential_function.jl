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