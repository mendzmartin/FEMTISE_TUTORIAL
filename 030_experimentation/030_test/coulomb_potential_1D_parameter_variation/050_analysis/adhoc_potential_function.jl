distance(coord_1,coord_2) = abs(coord_2-coord_1)


export coulomb_potential_interaction_1D
function coulomb_potential_interaction_1D(x,params::Tuple)
    charge_1,charge_2,minimun_interparticle_distance,reference_point=params
    (x[1] ≠ reference_point) ? interparticle_coord=distance(x[1],reference_point) : interparticle_coord=minimun_interparticle_distance
    return charge_1*charge_2*(1.0/interparticle_coord)
end

export coulomb_potential_interaction_2D
function coulomb_potential_interaction_2D(x,params::Tuple)
    charge_1,charge_2,minimun_interparticle_distance=params
    (x[1] ≠ x[2]) ? interparticle_coord=distance(x[1],x[2]) : interparticle_coord=minimun_interparticle_distance
    return charge_1*charge_2*(1.0/interparticle_coord)
end