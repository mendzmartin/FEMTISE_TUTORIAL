export coulomb_potential_interaction_1d
function coulomb_potential_interaction_1D(x,params::Tuple)
    charge_1,charge_2,minimun_interparticle_distance=params
    return charge_1*charge_2*(1.0/max(x[1],minimun_interparticle_distance))
end