# Writes `pec.dat`, the potential energy curve of the tutorial, as a table of values.
#
# In practice such a table comes from an electronic-structure calculation or from
# experiment. Here it is generated from the Morse potential of H₂,
#     V(R) = D [1 - exp(-a (R - Rₑ))]²,
# so that the vibrational levels computed from the table can be compared with the
# exact Morse levels.
#
# Usage: julia make_potential_table.jl   (run from this folder)

const D = 0.1745      # dissociation energy (hartree)
const a = 1.0282      # range parameter (1/bohr)
const Rₑ = 1.4011     # equilibrium distance (bohr)

morse(R) = D * (1 - exp(-a * (R - Rₑ)))^2

open("pec.dat", "w") do io
    println(io, "# Morse potential energy curve of H2 (atomic units)")
    println(io, "# R (bohr)    V(R) (hartree)")
    for R in 0.0:0.02:8.0
        println(io, round(R, digits = 4), "    ", morse(R))
    end
end
println("wrote pec.dat")
