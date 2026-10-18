"""
    double_well(x, params)

Tilted double-well potential ``V(x) = λ (x² − 1)² + μ x`` (atomic units).

# Arguments
- `x`: mesh point (`x[1]` is the coordinate).
- `params::Tuple = (λ, μ)`: barrier height `λ` (the barrier at `x = 0` is `λ` above
  the minima at `x = ±1` when `μ = 0`) and tilt `μ` (energy difference ≈ `2μ`
  between the two wells).
"""
function double_well(x, params::Tuple)
    λ, μ = params
    return λ * (x[1]^2 - 1)^2 + μ * x[1]
end
