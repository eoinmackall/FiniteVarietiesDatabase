#####################################################################
#
#   Iterator for homogeneous monomials
#
#####################################################################

#This is faster than monomial_basis on homogeneous rings

function _homogeneous_monomial_exponents(n::Integer, d::Integer)
    HomogeneousExponents(n, d)
end

#Struct for generating exponent vectors in reverse lexicographical order
struct HomogeneousExponents
    num::Int
    deg::Int

    function HomogeneousExponents(n::Integer, d::Integer)
        if n < 0 || d < 0
            throw(ArgumentError("Both n and d must be non-negative integers"))
        end
        new(n, d)
    end
end

function Base.length(exps::HomogeneousExponents)
    if exps.num == 0
        return exps.deg == 0 ? 1 : 0
    end
    return binomial(exps.num + exps.deg - 1, exps.deg)
end

function Base.eltype(::Type{HomogeneousExponents})
    return Vector{Int}
end

function Base.iterate(exps::HomogeneousExponents)
    # Handle the edge case of n=0.
    # If n=0 and d=0, return an empty vector.
    # If n=0 and d>0, return nothing.
    if exps.num == 0
        return exps.deg == 0 ? (Int[0], nothing) : nothing
    end

    # The first exponent vector will be [d, 0, ..., 0].
    # This vector serves as the initial state for the next iteration.
    exps_init = zeros(Int, exps.num)
    exps_init[1] = exps.deg
    return (exps_init, exps_init)
end


"""
Generates vectors of length exps.num with nonnegative integer components
summing to exps.deg.
These vectors are generated in lexicographical order ``[d,0,...,0]`` -> ``[d-1,1,...,0]``.
"""
function Base.iterate(exps::HomogeneousExponents, state::Vector{Int})

    exps_vec = state

    #Stop if we reach [0,...,0,d]
    if exps.num > 0 && exps_vec[end] == exps.deg
        return nothing
    end

    # To find the next vector, we check from right to left
    i = 0
    for k in (exps.num-1):-1:1
        if exps_vec[k] > 0
            i = k
            break
        end
    end

    last = exps_vec[end]
    exps_vec[end] = 0

    # Decrement the pivot element.
    exps_vec[i] -= 1
    # Add the collected remainder (plus the one from the decrement) to the next element.
    exps_vec[i+1] = last + 1

    return (exps_vec, exps_vec)
end

function homogeneous_monomial_basis(S::MPolyRing, d::Int)

    x = gens(S)
    n = length(x)

    monomial_basis = []
    for v in HomogeneousExponents(n, d)
        f = S(0)
        f = setcoeff!(f, v, 1)
        push!(monomial_basis, f)
    end
    return monomial_basis
end

###################
#
#
#
###################

"""
Takes a nonzero vector with coefficients in the field F.
Checks if the leading nonzero coefficient is 1.
"""
function is_projective_rep(F, v)

    leading_coeff = findfirst(!iszero, v)

    return v[leading_coeff] == F(1)

end

function is_PGL_rep(A)

    first_column = @view A[:, 1]
    leading_coeff = findfirst(!iszero, first_column)

    return det(A) != 0 && leading_coeff == 1

end

function is_smooth_poly(f)

    mat = jacobian_matrix(parent(f), f)

end

