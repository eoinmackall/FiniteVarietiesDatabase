# This file contains code that will compute most, but potentially not all,
# projective equivalence classes when ran.

using Oscar
using Base.Threads
import Base.product

#####################################################################
#
#   Projective coefficient vector iterator 
#
#####################################################################

struct ProjectiveCoefficients
    F::FqField
    n::Int

    function ProjectiveCoefficients(F::FqField, n::Int)
        if n <= 1
            throw(ArgumentError("n should be the dimension of some projective space of positive dimension"))
        end
        new(F, n)
    end
end

function Base.length(coeffs::ProjectiveCoefficients)
    q = BigInt(size(coeffs.F))
    n = coeffs.n

    points = 0
    for i = 0:n
        points += q^i
    end

    return points
end

function Base.eltype(::Type{ProjectiveCoefficients})
    return Vector{Int}
end

function Base.iterate(coeffs::ProjectiveCoefficients)
    coeffs_init = zeros(Int, coeffs.n + 1)
    coeffs_init[1] = 1
    return (coeffs_init, (1, 0))
end

function Base.iterate(coeffs::ProjectiveCoefficients, state)

    k, counter = state
    q = BigInt(size(coeffs.F))
    n = coeffs.n

    if k > n
        return nothing
    elseif counter == q^(n + 1 - k) - 1
        coeffs_vec = zeros(Int, n + 1)
        coeffs_vec[k+1] = 1
        return (coeffs_vec, (k + 1, 0))
    else
        counter += 1
        coeffs_vec = zeros(Int, n + 1)
        coeffs_vec[k] = 1
        temp_counter = counter
        for i = 1:n+1-k
            temp_counter, coeffs_vec[k+i] = divrem(temp_counter, q)
        end
        return (coeffs_vec, (k, counter))
    end
end

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

#####################################################################
#
#   Projective equivalence classes
#
#####################################################################


"""
Produces a ?? containing a unique isomorphism class for each curve of degree d
over a finite field of p^n elements
"""
function projective_hypersurface_isomorphism_classes(p, r, n, d)

    #Order the vector space Sym^n(V) where V is k{x,y,z} by lexicographic order.
    #For each vector v in Sym^n, scale v so that first nonzero entry is 1
    #For each matrix M in the iterator, compute Mv and scale so that first nonzero entry is 1.
    #Keep a list/set of all outputs Mv, vary over all M.
    #vary over all v -- not in any previous set {Mw} for any w -- until #P(Sym^n(V))-elements.
    #
    #Can think about this as being symmetric pairs (v,w) -- there exists an M such that Mv=w then there is an N with Nw=v.
    #So can add symmetric pairs (v,w) to a graph in parallel, then count the number of connected components.
    #
    #Okay so here's how this is going to work.
    #1. Make a list of vertices, each representing one polynomial.
    #2. Make function that generates from a vertix v the graph Gv
    #3. In parallel, have multiple threads pick vertices v_1,...,v_n and make Gv_1,...,Gv_n
    #4. Maintain a list "Classes" that keeps track of Gv_i. Maintain a set hash(vertices(Gv_i)), and only add Gv_j to Classes if hash(vertices(Gv_j)) is not seen yet.

    #Set-up
    F = GF(p^r)
    M = matrix_space(F, n + 1, n + 1)
    R, x = polynomial_ring(F, ["x$i" for i = 0:n])
    S, x = grade(R)
    monomials_list = monomial_basis(S,d)

    #Create an iterator for representatives of all nonzero homogeneous degree d polynomials over F
    homogeneous_polynomials = (sum(c * m for (c, m) in zip(coeffs, monomials_list))
                               for coeffs in ProjectiveCoefficients(F,binomial(n+d,d)-1))

    #Create an itertator for representatives of PGL_n
    PGL = Iterators.filter(A -> is_PGL_rep(A), M)

    representatives = Set()
    isom_classes = Set()
    for f in homogeneous_polynomials
        equivalence_class = Set()
        for A in PGL
            y = A * x
            g = f(y...)
            c = leading_coefficient(term(g, 1))
            if c != 1
                g = g / c
            end
            push!(equivalence_class, g)
        end
        if !in(equivalence_class, isom_classes)
            push!(representatives, f)
            push!(isom_classes, equivalence_class)
        end
    end
    print(length(representatives))
    return representatives
end
