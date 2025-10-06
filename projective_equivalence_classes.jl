using Oscar
using Graphs
using Base.Threads
import Base.product

###################################################################
#
#      Initialization of inputs
#
###################################################################

# TODO: 
# Should add logic so that this doesn't crash if given integers instead of strings.
# Should also just make all of the integers BigInts from the start.
# Remove redundancy in calculating orbits (only need to check polynomials with leading coeff == 1)
# Make it parallel
# Maybe make it into a graph also


if length(ARGS) != 4
    error("Please provide a prime p, a power r>0, the dimension of ambient projective space n>1, and a degree d>0.")
end

p = parse(BigInt, ARGS[1])
r = parse(BigInt, ARGS[2])
n = parse(BigInt, ARGS[3])
d = parse(BigInt, ARGS[4])

#Add check that p is prime, n>0

function reduce(x::BigInt)

    if x < typemax(Int128)
        x = Int128(x)
    end

    if x < typemax(Int64)
        x = Int64(x)
    end

    if x < typemax(Int)
        x = Int(x)
    end

    return x
end

p = reduce(p)
r = reduce(r)
n = reduce(n)
d = reduce(d)


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

    return leading_coeff == F(1)

end

function is_PGL_rep(A)
    
    first_column = @view A[:, 1]
    leading_coeff = findfirst(!iszero, first_column)

    return det(M) !=0 && leading_coeff == 1

end

#####################################################################
#
#    Graph making function below
#
#####################################################################


#Create an iterator for GL_n(F)
#Create an iterator for PGL_n(F)
#Make function that takes matrix in GL_n(F) or PGL_n(F) and produces matrix in GL(Sym^m(F)) or PGL(Sym^m(F)) 

#Small question whether threads should access P^n(F) via an iterator, or whether this should be made and stored from the beginning

#Might just rely on the iterators constructed by collect(matrix_space(F,n,n)).
#Namely, call collect(matrix_space()) and det.collect(matrix_space()) then do findnext(x-> x!=0, det.collect(())) or something and then give matrix and index
#For the projective iterator, do the same, but skip the matrices where the first column has a non-1 leading term.

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
    S, x = polynomial_ring(F, ["x$i" for i = 0:n])

    #Create a list of all homogeneous monomials of degree d
    monomial_basis = homogeneous_monomial_basis(S, d)

    #Create an iterator for all nonzero homogeneous degree d polynomials over F
    coefficient_iterator = product(ntuple(_ -> F, length(monomial_basis))...)
    nonzero_coefficient_iterator = Iterators.drop(coefficient_iterator, 1)
    homogeneous_polynomials = (sum(c * m for (c, m) in zip(coeffs, monomial_basis))
                               for coeffs in nonzero_coefficient_iterator
                               if is_projective_rep(F, coeffs))
    

    #Create an itertator for representatives of PGL_n
    PGL = Iterators.filter(A->is_PGL_rep(A), M)

    isom_classes = Set()
    for p in homogeneous_polynomials
        equivalence_class = Set()
        for A in PGL
            y = A * x
            q = p(y...)
            push!(equivalence_class, q)
        end
    push!(isom_classes, equivalence_class)
    end

    print(length(isom_classes))
    return isom_classes
end

projective_hypersurface_isomorphism_classes(p, r, n, d)
