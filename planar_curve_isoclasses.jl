using AbstractAlgebra
using Nemo
using Graphs
using Base.Threads


###################################################################
#
#      Initialization of inputs
#
###################################################################


#Should add logic so that this doesn't crash if given integers instead of strings.
#Should also just make all of the integers BigInts from the start.

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
#      Constructing iterators for projective transformations
#
#####################################################################



function _is_colinear(v::Vector{T}, w::Vector{T}) where {T<:FqFieldElem}

    if iszero(v)
        return true
    else
        leading_entry = findfirst(!iszero, v)
        scalar = w[leading_entry] / v[leading_entry]

        return w == scalar .* v
    end
end

"""

"""
struct GLiterator{T<:FqField}
    F::T
    n::Int

    function GLiterator(F::T, n::Int) where {T<:FqField}
        new{T}(F, n)
    end

end

GL(n::Int, F::FqField) = GLiterator(F,n)


struct PGLiterator{T<:FqField}
    F::T
    n::Int

    function PGLiterator(F::T, n::Int) where {T<:FqField}
        new{T}(F, n)
    end
end

PGL(n::Int, F::FqField) = PGLiterator(F,n)

"""

"""
function matrix_iterator



end



"""
For each matrix M in the above iterator, compute the symmetric representation of M, by acting on a basis.
"""
function symemtric_matrix_iterator

end

"""

"""
function projective_iterator

end




"""
Produces a ?? containing a unique isomorphism class for each curve of degree d
over a finite field of p^n elements
"""
function plane_curve_isomorphism_classes(p, n, d)

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
end
