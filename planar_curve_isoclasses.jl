using AbstractAlgebra
using Nemo
using Graphs
using Base.Threads


###################################################################
#
#      Initialization of inputs
#
###################################################################


if length(ARGS) != 3
    error("Please provide a prime p, a power n>0, and a degree d.")
end

p = parse(BigInt, ARGS[1])
n = parse(BigInt, ARGS[2])
d = parse(BigInt, ARGS[3])

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
n = reduce(n)
d = reduce(d)


#####################################################################
#
#      Constructing iterators for projective transformations
#
#####################################################################



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
