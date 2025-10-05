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
# Edit Gemini's iterator so that I understand it
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

function compositions(n::Integer, d::Integer)
    IntegerCompositions(n, d)
end

# Internal struct representing the iterator.
# Users will typically create this via the `compositions` function.
struct IntegerCompositions
    n::Int
    d::Int

    function IntegerCompositions(n::Integer, d::Integer)
        if n < 0 || d < 0
            throw(ArgumentError("n and d must be non-negative integers"))
        end
        new(n, d)
    end
end

# --- Julia Iterator Interface Implementation ---

# Define the total number of compositions that will be generated.
function Base.length(it::IntegerCompositions)
    if it.n == 0
        return it.d == 0 ? 1 : 0
    end
    # This is the "stars and bars" formula
    return binomial(it.d + it.n - 1, it.n - 1)
end

# Define the type of element that the iterator yields.
function Base.eltype(::Type{IntegerCompositions})
    return Vector{Int}
end

# The initial call to the iterator.
function Base.iterate(it::IntegerCompositions)
    # Handle the edge case of n=0.
    # If n=0 and d=0, the single composition is an empty vector.
    # If n=0 and d>0, there are no compositions.
    if it.n == 0
        return it.d == 0 ? (Int[], nothing) : nothing
    end

    # The first composition in lexicographical order is [d, 0, ..., 0].
    # This vector also serves as the initial state for the next iteration.
    first_v = zeros(Int, it.n)
    first_v[1] = it.d
    return (first_v, first_v)
end

# Subsequent calls to the iterator.
function Base.iterate(it::IntegerCompositions, state::Vector{Int})
    # The `state` is the vector from the previous iteration.
    # We will mutate it to produce the next composition.
    v = state

    # The last composition is [0, 0, ..., d]. If we've reached it, stop.
    if it.n > 0 && v[end] == it.d
        return nothing
    end
    
    # To find the next lexicographical composition, we scan from right to left.
    # Find the rightmost element `v[i]` (not in the last position) that can be "moved".
    i = 0
    for k in (it.n - 1):-1:1
        if v[k] > 0
            i = k
            break
        end
    end

    # This condition should not be met due to the `v[end] == it.d` check, but acts as a safeguard.
    if i == 0
        return nothing
    end

    # Collect the sum of all elements to the right of v[i].
    # These values will be consolidated and moved to position i+1.
    remainder = 0
    for k in (i + 1):it.n
        remainder += v[k]
        v[k] = 0 # Reset these positions.
    end
    
    # Create the next composition by modifying the vector in place.
    # 1. Decrement the pivot element.
    v[i] -= 1
    # 2. Add the collected remainder (plus the one from the decrement) to the next element.
    v[i+1] = remainder + 1
    
    # Yield the mutated vector and pass it as the state for the next iteration.
    return (v, v)
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

    F=GF(p^r)
    M=matrix_space(F, n+1, n+1)
    S, x = polynomial_ring(F, ["x$i" for i=0:n])

    monomial_basis = []
    for v in compositions(n+1,d)
        f=S(0)
        f=setcoeff!(f, v, 1)
        push!(monomial_basis, f)
    end
    
    coefficient_iterator = product(ntuple(_ -> F, length(monomial_basis))...)
    homogeneous_polynomials = (sum(c*m for (c,m) in zip(coeffs, monomial_basis)) for coeffs in coefficient_iterator)

    isom_classes = Set()
    for p in homogeneous_polynomials
        equivalence_class=Set()
        for A in M
            if det(A) !=0
                first_column = @view A[:,1]
                leading_coeff = findfirst(!iszero, first_column)
                if leading_coeff == 1
                    y=A*x
                    q=p(y...)
                    push!(equivalence_class, q)
                end
            end
        end
        push!(isom_classes, equivalence_class)
    end
    print(length(isom_classes))
    return isom_classes
end

projective_hypersurface_isomorphism_classes(p,r,n,d)
