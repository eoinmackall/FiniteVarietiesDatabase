# This file contains code that will compute most, but potentially not all,
# projective equivalence classes when ran.

using Oscar
using Base.Threads

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
Inputs: finite field F, integer n>0, integer i>=0.
Outputs: nxn matrix with coefficients in F gotten by writing
the i in base-q where q is the size of F.
"""
function enumerate_matrix(F,n,i)

    q=order(F)
    A=Matrix{Int}(undef,n,n)

    quo=i
    for j=1:n
        for k=1:n
            quo,A[j,k]=divrem(quo,q)
        end
    end
    return A
end

"""
Takes a nonzero vector with coefficients in the field F.
Checks if the leading nonzero coefficient is 1.
"""
function is_projective_rep(F, v)

    leading_coeff = findfirst(!iszero, v)

    return v[leading_coeff] == F(1)

end

@doc raw"""
Takes a square ``n\times n`` matrix and checks if both the determinant
is nonzero and the first nonzero entry of the first column, going from
top to bottom, is equal to 1.
"""
function is_PGL_rep(A)

    first_column = @view A[:, 1]
    leading_coeff = findfirst(!iszero, first_column)

    return det(A) != 0 && first_column[leading_coeff] == 1

end

#####################################################################
#
#   Projective equivalence classes
#
#####################################################################


@doc raw"""
Produces a set of polynomial representatives for projective equivalence
classes of hypersurfaces of degree d in ``\mathbb{P}^n`` over a finite field
of p^r elements.

Uses multithreaded union-find method. Can be interrupted for a partial set.
"""
function projective_hypersurface_equivalence_classes(p, r, n, d)

    #Set-up
    F = GF(p^r)
    M = matrix_space(F, n + 1, n + 1)
    R, x = polynomial_ring(F, ["x$i" for i = 0:n])
    monomial_basis = homogeneous_monomial_basis(R, d)

    #Create an iterator for representatives of nonzero homogeneous degree d polynomials over F
    homogeneous_polynomials = (sum(c * m for (c, m) in zip(coeffs, monomial_basis))
                               for coeffs in ProjectiveCoefficients(F, binomial(n + d, d) - 1))

    #Create an itertator for representatives of PGL_n
    PGL_iter = Iterators.filter(A -> is_PGL_rep(A), M)
    PGL = vec(collect(PGL_iter))


    representatives = Set()
    seen_polynomials = Set()
    equiv_class_lock = ReentrantLock()

    try
        for f in homogeneous_polynomials
            if in(f, seen_polynomials)
                continue
            else
                push!(representatives, f)
                equivalence_class = Set()
                @threads for A in PGL
                    y = A * x
                    g = f(y...)
                    c = leading_coefficient(term(g, 1))
                    if c != 1
                        g = g / c
                    end
                    lock(equiv_class_lock) do
                        push!(equivalence_class, g)
                    end
                end
                union!(seen_polynomials, equivalence_class)
            end
        end
    catch err
        if err isa InterruptException
            println("Terminating computation early")
        end
    finally
        return representatives
    end
end

#####################################################################
#
#   Burnside's Theorem for Counting Orbits
#
#####################################################################


"""
Input: an nxn-matrix A with n>0 and entries in a field F; an integer d>0.
Output: the matrix corresponding to the symmetric representation of A of
degree d.
"""
function symmetric_representation(A, d)

    F = base_ring(A)
    n = ncols(A) - 1
    R, x = polynomial_ring(F, ["x$i" for i = 0:n])
    monomial_basis = homogeneous_monomial_basis(R, d)
    monomials = vec(collect(monomial_basis))

    L = []
    for m in monomials
        y = A * x
        g = m(y...)

        column = zeros(F, binomial(n + d, d))
        for i in 1:length(g)

            t = term(g, i)
            c = leading_coefficient(t)
            new_term = t / c
            ind = findfirst(x -> (x == new_term), monomials)
            column[ind] = c
        end
        push!(L, column)
    end

    R = transpose(hcat(L...))
    return R
end

function _count_normals(p,r,n,chunk)
    
    F=GF(p^r)
    M=matrix_space(F,n+1,n+1)

    normal_forms = Dict()
    for i in chunk
        A=M(enumerate_matrix(F,n+1,i))
        if det(A) != 0
            N, _ = rational_canonical_form(A)
            normal_forms[N] = get(normal_forms,N,0)+1
        end
        
        if i % 10000 == 0
            GC.gc()
        end
    end
    GC.gc()
    return normal_forms
end

function finite_normal_forms(p, r, n)
     
    q=p^r    

    normal_forms = Dict()
    
    N=(BigInt(q)^(n+1)^2)
    enums_ind = 0:(N-1)

    chunks = Iterators.partition(enums_ind, cld(N, Threads.nthreads()))
    tasks = map(chunks) do chunk
        Threads.@spawn _count_normals(p,r,n,chunk)
    end
    chunk_dicts = fetch.(tasks)
    
    return merge(+,chunk_dicts)
end

function _conjugacy_number(list_of_matrices, q, d, n)
    
    F = base_ring(first(keys(normal_forms)))
    binom = binomial(n + d, d)
    M = matrix_space(F, binom, binom)
    I = identity_matrix(M) 
 
    j=0 
    for A in list_of_matrices
        RA = symmetric_representation(A, d)
        eigen = M(RA) - I
        mult = binomial(n + d, d) - rank(eigen)
        j += normal_forms[A] * q^mult
    end
    
end

function number_of_GL_orbits(normal_forms, d)

    F = base_ring(first(keys(normal_forms)))
    q = order(F)
    n = ncols(first(keys(normal_forms))) - 1
    
    keys_n = collect(keys(normal_forms))
    chunks = Iterators.partition(keys_n, cld(length(keys_n), Threads.nthreads()))
    tasks = map(chunks) do chunk
        Threads.@spawn _conjugacy_number(chunk, q, d, n)
    end
    chunk_sums = fetch.(tasks)
    j=sum(chunk_sums)

    ord = 1
    for i = 0:n
        ord *= (q^(n + 1) - q^i)
    end

    return j / ord - 1

end
