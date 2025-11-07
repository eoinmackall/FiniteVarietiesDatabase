
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
#   Tools for iterating and creating projective representatives
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
enumerate_matrix(F::FqField, n::Int, i::Int)

Inputs: finite field F, integer n>0, integer i>=0.
Outputs: nxn matrix with coefficients in F gotten by writing
the integer i in base-q where q is the size of F.
"""
function enumerate_matrix(F::FqField, n::Int, i::Int)

    q = order(F)
    A = Matrix{Int}(undef, n, n)

    quo = i
    for j = 1:n
        for k = 1:n
            quo, A[j, k] = divrem(quo, q)
        end
    end
    A = matrix(F, A)
    return A
end

"""
is_projective_rep(F::FqField, v::Vector{FqFieldElem})

Takes a nonzero vector with coefficients in the field F.
Checks if the leading nonzero coefficient is 1.
"""
function is_projective_rep(F::FqField, v::Vector{FqFieldElem})

    leading_coeff = findfirst(!iszero, v)

    return v[leading_coeff] == F(1)

end

@doc raw"""
is_PGL_rep(A::FqMatrix)

Takes a square ``n\times n`` matrix and checks if both the determinant
is nonzero and the first nonzero entry of the first column, going from
top to bottom, is equal to 1.
"""
function is_PGL_rep(A::FqMatrix)

    first_column = @view A[:, 1]
    leading_coeff = findfirst(!iszero, first_column)

    return det(A) != 0 && first_column[leading_coeff] == base_ring(A)(1)

end

#####################################################################
#
#   Projective equivalence classes (orbit find)
#
#####################################################################

function _add_to_equiv_class(PGL_list, R::FqMPolyRing, f::FqMPolyRingElem)

    F = base_ring(R)
    x = gens(R)
    partial_equiv_class = Set{FqMPolyRingElem}()
    for A in PGL_list
        y = A * x
        g = f(y...)
        c = leading_coefficient(term(g, 1))
        if c != F(1)
            g = g / c
        end
        push!(partial_equiv_class, g)
    end
    return partial_equiv_class
end

function _add_to_equiv_class2(PGL_list, R::FqMPolyRing, f::FqMPolyRingElem, representatives::Set{FqMPolyRingElem}, seen::Atomic{Bool})

    F = base_ring(R)
    x = gens(R)
    for A in PGL_list
        if seen[]
            return seen[]
        end
        y = A * x
        g = f(y...)
        c = leading_coefficient(term(g, 1))
        if c != F(1)
            g = g / c
        end
        if in(g, representatives)
            seen[] = true
            return seen[]
        end
    end
    return seen[]
end


@doc raw"""
projective_hypersurface_equivalence_classes(p::Int, r::Int, n::Int, d::Int; cached::Bool=false)

Produces a set of polynomial representatives for projective equivalence
classes of hypersurfaces of degree d in ``\mathbb{P}^n`` over a finite field
of p^r elements.

Uses multithreaded union-find method. Can be interrupted for a partial set.

Setting cached=true will preallocate a vector of representatives for ``\mathbb{P}(V)``
where ``V`` is the vector space of homogeneous polynomials of degree d in n+1 variables.
This is a faster method, but can easily cause the system to run out of memory.
"""
function projective_hypersurface_equivalence_classes(p::Int, r::Int, n::Int, d::Int; cached::Bool=false)

    if cached
        return _projective_hypersurface_equivalence_classes1(p, r, n, d)
    else
        return _projective_hypersurface_equivalence_classes2(p, r, n, d)
    end
end

function _projective_hypersurface_equivalence_classes1(p::Int, r::Int, n::Int, d::Int)

    #Set-up
    q = p^r
    F = GF(q)
    M = matrix_space(F, n + 1, n + 1)
    R, x = polynomial_ring(F, ["x$i" for i = 0:n])
    monomial_basis = homogeneous_monomial_basis(R, d)


    println("Beginning memory allocation")
    #Create an iterator for representatives of nonzero homogeneous degree d polynomials over F
    homogeneous_polynomials_iterator = (sum(c * m for (c, m) in zip(coeffs, monomial_basis))
                                        for coeffs in ProjectiveCoefficients(F, binomial(n + d, d) - 1))
    homogeneous_polynomials = Set{FqMPolyRingElem}()
    sizehint!(homogeneous_polynomials, q^(binomial(n + d, d) - 1))
    union!(homogeneous_polynomials, homogeneous_polynomials_iterator)
    println("Polynomials built")

    #Create an itertator for representatives of PGL_(n+1)
    order = prod([q^(n + 1) - q^i for i = 0:n]) / (q - 1)
    PGL_iter = Iterators.filter(A -> is_PGL_rep(A), M)
    PGL = Vector{FqMatrix}()
    sizehint!(PGL, Int(order))
    for A in PGL_iter
        push!(PGL, A)
    end
    PGL_chunks = Iterators.partition(PGL, cld(Int(order), Threads.nthreads()))

    println("Starting collection process")
    representatives = Set{FqMPolyRingElem}()
    try
        while !isempty(homogeneous_polynomials)
            f = first(homogeneous_polynomials)
            push!(representatives, f)
            tasks = map(PGL_chunks) do chunk
                Threads.@spawn _add_to_equiv_class(chunk, R, f)
            end
            partial_equiv_classes = fetch.(tasks)
            setdiff!(homogeneous_polynomials, partial_equiv_classes...)
        end
    catch err
        if err isa InterruptException
            println("Terminating computation early")
        else
            rethrow(err)
        end
    finally
        return representatives
    end
end

function _projective_hypersurface_equivalence_classes2(p::Int, r::Int, n::Int, d::Int)

    #Set-up
    q = p^r
    F = GF(q)
    M = matrix_space(F, n + 1, n + 1)
    R, x = polynomial_ring(F, ["x$i" for i = 0:n])
    monomial_basis = homogeneous_monomial_basis(R, d)

    #Create an iterator for representatives of nonzero homogeneous degree d polynomials over F
    homogeneous_polynomials = (sum(c * m for (c, m) in zip(coeffs, monomial_basis))
                               for coeffs in ProjectiveCoefficients(F, binomial(n + d, d) - 1))

    println("Beginning memory allocation")
    #Create an itertator for representatives of PGL_(n+1)
    order = prod([q^(n + 1) - q^i for i = 0:n]) / (q - 1)
    PGL_iter = Iterators.filter(A -> is_PGL_rep(A), M)
    PGL = Vector{FqMatrix}()
    sizehint!(PGL, Int(order))
    for A in PGL_iter
        push!(PGL, A)
    end
    PGL_chunks = Iterators.partition(PGL, cld(Int(order), Threads.nthreads()))
    println("Starting collection process")
    representatives = Set{FqMPolyRingElem}()

    try
        for f in homogeneous_polynomials
            seen = Atomic{Bool}(false)
            @sync begin
                for chunk in PGL_chunks
                    @spawn _add_to_equiv_class2(chunk, R, f, representatives, seen)
                end
            end
            if !seen[]
                push!(representatives, f)
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
