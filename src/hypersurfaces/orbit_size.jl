
#####################################################################
#
#   Burnside's Theorem for Counting Orbits
#
#####################################################################


"""
symmetric_representation(A::FqMatrix, d::Int)

Input: an nxn-matrix A with n>0 and entries in a field F; an integer d>0.
Output: the matrix corresponding to the symmetric representation of A of
degree d.
"""
function symmetric_representation(A::FqMatrix, d::Int)

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

    T = transpose(hcat(L...))
    return T
end

"""
normal_forms(F::FqField, n::Int)

Reutrns a dictionary with keys a representative of each conjugacy class
of GL(n+1,F) and with values the size of the corresponding class.
"""
function normal_forms(F::FqField, n::Int)

    q = Int(order(F))
    G = GAP.Globals.GL(n + 1, q)
    C = GAP.Globals.ConjugacyClasses(G)
    len = GAP.Globals.Length(C)

    N = Dict{FqMatrix,BigInt}()
    for i = 1:len
        a = GAP.Globals.Representative(C[i])
        num = BigInt(GAP.Globals.Size(C[i]))
        mat = matrix(F, a)
        N[mat] = num
    end

    return N
end

"""
Computes the contribution to the orbit size coming from an individual
conjugacy class, by way of Burnside's theorem.
"""
function _conjugacy_number(normal_forms::Dict{FqMatrix,BigInt}, q::ZZRingElem, d::Int, n::Int)

    F = base_ring(first(keys(normal_forms)))
    binom = binomial(n + d, d)
    M = matrix_space(F, binom, binom)
    I = identity_matrix(F, binom)

    j = 0
    for A in keys(normal_forms)
        RA = symmetric_representation(A, d)
        eigen = M(RA) - I
        mult = binomial(n + d, d) - rank(eigen)
        j += normal_forms[A] * q^mult
    end
    return j
end


@doc raw"""
orbit_size(normal_forms::Dict{FqMatrix,BigInt}, d::Int})

Inputs: a dictionary with keys a reprsentative of a given conjugacy class,
and with values the size of the conjugacy class, from GL(n+1,p^r), and an
integer d>0.

Outputs: the number of nonzero orbits of GL(n+1,p^r) on the dth symmetric
power of the canonical n+1 vector space that GL(n+1,p^r) acts on.

If every element of GF(p^r) is a dth power, e.g. if gcd(p^r-1,d)=1, then
this is the number of projective equivalence classes of hypersurfaces of
degree d in ``\mathbb{P}^n`` over GF(p^r).
"""
function orbit_size(normal_forms::Dict{FqMatrix,BigInt}, d::Int)

    F = base_ring(first(keys(normal_forms)))
    q = order(F)
    n = ncols(first(keys(normal_forms))) - 1

    chunks = Iterators.partition(normal_forms, cld(length(normal_forms), Threads.nthreads()))
    tasks = map(chunks) do chunk
        chunk_dict = Dict(chunk)
        Threads.@spawn _conjugacy_number(chunk_dict, q, d, n)
    end
    chunk_sums = fetch.(tasks)
    j = sum(chunk_sums)

    ord = 1
    for i = 0:n
        ord *= (q^(n + 1) - q^i)
    end

    return j / ord - 1

end
