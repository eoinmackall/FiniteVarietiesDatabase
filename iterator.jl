using Nemo

# A helper function to check if a vector is the zero vector.
_is_zero(v::Vector) = all(iszero, v)

# A helper to check if v2 is a scalar multiple of v1.
# Assumes v1 is not the zero vector.
function _is_collinear(v1::Vector{T}, v2::Vector{T}) where {T<:FinFieldElem}
    # Find the first non-zero component of v1 to determine the potential scalar.
    pivot_idx = findfirst(!iszero, v1)
    isnothing(pivot_idx) && return _is_zero(v2) # Should not happen if v1 is non-zero

    scalar = v2[pivot_idx] / v1[pivot_idx]
    return v2 == scalar .* v1
end

"""
    GL3Iterator{T}
An iterator for the General Linear Group GL(3, F) over a finite field F.

# Fields
- `F::T`: The finite field, e.g., GF(q).
- `elements::Vector`: A cached vector of all elements in the field F.
- `q::Int`: The order of the field F.
"""
struct GL3Iterator{T<:Nemo.GaloisField}
    F::T
    elements::Vector
    q::Int

    function GL3Iterator(F::T) where {T<:Nemo.GaloisField}
        elements = collect(F)
        q = length(elements)
        new{T}(F, elements, q)
    end
end

# Define the type of elements the iterator will produce
Base.eltype(::GL3Iterator) = Matrix{elem_type(Nemo.GaloisField)}

# Define the total number of elements the iterator can produce
function Base.length(iter::GL3Iterator)
    q = iter.q
    return (q^3 - 1) * (q^3 - q) * (q^3 - q^2)
end


"""
Core iteration function. Given a state (the last generated matrix),
it finds and returns the next valid matrix in GL(3, F).
"""
function Base.iterate(iter::GL3Iterator, state::Tuple{Vector,Vector,Vector})
    F = iter.F
    all_vectors = Iterators.product(iter.elements, iter.elements, iter.elements)

    # Unpack the state (rows of the last matrix)
    r1_last, r2_last, r3_last = state

    # Start searching from the next potential third row
    for (i3, v3_tuple) in enumerate(all_vectors)
        # Fast-forward to the state we left off at
        if i3 <= (r3_last[1].d * iter.q^2 + r3_last[2].d * iter.q + r3_last[3].d)
            continue
        end
        r3 = [F(x.d) for x in v3_tuple] # Convert tuple to vector

        # Check if r3 is linearly independent of r1_last and r2_last
        # This is equivalent to det([r1; r2; r3]) != 0
        mat = matrix(F, [r1_last r2_last r3])'
        if !iszero(det(mat))
            return (mat, (r1_last, r2_last, r3))
        end
    end

    # If we exhausted all r3 for the given (r1, r2), find the next (r1, r2) pair
    # and then find the *first* valid r3 for it.

    # This outer loop structure is complex. A more robust way is to just search
    # from the day after the current state. The following code is a bit verbose
    # but logically follows the triple-nested loop.

    current_r1_idx = r1_last[1].d * iter.q^2 + r1_last[2].d * iter.q + r1_last[3].d
    current_r2_idx = r2_last[1].d * iter.q^2 + r2_last[2].d * iter.q + r2_last[3].d

    # Iterate through all possible rows
    for (i1, v1_tuple) in enumerate(all_vectors)
        if i1 < current_r1_idx
            continue
        end # Fast forward r1
        r1 = [F(x.d) for x in v1_tuple]
        if _is_zero(r1)
            continue
        end

        for (i2, v2_tuple) in enumerate(all_vectors)
            if i1 == current_r1_idx && i2 <= current_r2_idx
                continue
            end # Fast forward r2
            r2 = [F(x.d) for x in v2_tuple]
            if _is_collinear(r1, r2)
                continue
            end

            for (i3, v3_tuple) in enumerate(all_vectors)
                # No fast-forward for r3, we start from the beginning for a new (r1, r2)
                r3 = [F(x.d) for x in v3_tuple]
                mat = matrix(F, [r1 r2 r3])'
                if !iszero(det(mat))
                    return (mat, (r1, r2, r3))
                end
            end
        end
    end

    return nothing # Iteration is finished
end


"""
Initial iteration call. Finds the very first element of GL(3, F).
"""
function Base.iterate(iter::GL3Iterator)
    F = iter.F

    # Find the first valid matrix
    for r1_tuple in Iterators.product(iter.elements, iter.elements, iter.elements)
        r1 = [F(x.d) for x in r1_tuple]
        if _is_zero(r1)
            continue
        end

        for r2_tuple in Iterators.product(iter.elements, iter.elements, iter.elements)
            r2 = [F(x.d) for x in r2_tuple]
            if _is_collinear(r1, r2)
                continue
            end

            for r3_tuple in Iterators.product(iter.elements, iter.elements, iter.elements)
                r3 = [F(x.d) for x in r3_tuple]

                # Form matrix and check determinant
                mat = matrix(F, [r1 r2 r3])'
                if !iszero(det(mat))
                    return (mat, (r1, r2, r3))
                end
            end
        end
    end

    return nothing # Should not happen for q > 1
end

# --- Example Usage ---
println("Example for F = GF(2)")
F, _ = GaloisField(2, 1, "a")
gl3_f2 = GL3Iterator(F)

println("The order of GL(3, 2) is: $(length(gl3_f2))") # Should be 168
println("First 5 elements of GL(3, 2):")
for (i, M) in enumerate(gl3_f2)
    println(M)
    println("-"^10)
    if i == 5
        break
    end
end

println("\n" * "="^20 * "\n")

println("Example for F = GF(3)")
F3, _ = GaloisField(3, 1, "b")
gl3_f3 = GL3Iterator(F3)
println("The order of GL(3, 3) is: $(length(gl3_f3))") # Should be 11232
println("First 3 elements of GL(3, 3):")
# Note: Iterating the whole group would be slow!
counter = 0
for M in gl3_f3
    global counter
    counter += 1
    if counter <= 3
        println(M)
        println("-"^10)
    end
end
println("Total elements iterated: $counter")
