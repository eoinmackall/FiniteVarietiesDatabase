#####################################################################
#
#   Projective equivalence classes (filtration method)
#
#####################################################################

# In this file you'll find functions for constructing a set of representatives
# for all (nonzero) homogeneous polynomials of degree d in n+1 variables over,
# a finite field F up to a linear change in variables.
#
# If F=GF(q) for some prime power q, and if gcd(q-1,d)=1, then this is also
# a set of representatives for hypersurfaces of degree d in ``\mathbb{P}^n``.
#
# There is one primary function that accomplishes this below:
#   (1) projective_hypersurface_equivalence_classes_from_filtration
#   
# The general idea for this function is the filtration algorithm that's
# specified in (https://arxiv.org/abs/2306.09908). This algorithm considers
# an equivariant filtration on the space of polynomials, and lifts up a
# chain of quotients from smallest to largest. The logic that handles the
# construction of an equivariant chain can be found in ../misc/.

#####################################################################
#
#   Projective equivalence classes maker (helper functions)
#
#####################################################################

# Stores a vector (of type FqMatrix) as an integer. 
# Can recover the vector by looking
# at the base-q representation of the corresponding integer.
function pack_vector(v::FqMatrix, q::Int, buffer_elem::FqFieldElem)   
    
    h = zero(UInt64)
    n = length(v)
    if q == 2
        @inbounds for i in 1:n
            if !iszero(getindex!(buffer_elem, v, 1, i)) 
                h |= (one(UInt64) << (i-1))
            end
        end
    else
        p = one(UInt64)
        @inbounds for i in 1:n
            getindex!(buffer_elem, v, 1, i)
            h = h + (UInt64(lift(ZZ, buffer_elem)) * p) 
            p *= q
        end
    end
    return h
end


# Given a set "Space" and a point in the space, this function will calculate
# the orbit of a point from "Space" and return the generators for the stabilizer of point.
# Uses a collection of integer hashes to keep track of "Space" elements.
function orbit_stabilizer!(buffers, hash_space::Set{UInt64}, G_gens, point, action_dict)

    empty!(buffers.orbit_queue)
    empty!(buffers.transversal)
    empty!(buffers.stab_gens)

    push!(buffers.orbit_queue, point)
    
    F = base_ring(point)
    q = Int(size(F))
    buffer_elem = F()

    init_key = pack_vector(point.v, q, buffer_elem)
    buffers.transversal[init_key] = buffers.identity
    push!(hash_space, init_key)

    temp_coords = buffers.temp_coords
    M_parent = parent(point)

    while !isempty(buffers.orbit_queue)
        vec = pop!(buffers.orbit_queue)
        vec_coords = vec.v
        packed_v = pack_vector(vec_coords, q, buffer_elem)
        trans_vec = buffers.transversal[packed_v]

        for g_original in G_gens
            M_linear = action_dict[g_original]

            mul!(temp_coords, vec_coords, M_linear)
            packed_key = pack_vector(temp_coords, q, buffer_elem)

            if !haskey(buffers.transversal, packed_key)
                g_next_wrapped = M_parent(deepcopy(temp_coords))
                push!(buffers.orbit_queue, g_next_wrapped)
                buffers.transversal[packed_key] = trans_vec * g_original
                push!(hash_space, packed_key)
            else
                trans_next = buffers.transversal[packed_key]
                s = trans_vec * g_original * inv(trans_next)
                push!(buffers.stab_gens, s)
            end
        end
    end
    return collect(buffers.stab_gens)
end


function final_orbit_stabilizer!(buffers, hash_space::Set{UInt64}, G_gens, point, action_dict)

    empty!(buffers.orbit_queue)
 
    push!(buffers.orbit_queue, point)

    F = base_ring(point)
    q = Int(size(F))
    buffer_elem = F()

    init_key = pack_vector(point.v, q, buffer_elem)
    push!(hash_space, init_key)

    temp_coords = buffers.temp_coords
    M_parent = parent(point)

    while !isempty(buffers.orbit_queue)
        vec = pop!(buffers.orbit_queue)
        vec_coords = vec.v

        for g_original in G_gens
            M_linear = action_dict[g_original]
            mul!(temp_coords, vec_coords, M_linear)

            packed_key = pack_vector(temp_coords, q, buffer_elem)
            
            if !(packed_key in hash_space)
                push!(hash_space, packed_key)
                g_next_wrapped = M_parent(deepcopy(temp_coords))
                push!(buffers.orbit_queue, g_next_wrapped)
            end
        end
    end
    return
end


# Instead of passing matrices, and determining their action on the
# fly, we instead determine the linear action induced by these matrices
# on our corresponding filtration quotients.
function action_morphism(A, quotient, poly, inv_poly, x)

    W = quotient[1]
    pi = quotient[2]
    k = dim(W)

    # New coordinates
    y = A * x

    action_morphism = Vector{elem_type(W)}()
    sizehint!(action_morphism, k)
    for v in gens(W)
        v_vec = preimage(pi, v)
        f = poly(v_vec)
        g = evaluate(f,y)
        w_vec = inv_poly(g)
        push!(action_morphism, pi(w_vec))
    end
    return ModuleHomomorphism(W, W, action_morphism)
end


# For fast lookups
function action_dictionary_parallel(G_gens, quotient, poly, inv_poly, x)

    results = Vector{Tuple{FqMatrix, FqMatrix}}(undef, length(G_gens))

    Threads.@threads for i in eachindex(G_gens)
        A = G_gens[i]
        results[i] = (A, matrix(action_morphism(A, quotient, poly, inv_poly, x)))
    end
    return Dict(results)
end


# The main driver behind lifting orbit-reps up our chain.
function lift_orbit_reps(x, ImV_i, poly, q, inv_poly, f_i, orbits_stabilizers_dict)
    # W_i=V/V_i, f_i = V/V_{i+1}-> V/V_i, pi = V->W_i 
    # inc and inv_inc go between V and R_d = homogeneous_component of degree d
    # orbits = set in V/V_i

    W = domain(f_i)
    F = base_ring(W)
    dim_W = dim(W)
    Im_i = ImV_i[2].((ImV_i)[1])

    unique_matrices_vec = Vector{FqMatrix}()
    for gens in values(orbits_stabilizers_dict)
        for g in gens
            if !(g in unique_matrices_vec)
                push!(unique_matrices_vec, g)
            end
        end
    end
    
    field_q = Int(size(F))
    temp_n = ncols(unique_matrices_vec[end])

    if BigInt(field_q)^temp_n >= BigInt(2)^63
        @warn """Due to the way the code is written, the size of the spaces under
        consideration is too large to garauntee that results are correct.
        
        Proceed with caution.""" maxlog=1
    end

    identity_mat = identity_matrix(F, temp_n)

    global_action_cache = Dict{FqMatrix, FqMatrix}()
    
    cache_lock = ReentrantLock()
    
    Threads.@threads for A in unique_matrices_vec
        morph = action_morphism(A, (W,q), poly, inv_poly, x)
        M = matrix(morph)
        lock(cache_lock) do
            global_action_cache[A] = M
        end
    end

    orbits = collect(keys(orbits_stabilizers_dict))
    chunks = Iterators.partition(orbits, cld(length(orbits), Threads.nthreads()))

    tasks = map(chunks) do chunk
        Threads.@spawn begin
            partial_orbit_stab_dict = Dict{elem_type(W), Vector{FqMatrix}}()

            buffers = (
                transversal = Dict{UInt64, FqMatrix}(),
                orbit_queue = Vector{elem_type(W)}(),
                stab_gens = Set{FqMatrix}(),
                temp_coords = zero_matrix(F, 1, dim_W),
                identity = identity_mat,
            )
            sizehint!(buffers.transversal, 256)

            buffer_elem = F()

            for vec in chunk
                lift = preimage(f_i, vec) # preimage in V/V_{i+1}
                hash_space = Set{UInt64}()

                for v in Im_i
                    cand_vec = lift + v
                    cand_hash = pack_vector(cand_vec.v, field_q, buffer_elem)

                    if !(cand_hash in hash_space)
                        stab_vec_lift = orbit_stabilizer!(buffers, hash_space, orbits_stabilizers_dict[vec], cand_vec, global_action_cache)
                        partial_orbit_stab_dict[cand_vec] = stab_vec_lift
                    end
                end
            end
            return partial_orbit_stab_dict
        end
    end

    new_orbits_stab_dict = fetch(tasks[1])
    for i in 2:length(tasks)
        merge!(new_orbits_stab_dict, fetch(tasks[i]))
    end

    return new_orbits_stab_dict
end


# A final function that returns a set of polynomials, rather than vectors.
function final_lift_orbit_reps(x, ImV_i, poly, inv_poly, q, f_i, orbits_and_stabilizers, orbit_size)
    # W_i=V/V_i, f_i = V/V_{i+1}-> V/V_i, pi = V->W_i 
    # inc and inv_inc go between V and R_d = homogeneous_component of degree d
    # orbits = set in V/V_i

    W = domain(f_i)
    F = base_ring(W)
    dim_W = dim(W)
    final_orbits = Set{FqMPolyRingElem}()
    sizehint!(final_orbits, orbit_size)
    q_int = Int(size(F))

    Im_i = ImV_i[2].((ImV_i)[1])

    unique_matrices_vec = Vector{FqMatrix}()
    for gens in values(orbits_and_stabilizers)
        for g in gens
            if !(g in unique_matrices_vec)
                push!(unique_matrices_vec, g)
            end
        end
    end
    
    global_action_cache = Dict{FqMatrix, FqMatrix}()
    
    cache_lock = ReentrantLock()
    
    Threads.@threads for A in unique_matrices_vec
        morph = action_morphism(A, (W,q), poly, inv_poly, x)
        M = matrix(morph)
        lock(cache_lock) do
            global_action_cache[A] = M
        end
    end

    orbits = collect(keys(orbits_and_stabilizers))
    chunks = Iterators.partition(orbits, cld(length(orbits), 4*Threads.nthreads()))

    tasks = map(chunks) do chunk
        Threads.@spawn begin
            partial_orbit = Set{FqMPolyRingElem}()
            
            buffers = (
                orbit_queue = Vector{elem_type(W)}(),
                temp_coords = zero_matrix(F, 1, dim_W),
                gens = Vector{FqMatrix}()
            )

            buffer_elem = F()

            for vec in chunk
                lift = preimage(f_i, vec) # preimage in V/V_{i+1}
                hash_space = Set{UInt64}()

                for v in Im_i
                    cand_vec = lift + v
                    cand_hash = pack_vector(cand_vec.v, Int(q_int), buffer_elem)

                    if !(cand_hash in hash_space)
                        push!(partial_orbit, forget_grading(poly(cand_vec)))    
                        final_orbit_stabilizer!(buffers, hash_space, orbits_and_stabilizers[vec], cand_vec, global_action_cache)
                    end
                end
            end
            return partial_orbit
        end
    end

    for t in tasks
        union!(final_orbits, fetch(t))
    end
    
    return final_orbits
end


#####################################################################
#
#   Projective equivalence classes maker (main function)
#
#####################################################################


@doc raw"""

    projective_hypersurface_equivalence_classes_from_filtration(F::FqField, n::Int, d::Int;
        waring_samples=48, basis_samples=100, verbose=false, interactive=false)

Returns representatives for the set of equivalence classes of nonzero homogeneous polynomials of 
degree ``d`` in ``n+1`` variables over ``F``, up to linear changes of the variables.

If ``\gcd(\#F-1,d)=1``, then this gives a set of equivalence classes for hypersurfaces of degree ``d``
in ``\mathbb{P}^n``.

This function first produces a chain of invariant subspaces inside the appropriate space of polynomials,
and lifts representatives along a sequence of quotients. This part of the algorithm is probabilistic,
and can affect the duration of the algorithms run-time.

    `waring_samples` sets the number of "samples" of invariant subspaces. Increasing this may produce
    a chain of invariant subspaces of longer length.

    `basis_samples` sets the number of attempts to initially check whether a subspace is invariant.

    `verbose` provides more information about the stages of the functions runtime.

    `interactive` allows the user to choose a chain manually (from top to bottom, ordered from 1 to the 
    number of chains returned to the user).

# Example
```julia-repl
julia> projective_hypersurface_equivalence_classes_from_filtration(GF(2), 2, 3)
Set{FqMPolyRingElem} with 21 elements:
  x0^3 + x0^2*x1 + x0^2*x2 + x0*x1*x2 + x0*x2^2 + x1^3
  x0^3 + x0^2*x1 + x0^2*x2 + x0*x1^2 + x0*x2^2 + x1^3
  x0^3 + x0^2*x1 + x0^2*x2 + x0*x1*x2 + x1^3 + x1^2*x2
  x0^3 + x0^2*x2 + x0*x2^2 + x1^3 + x1^2*x2 + x1*x2^2
  ⋮
```
"""
function projective_hypersurface_equivalence_classes_from_filtration(F::FqField, n::Int, d::Int;
    waring_samples=240, basis_samples=100, verbose=false, interactive=false)
    
    if verbose == true
        println("Starting chain collection...")
    end

    #The Set-Up: Finding a suitable chain of subspaces.
    head = _chain_constructor(F, n, d; waring_samples, basis_samples, verbose)
    if interactive == true
        chains = collect_chains(head)
        println("Choose a chain to use as filtration:")
        filtration = chains[parse(Int, readline())]
    else
        rel_dim, position, filtration = _chain_finder(head)
        # Considering the filtration ordered as V=V_1 > V_2 > V_3 > ... > V_{end-2} > V_{end-1} > 0

        j = 0
        while length(filtration) == 2 || dim((filtration[end].object)[1]) != 0
            head = _chain_constructor(F, n, d; waring_samples, basis_samples, verbose)
            rel_dim, position, filtration = _chain_finder(head)
            if verbose == true
                println("Failed to find a good chain, trying again.")
            end
            if j == 10
                println("Failed to find invariant submodules. Is this module irreducible?")
                return
            end
            j += 1
        end
    end

    # Initializing some variables
    V, poly = head.object
    R = codomain(poly)
    x = gens(R)
    id = identity_map(V)
    inv_poly = inv(poly)

    if verbose == true
        if interactive == false
            println("Found chain at position #", position, " with maximal relative dimension = ", rel_dim)
        end
        println("Beginning orbit collection")
    end

    quotients = []
    for i = 1:length(filtration)-2
        push!(quotients, quo(V, (filtration[end-i].object)[1]))
        # Storing V/V_{end-i} at index i
        # So V/V_{end-1} to V/V_2
    end

    projection_maps = []
    for i = 1:length(quotients)-1
        Source_quo, quo_source = quotients[i]
        Target_quo, quo_target = quotients[i+1]

        gen_images = [quo_target(preimage(quo_source, v)) for v in gens(Source_quo)]
        pi = ModuleHomomorphism(Source_quo, Target_quo, gen_images)
        push!(projection_maps, pi)
        # Storing f_i : V/V_{end-i} ->> V/V_{end-i-1} at index i
    end

    images = []
    for i = 2:length(filtration)-2
        push!(images, image(compose((filtration[end-i].object)[2], quotients[i-1][2])))
        # Storing V_{end-i}/V_{end-i+1} at index i
    end

    # Starting the actual algorithm
    if verbose == true
        println("Starting stage #1 -- finding orbits in V/V_1")
    end

    G = GL(n + 1, F)
    GL_gens = [matrix(g) for g in gens(G)]

    # Find orbits in V/V_2, add to orbits.
    W_2, _ = quotients[end]
    dim_W2 = dim(W_2)
    buffers = (
        transversal = Dict{UInt64, FqMatrix}(),
        orbit_queue = Vector{elem_type(W_2)}(),
        stab_gens = Set{FqMatrix}(),
        temp_coords = zero_matrix(F,1,dim_W2),
        identity = identity_matrix(F, n+1),
    )
    sizehint!(buffers.transversal, 256)

    hash_space = Set{UInt64}()

    action_dict = action_dictionary_parallel(GL_gens, quotients[end], poly, inv_poly, x)
    orbits_and_stabilizers = Dict{elem_type(W_2), Vector{FqMatrix}}()
    q_int = Int(size(F))
    buffer_elem = F()

    for vec in W_2
        hash_vec = pack_vector(vec.v, q_int, buffer_elem)

        if !(hash_vec in hash_space)
            stab_gens = orbit_stabilizer!(buffers, hash_space, GL_gens, vec, action_dict)
            orbits_and_stabilizers[vec]=stab_gens
        end
    end

    if verbose == true
        println("Starting stage #2 -- lifting orbits along chain")
    end

    step = 1
    for i in (length(projection_maps)):-1:1

        if verbose == true
            println("Lifting step: #", step)
            step+=1
        end

        ImV_i = images[i]
        f_i = projection_maps[i]
        q = quotients[i][2]
        orbits_and_stabilizers = lift_orbit_reps(x, ImV_i, poly, q, inv_poly, f_i, orbits_and_stabilizers)
    end

    if verbose == true
        println("Starting stage #3 -- finding representatives in V")
    end

    orbit_size = Int(hypersurface_rep_size(normal_forms(F, n), d))
    pi = quotients[1][2]
    V_fin = (filtration[end-1].object)
    poly_reps = final_lift_orbit_reps(x, V_fin, poly, inv_poly, id, pi, orbits_and_stabilizers, orbit_size)
    delete!(poly_reps, forget_grading(R(0)))

    return poly_reps
end
