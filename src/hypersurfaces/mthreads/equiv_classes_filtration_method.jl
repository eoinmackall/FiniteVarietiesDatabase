#####################################################################
#
#   Projective equivalence classes (filtration method)
#
#####################################################################

# In this file you'll find functions for constructing a set of representatives
# for all (nonzero) homogeneous polynomials of degree d in n+1 variables over,
# a finite field F up to a linear change in variables.

# If F=GF(q) for some prime power q, and if gcd(q-1,d)=1, then this is also
# a set of representatives for hypersurfaces of degree d in ``\mathbb{P}^n``.

# There are two functions that accomplish this below:
#   (1) projective_hypersurface_equivalence_classes_from_filtration
#   (2) projective_hypersurface_equivalence_classes_from_filtration_GAP.
#
# The general idea for both functions is the same, and follows the algorithm
# specified in (https://arxiv.org/abs/2306.09908). However, their implementations
# are technically different, due to the thread unsafe nature of Oscar's GAP
# based stabilizer function.
#
# In general, the function (1) should be much faster. The function in (2) is
# included as an artifact of the development process (and just in case).

#####################################################################
#
#   Projective equivalence classes maker
#
#####################################################################

# This is the general algorithm used for finding the orbit of a vector
# (called point) and producing the stabilizing subgroup of that vector.
function orbit_stabilizer(G_gens, point, action_dict)

    orbit = Set([point])
    orbit_to_process = Set([point])
    transversal = Dict(point => identity_matrix(base_ring(G_gens[1]), ncols(G_gens[1])))
    stab_gens = Set{FqMatrix}()

    while !isempty(orbit_to_process)
        vec = pop!(orbit_to_process)
        for g in G_gens
            g_next = action_dict[g](vec)
            if !(g_next in orbit)
                push!(orbit, g_next)
                push!(orbit_to_process, g_next)
                transversal[g_next]=transversal[vec]*g
            else
                push!(stab_gens, transversal[vec]*g*inv(transversal[g_next]))
            end
        end
    end
    return (orbit, collect(stab_gens))
end


# Passing a set containing the point, we remove the orbit instead of return it.
# This function is not actually used in the final function, see the following function.
function orbit_stabilizer_pruning!(space::Set, G_gens, point, action_dict)

    orbit_to_process = Set([point])
    transversal = Dict(point => identity_matrix(base_ring(G_gens[1]), ncols(G_gens[1])))
    stab_gens = Set{FqMatrix}()

    while !isempty(orbit_to_process)
        vec = pop!(orbit_to_process)
        for g in G_gens
            g_next = action_dict[g](vec)
            if !haskey(transversal, g_next)
                delete!(space, g_next)
                push!(orbit_to_process, g_next)
                transversal[g_next]=transversal[vec]*g
            else
                push!(stab_gens, transversal[vec]*g*inv(transversal[g_next]))
            end
        end
    end
    return collect(stab_gens)
end

# A version of the above, but optimized to reduce allocations.
function orbit_stabilizer_pruning_optimized!(buffers, space::Set, G_gens, point, action_dict)

    empty!(buffers.orbit_queue)
    empty!(buffers.transversal)
    empty!(buffers.stab_gens)
    empty!(buffers.gen_pairs)

    for g in G_gens
        push!(buffers.gen_pairs, (g, action_dict[g]))
    end

    push!(buffers.orbit_queue, point)

    buffers.transversal[point.v] = buffers.identity
    
    temp_coords = buffers.temp_coords
    M_parent = parent(point)

    while !isempty(buffers.orbit_queue)
        vec = pop!(buffers.orbit_queue)
        vec_coords = vec.v
        trans_vec = buffers.transversal[vec_coords]

        for (g_original, M_linear) in buffers.gen_pairs
            mul!(temp_coords, vec_coords, M_linear)

            if !haskey(buffers.transversal, temp_coords)
                g_next_wrapped = M_parent(deepcopy(temp_coords))
                delete!(space, g_next_wrapped)
                push!(buffers.orbit_queue, g_next_wrapped)
                buffers.transversal[g_next_wrapped.v] = trans_vec * g_original
            else
                trans_next = buffers.transversal[temp_coords]
                s = trans_vec * g_original * inv(trans_next)
                push!(buffers.stab_gens, s)
            end
        end
    end
    return collect(buffers.stab_gens)
end

# Instead of passing matrices, and determining their action on the
# fly, we isntead determine the linear action induced by these matrices
# on our corresponding filtration quotients.
function action_morphism(A, quotient, poly, inv_poly, x)

    W = quotient[1]
    pi = quotient[2]
    k = dim(W)

    action_morphism = Vector{elem_type(W)}()
    sizehint!(action_morphism, k)
    for v in gens(W)
        v_vec = preimage(pi, v)
        f = poly(v_vec)
        y = A * x
        g = f(y...)
        w_vec = inv_poly(g)
        w = pi(w_vec)
        push!(action_morphism, w)
    end
    return ModuleHomomorphism(W, W, action_morphism)
end

# For fast look-ups.
function action_dictionary_parallel(G_gens, quotient, poly, inv_poly, x)

    results = Vector{Tuple{FqMatrix, FqMatrix}}(undef, length(G_gens))

    Threads.@threads for i in eachindex(G_gens)
        A = G_gens[i]
        results[i] = (A, matrix(action_morphism(A, quotient, poly, inv_poly, x)))
    end
    return Dict(results)
end

function action_dictionary(G_gens, quotient, poly, inv_poly, x)

    action_dict = Dict{eltype(G_gens), AbstractAlgebra.Generic.ModuleHomomorphism{FqFieldElem}}()

    for g in G_gens
        action_dict[g] = action_morphism(g, quotient, poly, inv_poly, x)
    end
    return action_dict
end

# The main driver behind lifting orbit-reps up our chain.
function lift_orbit_reps(x, ImV_i, poly, q, inv_poly, f_i, orbits_stabilizers_dict)
    # W_i=V/V_i, f_i = V/V_{i+1}-> V/V_i, pi = V->W_i 
    # inc and inv_inc go between V and R_d = homogeneous_component of degree d
    # orbits = set in V/V_i

    W=domain(f_i)
    F = base_ring(W)
    dim_W = dim(W)
    Im_i = ImV_i[2].((ImV_i)[1])


    unique_matrices = Set{FqMatrix}()
    for gens in values(orbits_stabilizers_dict)
        push!(unique_matrices, gens...)
    end
    unique_matrices_vec = collect(unique_matrices)
    
    identity_mat = identity_matrix(F, ncols(unique_matrices_vec[end]))

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
                transversal = Dict{FqMatrix, FqMatrix}(),
                orbit_queue = Vector{elem_type(W)}(),
                stab_gens = Set{FqMatrix}(),
                temp_coords = zero_matrix(F, 1, dim_W),
                identity = identity_mat,
                gen_pairs = Vector{Tuple{FqMatrix, FqMatrix}}()
            )
            sizehint!(buffers.gen_pairs, 16)
            sizehint!(buffers.transversal, 256)

            for vec in chunk
                lift = preimage(f_i, vec) # preimage in V/V_{i+1}
                coset = Set{elem_type(W)}(lift + v for v in Im_i)

                while !isempty(coset)
                    vec_lift=pop!(coset)
                    
                    stab_vec_lift = orbit_stabilizer_pruning_optimized!(buffers, coset, orbits_stabilizers_dict[vec], vec_lift, global_action_cache)
                    partial_orbit_stab_dict[vec_lift] = stab_vec_lift
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
function final_lift_orbit_reps(x, ImV_i, poly, inv_poly, q, f_i, orbits_and_stabilizers, size)
    # W_i=V/V_i, f_i = V/V_{i+1}-> V/V_i, pi = V->W_i 
    # inc and inv_inc go between V and R_d = homogeneous_component of degree d
    # orbits = set in V/V_i

    W=domain(f_i)
    F=base_ring(W)
    dim_W=dim(W)
    final_orbits = Set{FqMPolyRingElem}()
    sizehint!(final_orbits, size)

    Im_i = ImV_i[2].((ImV_i)[1])

    unique_matrices = Set{FqMatrix}()
    for gens in values(orbits_and_stabilizers)
        push!(unique_matrices, gens...)
    end
    unique_matrices_vec = collect(unique_matrices)
    
    identity_mat = identity_matrix(F, ncols(unique_matrices_vec[end]))

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
    chunks = Iterators.partition(orbits, cld(length(orbits), Threads.nthreads()))

    tasks = map(chunks) do chunk
        Threads.@spawn begin
            partial_orbit = Set{FqMPolyRingElem}()
            
            buffers = (
                transversal = Dict{FqMatrix, FqMatrix}(),
                orbit_queue = Vector{elem_type(W)}(),
                stab_gens = Set{FqMatrix}(),
                temp_coords = zero_matrix(F, 1, dim_W),
                identity = identity_mat,
                gen_pairs = Vector{Tuple{FqMatrix, FqMatrix}}()
            )
            sizehint!(buffers.gen_pairs, 16)
            sizehint!(buffers.transversal, 256)

            for vec in chunk
                lift = preimage(f_i, vec) # preimage in V/V_{i+1}
                coset = Set(lift + v for v in Im_i)

                while !isempty(coset)
                    vec_lift=pop!(coset)
                    push!(partial_orbit, forget_grading(poly(vec_lift)))
                    orbit_stabilizer_pruning_optimized!(buffers, coset, orbits_and_stabilizers[vec], vec_lift, global_action_cache)
                end
            end
            return partial_orbit
        end
    end

    union!(final_orbits, (fetch.(tasks))...)
    return final_orbits
end

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
    waring_samples=48, basis_samples=100, verbose=false, interactive=false)

    
    #The Set-Up: Finding a suitable chain of subspaces.
    head = _chain_constructor(F, n, d; waring_samples, basis_samples, verbose)
    if interactive == true
        chains = collect_chains(head)
        println("Choose a chain to use as filtration:")
        filtration = chains[parse(Int, readline())]
    else
        rel_dim, filtration = _chain_finder(head)
        # Considering the filtration ordered as V=V_1 > V_2 > V_3 > ... > V_{end-2} > V_{end-1} > 0

        j = 0
        while length(filtration) == 2 || dim((filtration[end].object)[1]) != 0
            head = _chain_constructor(F, n, d; waring_samples, basis_samples, verbose)
            rel_dim, filtration = _chain_finder(head)
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
            println("Found chain with maximal relative dimension = ", rel_dim)
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
        transversal = Dict{FqMatrix, FqMatrix}(),
        orbit_queue = Vector{elem_type(W_2)}(),
        stab_gens = Set{FqMatrix}(),
        temp_coords = zero_matrix(F,1,dim_W2),
        identity = identity_matrix(F, n+1),
        gen_pairs = Vector{Tuple{FqMatrix, FqMatrix}}()
    )
    sizehint!(buffers.gen_pairs, 16)
    sizehint!(buffers.transversal, 256)

    starting_vecs = Set(W_2)
    action_dict = action_dictionary_parallel(GL_gens, quotients[end], poly, inv_poly, x)
    orbits_and_stabilizers = Dict{elem_type(W_2), Vector{FqMatrix}}()

    while !isempty(starting_vecs)
        vec = pop!(starting_vecs)
        
        stab_gens = orbit_stabilizer_pruning_optimized!(buffers, starting_vecs, GL_gens, vec, action_dict)
        orbits_and_stabilizers[vec]=stab_gens
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

    size = Int(orbit_size(normal_forms(F, n), d))
    pi = quotients[1][2]
    V_fin = (filtration[end-1].object)
    poly_reps = final_lift_orbit_reps(x, V_fin, poly, inv_poly, id, pi, orbits_and_stabilizers, size)
    delete!(poly_reps, forget_grading(R(0)))

    return poly_reps
end

#####################################################################
#
#   Projective equivalence classes maker (using GAP/Oscar stabilizer)
#
#####################################################################


function action_for_stabilizer(A, v, p, poly, inv_poly, x)

    A = matrix(A)
    f = poly(preimage(p, v))
    y = A * x
    g = f(y...)
    g_vec = inv_poly(g)
    return p(g_vec)
end

function stabilizer_maker(G, orbits, action)
    stabilizers = []
    for vec in orbits
        stab, _ = stabilizer(G, vec, action)
        push!(stabilizers, stab)
    end
    return stabilizers
end

function stabilizer_maker_recursive(stabilizers, orbits, orbits_and_stabilizers, f_i, action)
    new_stabilizers = []
    for i in eachindex(orbits)
        vec = orbits[i]
        j = findfirst(pair -> pair[1] == f_i(vec), orbits_and_stabilizers)
        G = stabilizers[j]
        stab, _ = stabilizer(G, vec, action)
        push!(new_stabilizers, stab)
    end
    return new_stabilizers
end


function lift_orbit_representatives(x, ImV_i, p, poly, q, inv_poly, f_i, orbits_and_stabilizers)
    # W_i=V/V_i, f_i = V/V_{i+1}-> V/V_i, pi = V->W_i 
    # inc and inv_inc go between V and R_d = homogeneous_component of degree d
    # orbits = set in V/V_i

    new_orbits = []

    Im_i = ImV_i[2].((ImV_i)[1])

    orbits_and_stabs_chunks = Iterators.partition(orbits_and_stabilizers, cld(length(orbits_and_stabilizers), nthreads()))

    tasks = map(orbits_and_stabs_chunks) do chunk
        Threads.@spawn begin
            partial_orbit_reps = []
            for pair in chunk
                vec = pair[1]
                stabilizing_subgroup = pair[2]

                lift = preimage(f_i, vec) # preimage in V/V_{i+1}
                coset = Set(lift + v for v in Im_i)

                while !isempty(coset)
                    g = pop!(coset)
                    push!(partial_orbit_reps, g)

                    g_orbits = Set()
                    orbits_to_process = Set([g])
                    while !isempty(orbits_to_process)
                        g_vec = pop!(orbits_to_process)
                        for f in stabilizing_subgroup
                            orbit_vec = f(g_vec)
                            if !(orbit_vec in g_orbits)
                                push!(g_orbits, orbit_vec)
                                push!(orbits_to_process, orbit_vec)
                            end
                        end
                    end
                    setdiff!(coset, g_orbits)
                end
            end
            return partial_orbit_reps
        end
    end

    for partial_reps in fetch.(tasks)
        append!(new_orbits, partial_reps)
    end
    return new_orbits
end

function final_lift_orbit_representatives(ImV_i, poly, f_i, orbits_and_stabilizers, size)
    # W_i=V/V_i, f_i = V/V_{i+1}-> V/V_i, pi = V->W_i 
    # inc and inv_inc go between V and R_d = homogeneous_component of degree d
    # orbits = set in V/V_i

    final_orbits = Set{FqMPolyRingElem}()
    sizehint!(final_orbits, size)

    Im_i = ImV_i[2].((ImV_i)[1])

    orbits_and_stabs_chunks = Iterators.partition(orbits_and_stabilizers, cld(length(orbits_and_stabilizers), nthreads()))

    tasks = map(orbits_and_stabs_chunks) do chunk
        Threads.@spawn begin
            partial_orbit_reps = Set{FqMPolyRingElem}()
            for pair in chunk
                vec = pair[1]
                stabilizing_subgroup = pair[2]

                lift = preimage(f_i, vec) # preimage in V/V_{i+1}
                coset = Set(lift + v for v in Im_i)

                while !isempty(coset)
                    g = pop!(coset)
                    push!(partial_orbit_reps, forget_grading(poly(g)))

                    g_orbits = Set()
                    orbits_to_process = Set([g])
                    while !isempty(orbits_to_process)
                        g_vec = pop!(orbits_to_process)
                        for f in stabilizing_subgroup
                            orbit_vec = f(g_vec)
                            if !(orbit_vec in g_orbits)
                                push!(g_orbits, orbit_vec)
                                push!(orbits_to_process, orbit_vec)
                            end
                        end
                    end
                    setdiff!(coset, g_orbits)
                end
            end
            return partial_orbit_reps
        end
    end
    union!(final_orbits, (fetch.(tasks))...)
    return final_orbits
end

function projective_hypersurface_equivalence_classes_from_filtration_GAP(F, n, d; waring_samples=48, basis_samples=100, verbose=false, interactive=false)


    head = _chain_constructor(F, n, d; waring_samples, basis_samples, verbose)
    if interactive == true
        chains = collect_chains(head)
        println("Choose a chain to use as filtration:")
        filtration = chains[parse(Int, readline())]
    else
        rel_dim, filtration = _chain_finder(head)
        # Considering the filtration ordered as V=V_1 > V_2 > V_3 > ... > V_{end-2} > V_{end-1} > 0

        j = 0
        while length(filtration) == 2 || dim((filtration[end].object)[1]) != 0
            head = _chain_constructor(F, n, d; waring_samples, basis_samples, verbose)
            rel_dim, filtration = _chain_finder(head)
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

    V, poly = head.object
    R = codomain(poly)
    x = gens(R)
    id = identity_map(V)
    inv_poly = inv(poly)

    if verbose == true
        if interactive == false
            println("Found chain with maximal relative dimension = ", rel_dim)
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


    if verbose == true
        println("Starting stage #1 -- finding orbits in V/V_1")
    end

    G = GL(n + 1, F)
    GL_gens = [action_morphism(matrix(g), quotients[end], poly, inv_poly, x) for g in gens(G)]
    orbits = []

    # Find orbits in V/V_2, add to orbits.
    W_2, f_2 = quotients[end]
    starting_vecs = Set(W_2)

    while !isempty(starting_vecs)
        vec = pop!(starting_vecs)
        push!(orbits, vec)

        orbit_of_vec = Set([vec])
        orbit_to_process = Set([vec])

        while !isempty(orbit_to_process)

            v = pop!(orbit_to_process)
            for f in GL_gens
                h_vec = f(v)
                if !(h_vec in orbit_of_vec)
                    push!(orbit_of_vec, h_vec)
                    push!(orbit_to_process, h_vec)
                end
            end
        end
        setdiff!(starting_vecs, orbit_of_vec)
    end

    stabilizers = stabilizer_maker(G, orbits, (v, A) -> action_for_stabilizer(A, v, f_2, poly, inv_poly, x))
    stabilizers_gens = Vector{Vector{FqMatrix}}(undef, length(stabilizers))
    for i in eachindex(stabilizers)
        stabilizers_gens[i] = [matrix(g) for g in gens(stabilizers[i])]
    end

    stabilizer_maps = Vector{Vector{AbstractAlgebra.Generic.ModuleHomomorphism{FqFieldElem}}}(undef, length(stabilizers))
    @threads for i in eachindex(stabilizers_gens)
        num_gens = length(stabilizers_gens[i])
        gens_vec = Vector{AbstractAlgebra.Generic.ModuleHomomorphism{FqFieldElem}}(undef, num_gens)

        for j in 1:num_gens
            gens_vec[j] = action_morphism(stabilizers_gens[i][j], quotients[end-1], poly, inv_poly, x)
        end
        stabilizer_maps[i] = gens_vec
    end
    orbits_and_stabilizers = collect(zip(orbits, stabilizer_maps))


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
        p = quotients[i+1][2]
        q = quotients[i][2]
        orbits = lift_orbit_representatives(x, ImV_i, p, poly, q, inv_poly, f_i, orbits_and_stabilizers)

        if verbose == true
            println("Computing stabilizers")
        end

        stabilizers = stabilizer_maker_recursive(stabilizers, orbits, orbits_and_stabilizers, f_i, (v, A) -> action_for_stabilizer(A, v, q, poly, inv_poly, x))
        stabilizer_maps = Vector{Vector{AbstractAlgebra.Generic.ModuleHomomorphism{FqFieldElem}}}(undef, length(stabilizers))
        stabilizers_gens = Vector{Vector{FqMatrix}}(undef, length(stabilizers))
        for j in eachindex(stabilizers)
            stabilizers_gens[j] = [matrix(g) for g in gens(stabilizers[j])]
        end

        if verbose == true
            println("Setting-up for next lift")
        end

        if i > 1
            @threads for j in eachindex(stabilizers_gens)
                num_gens = length(stabilizers_gens[j])
                gens_vec = Vector{AbstractAlgebra.Generic.ModuleHomomorphism{FqFieldElem}}(undef, num_gens)

                for k in 1:num_gens
                    gens_vec[k] = action_morphism(stabilizers_gens[j][k], quotients[i-1], poly, inv_poly, x)
                end
                stabilizer_maps[j] = gens_vec
            end
        else
            @threads for j in eachindex(stabilizers_gens)
                num_gens = length(stabilizers_gens[j])
                gens_vec = Vector{AbstractAlgebra.Generic.ModuleHomomorphism{FqFieldElem}}(undef, num_gens)

                for k in 1:num_gens
                    gens_vec[k] = action_morphism(stabilizers_gens[j][k], (V, id), poly, inv_poly, x)
                end
                stabilizer_maps[j] = gens_vec
            end
        end
        orbits_and_stabilizers = collect(zip(orbits, stabilizer_maps))

    end

    if verbose == true
        println("Starting stage #3 -- finding reprsentatives in V")
    end

    size = Int(orbit_size(normal_forms(F, n), d))
    pi = quotients[1][2]
    V_fin = (filtration[end-1].object)
    poly_reps = final_lift_orbit_representatives(V_fin, poly, id, orbits_and_stabilizers, size)
    delete!(poly_reps, forget_grading(R(0)))

    return poly_reps
end
