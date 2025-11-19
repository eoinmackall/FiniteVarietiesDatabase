#####################################################################
#
#   Projective equivalence classes (distributed filtration method)
#
#####################################################################

#TODO: add explanation of this code, and the differences between
#this and mthreads

#####################################################################
#
#   Functions to run on worker processes
#
#####################################################################

@everywhere begin
    function action_for_stabilizer(A, v, p, poly, inv_poly, x)

        A = matrix(A)
        f = poly(preimage(p, v))
        y = A * x
        g = f(y...)
        g_vec = inv_poly(g)
        return p(g_vec)
    end

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

    function compute_stabilizer_worker(orbit, parent_gens, f_i, p_map, poly, inv_poly, x)
    #TODO: write this
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

    function final_lift_orbit_representatives(x, ImV_i, p, poly, q, inv_poly, f_i, orbits_and_stabilizers, size)
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
end

#####################################################################
#
#   Functions to run on master process
#
#####################################################################

function stabilizer_maker_init(G, orbits, action)
    stabilizer_gens = []
    for vec in orbits
        stab, _ = stabilizer(G, vec, action)
        push!(stabilizers, [matrix(g) for g in gens(stab)])
    end
    return stabilizers
end

function stabilizer_maker_recursive_distributed(orbits, orbits_and_stabilizers, f_i, action)

    parent_map_dict = Dict(pair[1] => pair[2] for pair in orbits_and_stabilizers)

    stabilizers_gens = pmap(orbit -> begin
        parent_vec = f_i(orbit)
        parent_gens = parent_map_dict[parent_vec]

        return compute_stabilizer_worker(orbit, parent_gens, f_i, p_map, poly, inv_poly, x)
    end, orbits; batch_size = max(1, div(length(orbits), nworkers()*4)))

    return stabilizer_gens
end
    
function projective_hypersurface_equivalence_classes_from_filtration_mprocs(F, n, d;
    waring_samples=48, basis_samples=100, verbose=false, interactive=false)


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
        println("Found chain with maximal relative dimension = ", rel_dim)
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

    stabilizers = stabilizer_maker_init(G, orbits, (v, A) -> action_for_stabilizer(A, v, f_2, poly, inv_poly, x))
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
    poly_reps = final_lift_orbit_representatives(x, V_fin, pi, poly, id, inv_poly, pi, orbits_and_stabilizers, size)
    delete!(poly_reps, forget_grading(R(0)))

    return poly_reps
end
