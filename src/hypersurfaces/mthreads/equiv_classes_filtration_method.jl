#####################################################################
#
#   Projective equivalence classes (filtration method)
#
#####################################################################

#TODO: make an introduction in here.

#TODO: Want to create a G-invariant filtration of subspaces in Sym^d(V), V=k{x_0,...,x_n}
# 1. Need a way to create G-invariant subspaces. Maybe this is done already?
# You can create these from weighted partitions. Should exclude trivial cases.
# Basically, given a partition lambda= l_1>l_2>l_3>... and weights w=w_1,w_2,w_3,...
# such that l_1*w_1+l_2*w_2+l_3*w_3+... = d, you can consider forms sum f_i^l_i where deg(f_i)=w_i
#
# 2. I think that I want to make a struct that keeps track of "chains" of subspaces.
# This way, we can just call a function to build this struct from the above. Then, we can look among
# the chains to find the best chain for performing the filtration algorithm.
# Might as well make this accept an initial object and a terminal object, and allow for chains created
# from there...

#####################################################################
#
#   Projective equivalence classes maker
#
#####################################################################

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

function orbit_stabilizer_pruning_optimized!(space::Set, G_gens, point, action_dict)

    orbit_to_process = [point]
    transversal = Dict(point.v => identity_matrix(base_ring(G_gens[1]), ncols(G_gens[1])))
    stab_gens = Set{FqMatrix}()

    gen_pairs = Vector{Tuple{eltype(G_gens), eltype(G_gens)}}(undef, length(G_gens))
    for (k, g) in enumerate(G_gens)
        gen_pairs[k] = (g, matrix(action_dict[g]))
    end

    M_parent = parent(point)
    R = base_ring(M_parent)
    temp_coords = zero_matrix(R, 1, ngens(M_parent)) 

    while !isempty(orbit_to_process)
        vec = pop!(orbit_to_process)
        vec_coords = vec.v
        trans_vec = transversal[vec_coords]

        for (g_original, M_linear) in gen_pairs
            mul!(temp_coords, vec_coords, M_linear)

            if !haskey(transversal, temp_coords)
                g_next_wrapped = M_parent(deepcopy(temp_coords))
                delete!(space, g_next_wrapped)
                push!(orbit_to_process, g_next_wrapped)
                transversal[g_next_wrapped.v] = trans_vec * g_original
            else
                trans_next = transversal[temp_coords]
                s = trans_vec * g_original * inv(trans_next)
                push!(stab_gens, s)
            end
        end
    end
    return collect(stab_gens)
end

function action_dictionary_parallel(G_gens, quotient, poly, inv_poly, x)

    results = Vector{Tuple{eltype(G_gens), AbstractAlgebra.Generic.ModuleHomomorphism{FqFieldElem}}}(undef, length(G_gens))

    Threads.@threads for i in eachindex(G_gens)
        A = G_gens[i]
        results[i] = (A, action_morphism(A, quotient, poly, inv_poly, x))
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

function lift_orbit_representatives2(x, ImV_i, poly, q, inv_poly, f_i, orbits_stabilizers_dict)
    # W_i=V/V_i, f_i = V/V_{i+1}-> V/V_i, pi = V->W_i 
    # inc and inv_inc go between V and R_d = homogeneous_component of degree d
    # orbits = set in V/V_i

    W=domain(f_i)
    new_orbits_stab_dict = Dict{elem_type(W), Vector{FqMatrix}}()

    Im_i = ImV_i[2].((ImV_i)[1])

    orbits = collect(keys(orbits_stabilizers_dict))

    tasks = map(orbits) do vec
        Threads.@spawn begin
            partial_orbit_stab_dict = Dict{elem_type(W), Vector{FqMatrix}}()
            lift = preimage(f_i, vec) # preimage in V/V_{i+1}
            coset = Set{elem_type(W)}(lift + v for v in Im_i)
            action_dict = action_dictionary(orbits_stabilizers_dict[vec],(W,q), poly, inv_poly, x)

            while !isempty(coset)
                vec_lift=pop!(coset)
                
                stab_vec_lift = orbit_stabilizer_pruning_optimized!(coset, orbits_stabilizers_dict[vec], vec_lift, action_dict)
                partial_orbit_stab_dict[vec_lift] = stab_vec_lift
            end
            return partial_orbit_stab_dict
        end
    end

    for t in tasks
        merge!(new_orbits_stab_dict, fetch(t))
    end

    return new_orbits_stab_dict
end

function final_lift_orbit_representatives2(x, ImV_i, poly, inv_poly, q, f_i, orbits_and_stabilizers, size)
    # W_i=V/V_i, f_i = V/V_{i+1}-> V/V_i, pi = V->W_i 
    # inc and inv_inc go between V and R_d = homogeneous_component of degree d
    # orbits = set in V/V_i

    W=domain(f_i)

    final_orbits = Set{FqMPolyRingElem}()
    sizehint!(final_orbits, size)

    Im_i = ImV_i[2].((ImV_i)[1])

    orbits = collect(keys(orbits_and_stabilizers))

    tasks = map(orbits) do vec
        Threads.@spawn begin
            partial_orbit = Set{FqMPolyRingElem}()
            lift = preimage(f_i, vec) # preimage in V/V_{i+1}
            coset = Set(lift + v for v in Im_i)
            action_dict = action_dictionary(orbits_and_stabilizers[vec], (W,q), poly, inv_poly, x)

            while !isempty(coset)
                vec_lift=pop!(coset)
                push!(partial_orbit, forget_grading(poly(vec_lift)))
                orbit_stabilizer_pruning_optimized!(coset, orbits_and_stabilizers[vec], vec_lift, action_dict)
            end
            return partial_orbit
        end
    end

    union!(final_orbits, (fetch.(tasks))...)
    return final_orbits
end


function projective_hypersurface_equivalence_classes_from_filtration(F, n, d; waring_samples=48, basis_samples=100, verbose=false, interactive=false)


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




function projective_hypersurface_equivalence_classes_from_filtration2(F, n, d; waring_samples=48, basis_samples=100, verbose=false, interactive=false)


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
    GL_gens = [matrix(g) for g in gens(G)]
    orbits = []

    # Find orbits in V/V_2, add to orbits.
    W_2, _ = quotients[end]
    starting_vecs = Set(W_2)
    action_dict = action_dictionary_parallel(GL_gens, quotients[end], poly, inv_poly, x)
    orbits_and_stabilizers = Dict{elem_type(W_2), Vector{FqMatrix}}()

    while !isempty(starting_vecs)
        vec = pop!(starting_vecs)
        
        stab_gens = orbit_stabilizer_pruning_optimized!(starting_vecs, GL_gens, vec, action_dict)
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
        orbits_and_stabilizers = lift_orbit_representatives2(x, ImV_i, poly, q, inv_poly, f_i, orbits_and_stabilizers)
    end

    if verbose == true
        println("Starting stage #3 -- finding reprsentatives in V")
    end

    size = Int(orbit_size(normal_forms(F, n), d))
    pi = quotients[1][2]
    V_fin = (filtration[end-1].object)
    poly_reps = final_lift_orbit_representatives2(x, V_fin, poly, inv_poly, id, pi, orbits_and_stabilizers, size)
    delete!(poly_reps, forget_grading(R(0)))

    return poly_reps
end
