#####################################################################
#
#   Projective equivalence classes (distributed filtration method)
#
#####################################################################

# Warning: this should not be used. At the moment, it will fail
# (probably by a seg-fault).
#
# This code was a first attempt at changing some of the functions
# for finding hypersurface representatives to work distributed-ively.
# However, I'm not sure if there is functionality for passing Oscar
# objects to remote processes. This could potentially be bypassed
# by converting to elementary types, but I don't know if I'll come
# back to this for a bit.

#####################################################################
#
#   Functions to run on worker processes
#
#####################################################################

function worker_action_for_stabilizer(A, v, p, poly, inv_poly, x)
    A = matrix(A)
    f = poly(preimage(p, v))
    y = A * x
    g = f(y...)
    g_vec = inv_poly(g)
    return p(g_vec)
end

function compute_stabilizer_worker(orbit, parent_gens, p_map, poly, inv_poly, x)
    
    R = base_ring(parent_gens[1])
    n = nrows(parent_gens[1])
    G = matrix_group(R,n,parent_gens)

    stab, _ = stabilizer(G, orbit, (v,A) -> worker_action_for_stabilizer(A, v, p_map, poly, inv_poly, x))

    return [matrix(g) for g in gens(stab)]
end

function lift_orbit_worker(pair, Im_i, f_i)
    
    vec, induced_gens = pair

    W=domain(f_i)
    
    partial_reps = []
    
    lift = preimage(f_i, vec)
    coset = Set(lift + v for v in Im_i)

    while !isempty(coset)
        g = pop!(coset)
        push!(partial_reps, g)

        g_orbits = Set()
        orbits_to_process = Set([g])
        
        while !isempty(orbits_to_process)
            g_vec = pop!(orbits_to_process)
            for A_mat in induced_gens
                pre_orb = A_mat*g_vec.v
                orbit_vec = W(pre_orb)
                
                if !(orbit_vec in g_orbits)
                    push!(g_orbits, orbit_vec)
                    push!(orbits_to_process, orbit_vec)
                end
            end
        end
        setdiff!(coset, g_orbits)
    end
    return partial_reps
end

function final_lift_worker(pair, Im_i_vals, pi_lift, id, poly, inv_poly, x)
    vec = pair[1]
    stabilizing_subgroup_gens = pair[2]
    
    action_func = (v, A) -> worker_action_for_stabilizer(A, v, id, poly, inv_poly, x)
    
    local_final_orbits = Set{FqMPolyRingElem}()
    
    # Lift from quotient to V
    lift = preimage(pi_lift, vec)
    coset = Set(lift + v for v in Im_i_vals)

    while !isempty(coset)
        g = pop!(coset)
        
        # Store the polynomial representation
        push!(local_final_orbits, forget_grading(poly(g)))

        g_orbits = Set()
        orbits_to_process = Set([g])
        while !isempty(orbits_to_process)
            g_vec = pop!(orbits_to_process)
            for A_mat in stabilizing_subgroup_gens
                
                # Apply action in V
                orbit_vec = action_func(g_vec, A_mat)
                
                if !(orbit_vec in g_orbits)
                    push!(g_orbits, orbit_vec)
                    push!(orbits_to_process, orbit_vec)
                end
            end
        end
        setdiff!(coset, g_orbits)
    end
    return local_final_orbits
end

#####################################################################
#
#   Functions to run on master process
#
#####################################################################

function action_matrix_for_worker(A, quotient, poly, inv_poly, x)

    W = quotient[1]
    pi = quotient[2]
    F = base_ring(W)
    k = dim(W)

    action_matrix = Vector{Vector{FqFieldElem}}()
    sizehint!(action_matrix, k)
    for v in gens(W)
        v_vec = preimage(pi, v)
        f = poly(v_vec)
        y = A * x
        g = f(y...)
        w_vec = inv_poly(g)
        w = pi(w_vec)
        push!(action_matrix, vec(collect(w.v)))
    end
    return matrix(F, k, k, reduce(vcat,action_matrix))
end


function stabilizer_maker_init(G, orbits, action)
    stabilizers_gens = []
    for vec in orbits
        stab, _ = stabilizer(G, vec, action)
        push!(stabilizers_gens, [matrix(g) for g in gens(stab)])
    end
    return stabilizers_gens
end

function stabilizer_maker_recursive_distributed(orbits, orbits_and_stabilizers, f_i, p_map, poly, inv_poly, x)

    parent_map_dict = Dict(pair[1] => pair[2] for pair in orbits_and_stabilizers)

    stabilizers_gens = pmap(orbit -> begin
        parent_vec = f_i(orbit)
        parent_gens = parent_map_dict[parent_vec]
        return compute_stabilizer_worker(orbit, parent_gens, p_map, poly, inv_poly, x)
    end, orbits; batch_size = max(1, div(length(orbits), nworkers()*4)))

    return stabilizers_gens
end

function lift_orbit_reps_distributed(ImV_i, f_i, orbits_and_stabilizers, poly, p, inv_poly, x)

    W_i = codomain(p)

    Im_i = ImV_i[2].((ImV_i)[1])

    action_data = Vector{Tuple{elem_type(W_i), Vector{FqMatrix}}}(undef, length(orbits_and_stabilizers))
    
    for (i,pair) in enumerate(orbits_and_stabilizers)
        vec, stabilizer = pair
        morphisms = [action_matrix_for_worker(M, (W_i,p), poly, inv_poly, x) for M in stabilizer]
        action_data[i]=(vec,morphisms)
    end

    results = pmap(pair -> lift_orbit_worker(pair, Im_i, f_i),
                   action_data;
                   batch_size = max(1,div(length(action_data), nworkers()*4)))

    new_orbits = []
    for result in results
        append!(new_orbits, result)
    end
    return new_orbits
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
    GL_gens = [action_matrix_for_worker(matrix(g), quotients[end], poly, inv_poly, x) for g in gens(G)]
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

            vec = pop!(orbit_to_process)
            for A in GL_gens
                h_vecv = A*vec.v
                h_vec = W_2(h_vecv)
                if !(h_vec in orbit_of_vec)
                    push!(orbit_of_vec, h_vec)
                    push!(orbit_to_process, h_vec)
                end

            end
        end
        setdiff!(starting_vecs, orbit_of_vec)
    end

    stabilizers = stabilizer_maker_init(G, orbits, (v, A) -> worker_action_for_stabilizer(A, v, f_2, poly, inv_poly, x))
    orbits_and_stabilizers = collect(zip(orbits, stabilizers))

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
        orbits = lift_orbit_reps_distributed(ImV_i, f_i, orbits_and_stabilizers, poly, p, inv_poly, x)

        if verbose == true
            println("Computing stabilizers")
        end

        stabilizers_gens = stabilizer_maker_recursive_distributed(orbits, orbits_and_stabilizers, f_i, q, poly, inv_poly, x)
        orbits_and_stabilizers = collect(zip(orbits, stabilizers_gens))
    end

    if verbose == true
        println("Starting stage #3 -- finding reprsentatives in V")
    end

    size = Int(orbit_size(normal_forms(F, n), d))
    ImV_fin = (filtration[end-1].object)
    pi_fin = quotients[1][2] 
    
    Im_i_vals = ImV_fin[2].((ImV_fin)[1])

    final_results = pmap(pair -> final_lift_worker(pair, Im_i_vals, pi_fin, id, poly, inv_poly, x),
                         orbits_and_stabilizers;
                         batch_size=max(1, div(length(orbits_and_stabilizers), nworkers() * 4)))

    final_orbits = Set{FqMPolyRingElem}()
    sizehint!(final_orbits, size)
    union!(final_orbits, final_results...)
    delete!(final_orbits, forget_grading(R(0)))

    return final_orbits
end
