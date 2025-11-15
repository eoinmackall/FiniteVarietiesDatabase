#####################################################################
#
#   Projective equivalence classes (filtration method)
#
#####################################################################


#####################################################################
#
#   Projective equivalence classes maker
#
#####################################################################

function lift_orbit_representatives(x, ImV_i, p, poly, q, inv_poly, f_i, orbits, GL_vec)
    # W_i=V/V_i, f_i = V/V_{i+1}-> V/V_i, pi = V->W_i 
    # inc and inv_inc go between V and R_d = homogeneous_component of degree d
    # orbits = set in V/V_i
    
    new_orbits = Set()

    Im_i=ImV_i[2].((ImV_i)[1])

    orbits_chunks = Iterators.partition(orbits, cld(length(orbits), nthreads()))

    tasks = map(orbits_chunks) do chunk
        @spawn begin
            partial_orbit_reps = Set()
            for vec in chunk
                preim = preimage(p,vec)
                preim_poly = poly(preim)
                stabilizing_subgroup = Vector{FqMatrix}()
                for A in GL_vec
                    if is_stabilizing(A, vec, preim_poly, p, x, inv_poly)
                        push!(stabilizing_subgroup,A)
                    end
                end

                lift = preimage(f_i, vec) # preimage in V/V_{i+1}
                coset = Set(lift+v for v in Im_i)

                while !isempty(coset)
                    g=pop!(coset)
                    g_orbits = Set()
                    push!(partial_orbit_reps, g)
                    for A in stabilizing_subgroup
                        f=poly(preimage(q,g))
                        y=A*x
                        h=f(y...)
                        orbit_vec = q(inv_poly(h))
                        push!(g_orbits, orbit_vec)
                    end
                    setdiff!(coset,g_orbits)
                end
            end
            return partial_orbit_reps
        end
    end
    union!(new_orbits, (fetch.(tasks))...)
    return new_orbits
end

function final_lift_orbit_representatives(x, ImV_i, p, poly, q, inv_poly, f_i, orbits, GL_vec, size)
    # W_i=V/V_i, f_i = V/V_{i+1}-> V/V_i, pi = V->W_i 
    # inc and inv_inc go between V and R_d = homogeneous_component of degree d
    # orbits = set in V/V_i
    
    final_orbits = Set{FqMPolyRingElem}()
    sizehint!(final_orbits, size)

    Im_i=ImV_i[2].((ImV_i)[1])

    orbits_chunks = Iterators.partition(orbits, cld(length(orbits), nthreads()))

    tasks = map(orbits_chunks) do chunk
        @spawn begin
            partial_orbit_reps = Set()
            for vec in chunk
                preim = preimage(p,vec)
                preim_poly = poly(preim)
                stabilizing_subgroup = Vector{FqMatrix}()
                for A in GL_vec
                    if is_stabilizing(A, vec, preim_poly, p, x, inv_poly)
                        push!(stabilizing_subgroup,A)
                    end
                end

                lift = preimage(f_i, vec) # preimage in V/V_{i+1}
                coset = Set(lift+v for v in Im_i)

                while !isempty(coset)
                    g=pop!(coset)
                    g_orbits = Set()
                    push!(partial_orbit_reps, forget_grading(poly(g)))
                    for A in stabilizing_subgroup
                        f=poly(preimage(q,g))
                        y=A*x
                        h=f(y...)
                        orbit_vec = q(inv_poly(h))
                        push!(g_orbits, orbit_vec)
                    end
                    setdiff!(coset,g_orbits)
                end
            end
            return partial_orbit_reps
        end
    end
    union!(final_orbits, (fetch.(tasks))...)
    return final_orbits
end


function projective_hypersurface_equivalence_classes_from_filtration(F, a, n, d; waring_samples=48, basis_samples=100, verbose=false, interactive=false)
    

    head = _chain_constructor(F, a, n, d; waring_samples, basis_samples, verbose)
    if interactive == true
        chains = collect_chains(head)
        println("Choose a chain to use as filtration:")
        filtration = chains[parse(Int, readline())]
    else
        rel_dim, filtration = _chain_finder(head)
        # Considering the filtration ordered as V=V_1 > V_2 > V_3 > ... > V_{end-2} > V_{end-1} > 0

        j=0
        while length(filtration) == 2 || dim((filtration[end].object)[1]) != 0
            head = _chain_constructor(F, a, n, d; waring_samples, basis_samples, verbose)
            rel_dim, filtration = _chain_finder(head)
            if verbose == true
                println("Failed to find a good chain, trying again.")
            end
            if j == 10
                println("Failed to find invariant submodules. Is this module irreducible?")
                return
            end
            j+=1
        end
    end

    V, poly = head.object
    R=codomain(poly)
    x=gens(R)

    inv_poly = inv(poly)

    if verbose == true
        println("Found chain with maximal relative dimension = ", rel_dim)
        println("Beginning orbit collection")
    end
 
    quotients = []
    for i=1:length(filtration)-2
        push!(quotients,quo(V,(filtration[end-i].object)[1]))
        # Storing V/V_{end-i} at index i
        # So V/V_{end-1} to V/V_2
    end

    projection_maps = []
    for i=1:length(quotients)-1
        Source_quo, quo_source = quotients[i]
        Target_quo, quo_target = quotients[i+1]

        gen_images = [quo_target(preimage(quo_source, v)) for v in gens(Source_quo)]
        pi = ModuleHomomorphism(Source_quo, Target_quo, gen_images)
        push!(projection_maps, pi)
        # Storing f_i : V/V_{end-i} ->> V/V_{end-i-1} at index i
    end

    images = []
    for i=2:length(filtration)-2
        push!(images, image(compose((filtration[end-i].object)[2], quotients[i-1][2])))
        # Storing V_{end-i}/V_{end-i+1} at index i
    end


    if verbose == true
        println("Starting stage #1 -- finding orbits in V/V_1")
    end

    GL_vec = _GL(n + 1, F)
    orbits = Set()

    # Find orbits in V/V_2, add to orbits.
    W_2, f_2 = quotients[end]
    starting_vecs = Set(W_2)
    while !isempty(starting_vecs)
        vec = pop!(starting_vecs)
        push!(orbits,vec)
        orbit_of_vec = Set()
        
        lift = preimage(f_2, vec)
        f = poly(lift)
        
        GL_chunks = Iterators.partition(GL_vec, cld(length(GL_vec), nthreads()))

        tasks = map(GL_chunks) do chunk
            @spawn begin
                chunk_orbit = Set()
                for A in chunk
                    y = A*x
                    g = f(y...)
                    g_vec = inv_poly(g)
                    orbit_vec = f_2(g_vec)
                    push!(chunk_orbit, orbit_vec)
                end
                return chunk_orbit
            end
        end
        union!(orbit_of_vec, (fetch.(tasks))...)
        setdiff!(starting_vecs,orbit_of_vec)
    end

    if verbose == true
        println("Starting stage #2 -- lifting orbits along chain")
    end
    
    for i in (length(projection_maps)):-1:1
        ImV_i = images[i]
        f_i = projection_maps[i]
        p = quotients[i+1][2]
        q = quotients[i][2]
        orbits = lift_orbit_representatives(x, ImV_i, p, poly, q, inv_poly, f_i, orbits, GL_vec)
    end

    if verbose == true
        println("Starting stage #3 -- finding reprsentatives in V")
    end

    size = Int(orbit_size(normal_forms(F,n), d))
    pi = quotients[1][2]
    V_fin = (filtration[end-1].object)
    id = identity_map(V)
    poly_reps=final_lift_orbit_representatives(x, V_fin, pi, poly, id, inv_poly, pi, orbits, GL_vec, size)
    delete!(poly_reps, forget_grading(R(0)))

    return poly_reps
end
