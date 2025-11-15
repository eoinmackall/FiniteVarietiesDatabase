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

############################################################
#
#  Struct for descending chains
#
############################################################

mutable struct DChainNode{T}

    object::T
    subobjects::Set{DChainNode}

    function DChainNode(node::T) where {T}

        new{T}(node, Set{DChainNode}())
    end
end

function add_subnode!(supernode::DChainNode, subnode::DChainNode)

    push!(supernode.subobjects, subnode)
end

function delete_subnode!(supernode::DChainNode, subnode::DChainNode)

    delete!(supernode.subobjects, subnode)
end

Base.:(==)(a::DChainNode{T}, b::DChainNode{S}) where {T,S} = (a.object[1] == b.object[1])

Base.hash(node::DChainNode{T}, h::UInt) where {T} = hash(node.object[1], h)

function is_in(head::DChainNode, node::DChainNode)

    return (head == node) || any(x -> is_in(x, node), head.subobjects)
end

function collect_chains(head::DChainNode)

    all_chains = Vector{Vector{DChainNode}}()

    function _chain_builder(node::DChainNode, chain::Vector{DChainNode})

        if isempty(node.subobjects)
            push!(all_chains, chain)
        else
            for subnode in node.subobjects
                chain_path = [chain..., subnode]
                _chain_builder(subnode, chain_path)
            end
        end
    end

    list = Vector{DChainNode}()
    push!(list, head)
    _chain_builder(head, list)
    return all_chains
end

############################################################
#
#  Construction of filtrations
#
############################################################

function _complete_poset(head::DChainNode, nodes::Set{DChainNode})

    for node1 in nodes
        for node2 in nodes
            if node1 == node2
                continue
            end
            has_intermediate_term = false
            intersection, _ = intersect(node1.object[1], node2.object[1])

            if intersection == node2.object[1]
                for node3 in nodes
                    if (node3 == node1 || node3 == node2)
                        continue
                    end
                    intersection1, _ = intersect(node1.object[1], node3.object[1])
                    intersection2, _ = intersect(node2.object[1], node3.object[1])
                    if (intersection1 == node3.object[1] && intersection2 == node2.object[1])
                        has_intermediate_term = true
                        break
                    end
                end
                if has_intermediate_term == true
                    continue
                else
                    add_subnode!(node1, node2)
                end
            end
        end
    end
    return head
end

function _chain_finder(head::DChainNode)

    chains = collect_chains(head)

    d = dim(head.object[1])
    good_chain = chains[1]
    for chain in chains
        temp_d = 0
        for i = 1:(length(chain)-1)
            temp_d = max(temp_d, dim((chain[i].object)[1]) - dim((chain[i+1].object)[1]))
        end
        if dim((chain[end].object)[1]) != 0
            temp_d = max(temp_d, dim((chain[end].object)[1]))
        end
        if temp_d < d
            d = temp_d
            good_chain = chain
        end
        if temp_d == d
            good_chain = length(chain)>length(good_chain) ? chain : good_chain
        end
    end
    return (d, good_chain)
end

function _polynomial_sampler(components, num_samples)

    len = length(components)
    num_component_samples = ceil(Int, ((1 + log(len)) * (num_samples))^(1 / len))
    rand_R_polys = []
    for V in components
        component_samples = []
        for _ in 1:num_component_samples
            push!(component_samples, V[2](rand(V[1])))
        end
        push!(rand_R_polys, component_samples)
    end
    return rand_R_polys

end

function _random_polynomial(components)

    return [V[2](rand(V[1])) for V in components]

end

"""
GL_generators(F::FqField, a::FqFieldElem, n::Int)

Return a pair of generators for the finite group GL(n,F),
depending on ``a``, a generator for the unit group of F.

Note: Calling F,a = finite_field(p,r,"a") provides a generator
``a`` so long as r>1.

If r=1, a primitive root can be gotten by calling primitive_root(p)
and converting to an element of F (e.g. F(primitive_root(p))).
"""
function GL_generators(F::FqField, a::FqFieldElem, n::Int)

    q = order(F)
    if q == 2
        X = identity_matrix(F, n)
        X[1, 2] = 1

        Y = zero_matrix(F, n, n)
        for i = 1:n-1
            Y[i+1, i] = 1
        end
        Y[1, n] = 1
    else
        X = identity_matrix(F, n)
        X[1, 1] = a

        Y = zero_matrix(F, n, n)
        for i = 1:n-1
            Y[i+1, i] = -1
        end
        Y[1, 1] = -1
        Y[1, n] = 1
    end
    return X, Y
end

function _is_GL_invariant(X, Y, x, W, sub_map, inc, inv_inc)

    W_basis = gens(W)
    for vec_f in W_basis
        f = inc(sub_map(vec_f))
        y = X * x
        g = f(y...)
        z = Y * x
        h = f(z...)

        bool1, _ = haspreimage(sub_map, inv_inc(g))
        bool2, _ = haspreimage(sub_map, inv_inc(h))
        if (!bool1 || !bool2)
            return false
        end
    end
    return true
end

function _chain_constructor(F, a, n, d; waring_samples=48, basis_samples=100, verbose=false)

    q = order(F)
    R, x = graded_polynomial_ring(F, ["x$i" for i = 0:n])
    Rd, inc = homogeneous_component(R, d)
    inv_inc = inv(inc)

    head = DChainNode((Rd, inc))

    X, Y = GL_generators(F, a, n + 1)

    ALL = Set{DChainNode}()
    seen_weights = []
    for P in partitions(d)
        len_p = length(P)
        divisors_list = []
        for i = 1:len_p
            push!(divisors_list, collect(divisors(P[i])))
        end
        divisors_iterator = Iterators.product(divisors_list...)
        for divs in divisors_iterator

            if length(P) == 1
                if divs == (d,)
                    continue
                end
            end

            if !(divs in seen_weights)
                len_divs = length(divs)
                S, _ = graded_polynomial_ring(F, ["w$i" for i = 1:len_divs], [divs...])
                L, phi = homogeneous_component(S, d)

                components = []
                for j in divs
                    R_j, i_j = homogeneous_component(R, j)
                    push!(components, (R_j, i_j))
                end

                max_basis_samples = min(2 * dim(Rd), basis_samples)
                rand_R_polys = _polynomial_sampler(components, max_basis_samples)

                Waring_iterator = Iterators.product(rand_R_polys...)


                max_samples = Int(min(2 * BigInt(q^dim(L)), waring_samples))
                sample_range = collect(1:max_samples)
                sample_range_chunks = Iterators.partition(sample_range, cld(max_samples, nthreads()))

                tasks = map(sample_range_chunks) do chunk
                    Threads.@spawn begin
                        partial = Set()
                        gens = Vector{AbstractAlgebra.Generic.FreeModuleElem{FqFieldElem}}()
                        for _ in chunk
                            g = phi(rand(L))

                            for f in Waring_iterator
                                push!(gens, inv_inc(evaluate(g, [f...])))
                            end

                            W, sub_map = sub(Rd, gens)

                            while !(_is_GL_invariant(X, Y, x, W, sub_map, inc, inv_inc))
                                f = _random_polynomial(components)
                                push!(gens, inv_inc(evaluate(g, f)))
                                W, sub_map = sub(Rd, gens)
                            end
                            W_node = DChainNode((W, sub_map))
                            push!(partial, W_node)
                            empty!(gens)
                        end
                        return partial
                    end
                end
                partial_spaces = fetch.(tasks)
                union!(ALL, partial_spaces...)
            end
            push!(seen_weights, divs)
        end
    end

    push!(ALL, head)
    _complete_poset(head, ALL)

    if verbose == true
        chains = collect_chains(head)
        println("Found the following chains:")
        for chain in chains
            for i = 1:length(chain)
                print(dim(chain[i].object[1]), " ")
            end
            println()
        end
    end
    return head
end

#####################################################################
#
#   Projective equivalence classes maker
#
#####################################################################

function is_stabilizing(A, vec, preim_poly, p, x, inv_poly)
    # Function to check if a matrix A stabilizes a vector vec

    y = A*x
    h=preim_poly(y...)
    h_vec=inv_poly(h)
    if vec == p(h_vec)
        return true
    else
        return false
    end
end

function lift_orbit_representatives(x, ImV_i, p, poly, q, inv_poly, f_i, orbits, GL_vec)
    # W_i=V/V_i, f_i = V/V_{i+1}-> V/V_i, pi = V->W_i 
    # inc and inv_inc go between V and R_d = homogeneous_component of degree d
    # orbits = set in V/V_i
    
    new_orbits = Set()

    Im_i=ImV_i[2].((ImV_i)[1])

    orbits_chunks = Iterators.partition(orbits, cld(length(orbits), nthreads()))

    tasks = map(orbits_chunks) do chunk
        Threads.@spawn begin
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
        Threads.@spawn begin
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
            Threads.@spawn begin
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
