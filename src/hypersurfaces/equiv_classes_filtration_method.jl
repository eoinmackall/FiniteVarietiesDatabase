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

function complete_poset(head::DChainNode, nodes::Set{DChainNode})

    for node1 in nodes
        for node2 in nodes
            if node1==node2
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
            good_chain = max(length, good_chain, chain)
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

    X, Y = GL_generators(F, a, n+1)

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


                max_samples = Int(min(2*BigInt(q^dim(L)), waring_samples))
                sample_range = collect(1:max_samples)
                sample_range_chunks = Iterators.partition(sample_range, cld(max_samples, Threads.nthreads()))

                tasks = map(sample_range_chunks) do chunk
                    Threads.@spawn begin
                        partial = Set()
                        for _ in chunk
                            g = phi(rand(L))
                            gens = Vector{AbstractAlgebra.Generic.FreeModuleElem{FqFieldElem}}()

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
    complete_poset(head, ALL)
    
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

function quotient_matrix(M, V, p, W)

    #M is a matrix, representing PGL element
    #V is vector space, W is quotient V/U, and p is projection V->W
    #iterate over basis for W:
    #1. pick preimage of basis element
    #2. do action of matrix in V
    #3. project down to W, get new vector
    #4. Assemble these vectors into a matrix
    
end

function projective_hypersurface_equivalence_classes_from_filtration(filt)

    #could possibly get n,d from filt?
    #
    # I think the following works:
    # take filtration 0<V_k<V_(k-1)<...<V_2<V_1<V
    # Compute V->V/V_k->V/V_(k-1)->...->V/V_2->V/V_1
    # Maybe put in a list
    #
    # Make a set S of 0 in V/V_1
    # 
    # new func(set_of_orbit_reps, Maybe V/V_(i+1)->V/V_i)
    #   for element in orbit reps set
    #       make PGL stabilizing subgroup
    #
    #   do orbit find on coset of V/V_i with stabilizing subgroup
    #   (returns a set. Splat the set into a bigger set).
    #
    #   return set_of_orbit_reps
    #
    #   for i=1:k+1
    #       S=func(S,V/V_(i+1)->V/V_i)
    #
    #   return S
    #   

end

#TODO: 1) Need a function that makes an action matrix on a quotient vector space
# 2) Need a function that runs the algorithm "recursively" on a filtration
# 3) Need a function filters through PGL and finds the stabilizer of an action matrix
