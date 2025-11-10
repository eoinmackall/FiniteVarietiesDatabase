#####################################################################
#
#   Projective equivalence classes (filtration method)
#
#####################################################################

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

Base.:(==)(a::DChainNode{T}, b::DChainNode{S}) where {T,S} = (a.object == b.object)

Base.hash(node::DChainNode{T}, h::UInt) where {T} = hash(node.object, h)

function is_in(head::DChainNode, node::DChainNode)

    return (head == node) || any(x -> is_in(x, node), head.subobjects)
end

function _prune_and_collect_subspaces(
    head::Union{DChainNode{AbstractAlgebra.Generic.FreeModule{FqFieldElem}},DChainNode{AbstractAlgebra.Generic.Submodule{FqFieldElem}}},
    node::DChainNode{AbstractAlgebra.Generic.Submodule{FqFieldElem}})

    node_contains = Set{DChainNode{AbstractAlgebra.Generic.Submodule{FqFieldElem}}}()
    _prune_and_collect_subspaces!(node_contains, head, node)
    return node_contains
end

function _prune_and_collect_subspaces!(result::Set{DChainNode{AbstractAlgebra.Generic.Submodule{FqFieldElem}}},
    head::Union{DChainNode{AbstractAlgebra.Generic.FreeModule{FqFieldElem}},DChainNode{AbstractAlgebra.Generic.Submodule{FqFieldElem}}},
    new_node::DChainNode{AbstractAlgebra.Generic.Submodule{FqFieldElem}})

    temp = Set()
    for node in head.subobjects

        Intersection, _ = intersect(new_node.object, node.object)

        if Intersection == node.object
            push!(result, node)
            push!(temp, node)
        else
            _prune_and_collect_subspaces!(result, node, new_node)
        end
    end
    for subnode in temp
        delete_subnode!(head, subnode)
    end
end

function add_subspace_descending!(
    head::Union{DChainNode{AbstractAlgebra.Generic.FreeModule{FqFieldElem}},DChainNode{AbstractAlgebra.Generic.Submodule{FqFieldElem}}},
    new_node::DChainNode{AbstractAlgebra.Generic.Submodule{FqFieldElem}})

    if is_in(head, new_node)
        return
    end

    new_node_contains = _prune_and_collect_subspaces(head, new_node)
    for node in new_node_contains
        add_subnode!(new_node, node)
    end

    for node in head.subobjects
        Intersection, _ = intersect(new_node.object, node.object)
        if Intersection == new_node.object
            add_subspace_descending!(node, new_node)
        end
    end

    add_subnode!(head, new_node)
    return
end

function complete_poset(head::DChainNode, nodes::Set{DChainNode})

    for node1 in nodes
        for node2 in nodes
            if node1==node2
                continue
            end
            has_intermediate_term = false
            intersection, _ = intersect(node1.object, node2.object)

            if intersection == node2.object
                for node3 in nodes
                    if (node3 == node1 || node3 == node2)
                        continue
                    end
                    intersection1, _ = intersect(node1.object, node3.object)
                    intersection2, _ = intersect(node2.object, node3.object)
                    if (intersection1 == node3.object && intersection2 == node2.object)
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

    d = dim(head.object)
    good_chain = chains[1]
    for chain in chains
        temp_d = 0
        for i = 1:(length(chain)-1)
            temp_d += dim(chain[i].object) - dim(chain[i+1])
        end
        if dim(chain[end].object) != 0
            temp_d += dim(chain[end].object)
        end
        if temp_d < d
            d = temp_d
            good_chain = chain
        end
    end
    return (d, good_chain)
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

function _polynomial_sampler(components, num_samples)

    len = length(components)
    num_component_samples = ceil(Int, (1 + log(len)) * (num_samples)^(1 / len))
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

function GL_generators(F, a, n)

    q = order(F)
    if q == 2
        X = identity_matrix(F, n + 1)
        X[1, 2] = 1

        Y = zero_matrix(F, n + 1, n + 1)
        for i = 1:n
            Y[i+1, i] = 1
        end
        Y[1, n+1] = 1
    else
        X = identity_matrix(F, n + 1)
        X[1, 1] = a

        Y = zero_matrix(F, n + 1, n + 1)
        for i = 1:n
            Y[i+1, i] = -1
        end
        Y[1, 1] = -1
        Y[1, n+1] = 1
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

function _chain_constructor(F, a, n, d; waring_samples=48, basis_samples=100, verbose=false, complete=false)

    q = order(F)
    R, x = graded_polynomial_ring(F, ["x$i" for i = 0:n])
    Rd, inc = homogeneous_component(R, d)
    inv_inc = inv(inc)

    head = DChainNode(Rd)

    X, Y = GL_generators(F, a, n)

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


                max_samples = Int(min(BigInt(q^dim(L)), waring_samples))
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
                            W_node = DChainNode(W)
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
    if complete == false
        for node in ALL
            add_subspace_descending!(head, node)
        end
    else
        push!(ALL, head)
        complete_poset(head, ALL)
    end

    if verbose == true
        chains = collect_chains(head)
        println("Found the following chains:")
        for chain in chains
            for i = 1:length(chain)
                print(dim(chain[i].object), " ")
            end
            println()
        end
    end
    return head
end
