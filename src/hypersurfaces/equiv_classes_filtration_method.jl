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

==(a::DChainNode, b::DChainNode) = (a.object == b.object)

hash(node:: DChainNode, h::UInt) = hash(node.object, h)

function is_in(head::DChainNode, node::DChainNode)
    
    return (head==node) || any(x->is_in(x,node), head.subobjects)
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

        Intersection = intersect(new_node.object, node.object)

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

    added_later = false
    for node in head.subobjects
        Intersection = intersect(new_node.object, node.object)
        if Intersection == new_node.object
            add_subspace_descending!(node, new_node)
            added_later = true
        end
    end
    
    if !added_later
        add_subnode!(head, new_node)
    end
    return
end

function _print_chain_dimensions(head::DChainNode)
    
    println(dim(head.object))
    for node in head.subobjects
        _print_chain_dimensions(node)
    end
    return
end

function _chain_finder(head::DChainNode)

    

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

function _is_GL_invariant(X, Y, x, V, W, sub_map, inc, inv_inc)

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

function _chain_constructor(F, a, n, d; waring_samples=48, basis_samples=100, verbose = false)

    q = order(F)
    R, x = graded_polynomial_ring(F, ["x$i" for i = 0:n])
    Rd, inc = homogeneous_component(R, d)
    inv_inc = inv(inc)

    head = DChainNode(Rd)

    X, Y = GL_generators(F, a, n)

    ALL = Set()
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

                            while !(_is_GL_invariant(X, Y, x, Rd, W, sub_map, inc, inv_inc))
                                f = _random_polynomial(components)
                                push!(gens, inv_inc(evaluate(g, f)))
                                W, sub_map = sub(Rd, gens)
                            end
                            W_node = DChainNode(W)
                            push!(partial,W_node)
                        end
                        return partial
                    end
                end
                partial_spaces = fetch.(tasks)
                union!(ALL, partial_spaces...)
                for node in ALL
                    add_subspace_descending!(head,node)
                end
            end
            push!(seen_weights, divs)
        end
    end
    if verbose == true
        _print_chain_dimensions(head)
    end
    return head
end
