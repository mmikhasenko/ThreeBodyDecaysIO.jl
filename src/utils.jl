function update2values(x, ref)
    if x isa String && haskey(ref, x)
        return ref[x]
    end
    if x isa Dict
        _d = Dict{String,Any}()
        for (k, v) in x
            _d[k] = update2values(v, ref)
        end
        return _d
    end
    if x isa AbstractArray
        return update2values.(x, Ref(ref))
    end
    return x
end

string2complex(s) = eval(Meta.parse(replace(s, "i" => "im")))

flatten_topology(topology) =
    topology isa Array ? vcat(flatten_topology.(topology)...) : topology

function topology2k(topology::AbstractArray)
    indices = fully_flatten(topology)
    length(indices) != 3 && error("Topology with more that three particles is not implemented")
    !(topology ∈ ([[1, 2], 3], [[2, 3], 1], [[3, 1], 2])) &&
        error("Only regular topologies [[1, 2], 3], [[2, 3], 1], [[3, 1], 2] are implememnted")
    # 
    return topology[2]
end

function reorder(last_index)
    k = last_index
    i, j = ij_from_k(k)
    t -> invpermute!(collect(t), [i, j, k]) |> Tuple
end

function array2dict(a::AbstractArray, key_of_key)
    map(a) do p
        _p = copy(p)
        key = p[key_of_key]
        pop!(_p, key_of_key)
        key => _p
    end |> Dict
end