
function cascade_nodes(σs, ms; k)
    i, j = ij_from_k(k)
    σk = σs[k]
    _cosθij = cosθij(σs, ms^2; k)
    Dict(
        "kinematic_point" => [
            Dict(
                "node" => [[i, j], k],
                "phi" => 0.0,
                "theta" => 0.0,
                "mass" => ms.m0),
            Dict(
                "node" => [i, j],
                "phi" => 0.0,
                "theta" => _cosθij,
                "mass" => sqrt(σk)),
        ])
end

function fully_flatten(x)
    _x = copy(x)
    i = 0
    while !(_x isa AbstractVector{<:Number})
        i += 1
        _x = collect(Iterators.flatten(_x))
        i > 10 && break
    end
    return _x
end

function validation_section(model, dpp, topology)
    k = topology2k(topology)
    ms = masses(model)
    @unpack σs, two_λs = dpp
    # 
    amplitudes = map(model.chains) do chain
        value = amplitude(chain, σs, two_λs)
        replace(string(value), "im" => "i")
    end
    # 
    _validation = Dict{Symbol,Any}(
        :kinematic_point => Dict(
            :masses_angles => cascade_nodes(σs, ms; k),
            :spin_projections => values(two_λs) .|> d2,
        ),
        :chains_values => amplitudes,
    )
    return _validation
end


# reading validation

function angles_invariants(mass_angles_cascade, ms; k)
    _m0, _α, _cosβ = mass_angles_cascade[1]
    @assert _m0 == ms[4] "Error: The m0 value of `mass_angles_cascade` does not match m0 in `kinematics``."
    # 
    _mk, _γ, _cosθij = mass_angles_cascade[2]
    _σk = _mk^2
    # 
    _σi = σiofk(_cosθij, _σk, ms^2; k)
    _σj = sum(ms^2) - _σk - _σi
    σs = (_σi, _σj, _σk) |> reorder(k) |> MandestamTuple{Float64}
    (α=_α, cosβ=_cosβ, γ=_γ), σs
end
