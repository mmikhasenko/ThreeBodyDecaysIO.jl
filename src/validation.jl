
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

# reading validation

function angles_invariants(mass_angles_cascade, ms; k)
    _m0, _α, _cosβ = mass_angles_cascade[1]
    @assert _m0 == ms[4] "Error: The m0 value of `mass_angles_cascade` does not match m0 in `kinematics``."
    # 
    _mk, _γ, _cosθij = mass_angles_cascade[2]
    _σk = _mk^2
    # 
    _σj = σjofk(_cosθij, _σk, ms^2; k)
    _σi = sum(ms^2) - _σk - _σj
    σs = (_σi, _σj, _σk) |> reorder(k) |> MandestamTuple{Float64}
    (α=_α, cosβ=_cosβ, γ=_γ), σs
end





function validation_fields(model, point::NamedTuple{(:m, :cosθ)}; k, point_name, model_name)
    @unpack m, cosθ = point
    ms = masses(model)
    xσ = (m^2 - lims(k, ms)[1]) / (lims(k, ms)[2] - lims(k, ms)[1])
    σs = x2σs([xσ, (cosθ + 1) / 2], ms; k)
    validation_fields(model, σs; k, point_name, model_name)
end


function validation_fields(model, σs::MandestamTuple; k, point_name, model_name)
    i, j = ij_from_k(k)
    mk = sqrt(σs[k])
    # 
    ms = masses(model)
    cosθk = cosθij(σs, ms^2; k)
    value = unpolarized_intensity(model, σs)
    # 
    ij, ij_k = "$(i)$(j)", "$(i)$(j)_$(k)"
    point = LittleDict(
        "name" => point_name,
        "parameters" => [
            LittleDict("name" => "m_$(ij_k)", "value" => ms.m0),
            LittleDict("name" => "phi_$(ij_k)", "value" => 0.0),
            LittleDict("name" => "cos_theta_$(ij_k)", "value" => 0.0),
            LittleDict("name" => "m_$(ij)", "value" => mk),
            LittleDict("name" => "phi_$(ij)", "value" => 0.0),
            LittleDict("name" => "cos_theta_$(ij)", "value" => cosθk),
        ]
    )
    amplitude_model_checksums = LittleDict(
        "distribution" => model_name,
        "point" => point_name,
        "value" => value
    )
    (; point, amplitude_model_checksums)
end


function validation_section(model, validation_points; k, point_names, model_name)
    validation_dict = map(zip(point_names, validation_points)) do (point_name, point)
        validation_fields(model, point; k=3, point_name, model_name)
    end
    # 
    _dict = LittleDict()
    _dict[:misc] = Dict(
        :amplitude_model_checksums =>
            getfield.(validation_dict, :amplitude_model_checksums))
    # 
    _dict[:parameter_points] = []
    append!(_dict[:parameter_points],
        getfield.(validation_dict, :point))
    return _dict
end