# deserialization

@with_kw struct NamedArgFunc{T,P,D}
    f::T
    variable_names::D
    parameters::P = String[] # currently empty
end


result = replace(["a", "a"], "a" => 1)

function (obj::NamedArgFunc)(dict::AbstractDict)
    @unpack variable_names = obj
    variable_value = dict[first(variable_names)]
    obj.f(variable_value)
end

function dict2instance(::Type{BreitWigner}, dict)
    @unpack mass, width, ma, mb, l, d, x = dict
    bw = BreitWigner(mass, width, ma, mb, l, d)
    parameters = String[]
    variables = [x]
    NamedArgFunc(bw, variables, parameters)
end

function dict2instance(::Type{MultichannelBreitWigner}, dict)
    @unpack mass, channels, x = dict
    # convert dict channels into NamedTuple
    _channels = map(channels) do channel
        @unpack gsq, ma, mb, l, d = channel
        (; gsq, ma, mb, l, d)
    end
    bw = MultichannelBreitWigner(mass, _channels)
    parameters = String[]
    variables = [x]
    return NamedArgFunc(bw, variables, parameters)
end

function dict2instance(::Type{BlattWeisskopf}, dict)
    @unpack radius, l = dict
    return BlattWeisskopf{l}(radius)
end

function dict2instance(::Type{MomentumPower}, dict)
    @unpack l = dict
    return MomentumPower{l}()
end



# serialization

function serializeToDict(obj::NamedArgFunc{BreitWigner})
    type = "BreitWigner"
    @unpack f, variable_names = obj
    @unpack ma, mb, l, d = f
    x = first(variable_names)
    dict = LittleDict{Symbol,Any}(pairs((; type, mass=f.m, width=f.Γ, ma, mb, l, d, x)))
    appendix = Dict()
    return (dict, appendix)
end
function serializeToDict(obj::NamedArgFunc{<:MultichannelBreitWigner})
    type = "MultichannelBreitWigner"
    @unpack f, variable_names = obj
    @unpack m, channels = f
    x = first(variable_names)
    _channels = map(channels) do channel
        @unpack gsq, ma, mb, l, d = channel
        LittleDict{Symbol,Any}(pairs((; gsq, ma, mb, l, d)))
    end
    appendix = Dict()
    dict = LittleDict{Symbol,Any}(pairs((; type, mass=m, channels=_channels, x)))
    return (dict, appendix)
end


function serializeToDict(x::BlattWeisskopf)
    l = orbital_momentum(x)
    radius = x.d
    type = "BlattWeisskopf"
    dict = LittleDict{Symbol,Any}(pairs((; type, l, radius)))
    appendix = Dict()
    return (dict, appendix)
end
function serializeToDict(x::MomentumPower)
    l = orbital_momentum(x)
    type = "MomentumPower"
    dict = LittleDict{Symbol,Any}(pairs((; type, l)))
    appendix = Dict()
    return (dict, appendix)
end
