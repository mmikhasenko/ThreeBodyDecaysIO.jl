function dict2instance(::Type{BreitWigner}, dict)
    @unpack mass, width, ma, mb, l, d = dict
    BreitWigner(mass, width, ma, mb, l, d)
end

function dict2instance(::Type{BlattWeisskopf}, dict)
    @unpack radius, l = dict
    return BlattWeisskopf{l}(radius)
end

function dict2instance(::Type{MomentumPower}, dict)
    @unpack l = dict
    return MomentumPower{l}()
end

function dict2instance(::Type{MultichannelBreitWigner}, dict)
    @unpack mass, channels = dict
    # convert dict channels into NamedTuple
    _channels = map(channels) do channel
        @unpack gsq, ma, mb, l, d = channel
        (; gsq, ma, mb, l, d)
    end
    return MultichannelBreitWigner(mass, _channels)
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
function serializeToDict(x::BreitWigner)
    type = "BreitWigner"
    @unpack ma, mb, l, d = x
    dict = LittleDict{Symbol,Any}(pairs((; type, mass=x.m, width=x.Γ, ma, mb, l, d)))
    appendix = Dict()
    return (dict, appendix)
end