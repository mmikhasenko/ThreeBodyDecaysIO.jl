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
