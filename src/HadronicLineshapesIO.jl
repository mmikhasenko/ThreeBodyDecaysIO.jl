
struct HS3InputWrapper{T}
    nt::NamedTuple{T}
end

# constuctor
function HS3InputWrapper(d::Dict)
    _d = copy(d)
    haskey(d, "type") && pop!(_d, "type")
    haskey(d, "name") && pop!(_d, "name")
    pairs = [k => v for (k, v) in _d]
    sort!(sort!(pairs; by=first); by=x -> length(x[1]), rev=true)
    names = Tuple(Symbol.(first.(pairs)))
    nt = NamedTuple{names}(getindex.(pairs, 2))
    HS3InputWrapper(nt)
end


function HadronicLineshapes.BreitWigner(wrapper::HS3InputWrapper{(:channels, :mass)})
    @unpack mass, channels = wrapper.nt
    _channels = map(channels) do channel
        @unpack gsq, ma, mb, l, d = channel
        (; gsq, ma, mb, l, d)
    end
    BreitWigner(mass, _channels)
end

function HadronicLineshapes.BreitWigner(wrapper::HS3InputWrapper{(:width, :mass, :ma, :mb, :d, :l)})
    @unpack mass, width, ma, mb, l, d = wrapper.nt
    BreitWigner(mass, width, ma, mb, l, d)
end

function dict2lineshape(fn)
    @unpack type = fn
    if type == "Flatte1405"
        @unpack mass, width = fn
        return BreitWigner(mass, width)
    elseif type == "BreitWignerWidthExp"
        @unpack mass, width = fn
        return BreitWigner(mass, width)
    elseif type == "BlattWeisskopf"
        @unpack radius = fn
        return (s, m1, m2, L) -> BlattWeisskopf{L}(radius * breakup(sqrt(s), m1, m2))
    end
    try
        constructor = Symbol(type)
        wrapper = HS3InputWrapper(fn)
        eval(quote
            $constructor($wrapper)
        end)
    catch
        error("Construction or the type, $type failed:\n")
    end
end

