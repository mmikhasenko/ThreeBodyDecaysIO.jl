
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


function HadronicLineshapes.BlattWeisskopf(wrapper::HS3InputWrapper{(:radius, :l)})
    @unpack radius, l = wrapper.nt
    return BlattWeisskopf{l}(radius)
end

function dict2lineshape(fn)
    @unpack type = fn
    if type == "BreitWignerWidthExp"
        @unpack mass, width = fn
        return BreitWigner(mass, width)
    end
    # 
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


@with_kw struct BreitWignerWidthExp <: HadronicLineshapes.AbstractFlexFunc
    m::Float64
    Γ::Float64
    γ::Float64
end


function BreitWignerWidthExp(wrapper::HS3InputWrapper{(:slope, :width, :mass)})
    @unpack mass, width, slope = wrapper
    return BreitWignerWidthExp(mass, width, slope)
end


function (BW::BreitWignerWidthExp)(σ)
    mK = 0.493677
    mπ = 0.13957018
    σA = mK^2 - mπ^2 / 2
    m, Γ, γ = BW.pars
    Γt = (σ - σA) / (m^2 - σA) * Γ * exp(-γ * σ)
    1 / (m^2 - σ - 1im * m * Γt)
end