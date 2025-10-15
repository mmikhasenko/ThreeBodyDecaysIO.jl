# deserialization

@with_kw struct NamedArgFunc{M,D}
    f::M
    variable_names::D
end

function (obj::NamedArgFunc)(dict::AbstractDict)
    @unpack variable_names = obj
    variable_value = dict[first(variable_names)]
    obj.f(variable_value)
end

# ## Functions of s

function dict2instance(::Type{BreitWigner}, dict::AbstractDict)
    @unpack mass, width, ma, mb, l, d, x = dict
    bw = BreitWigner(mass, width, ma, mb, l, d)
    variables = [x]
    NamedArgFunc(bw, variables)
end

function dict2instance(::Type{MultichannelBreitWigner}, dict::AbstractDict)
    @unpack mass, channels, x = dict
    # convert dict channels into NamedTuple
    _channels = map(channels) do channel
        @unpack gsq, ma, mb, l, d = channel
        (; gsq, ma, mb, l, d)
    end
    bw = MultichannelBreitWigner(mass, _channels)
    variables = [x]
    return NamedArgFunc(bw, variables)
end

function dict2instance(::Type{Polynomial}, dict::AbstractDict)
    @unpack coefficients, x = dict
    P = Polynomial(coefficients, Symbol(x))
    arithmetic_P = WrapFlexFunction(P)
    variables = [x]
    NamedArgFunc(arithmetic_P, variables)
end

# ## Functions of p

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
    dict = LittleDict{String,Any}(
        "type" => type,
        "mass" => f.m,
        "width" => f.Γ,
        "ma" => ma,
        "mb" => mb,
        "l" => l,
        "d" => d,
        "x" => x,
    )
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
        LittleDict{String,Any}("gsq" => gsq, "ma" => ma, "mb" => mb, "l" => l, "d" => d)
    end
    appendix = Dict()
    dict = LittleDict{String,Any}(
        "type" => type,
        "mass" => m,
        "channels" => _channels,
        "x" => x,
    )
    return (dict, appendix)
end


function serializeToDict(x::BlattWeisskopf)
    l = orbital_momentum(x)
    radius = x.d
    type = "BlattWeisskopf"
    dict = LittleDict{String,Any}("type" => type, "l" => l, "radius" => radius)
    appendix = Dict()
    return (dict, appendix)
end
function serializeToDict(x::MomentumPower)
    l = orbital_momentum(x)
    type = "MomentumPower"
    dict = LittleDict{String,Any}("type" => type, "l" => l)
    appendix = Dict()
    return (dict, appendix)
end
function serializeToDict(H::NoFormFactor)
    appendix = Dict()
    ("", appendix)
end





# For lineshape is not just serializeToDict,
# but it's split into three parts:
# - scattering
# - FF_production
# - FF_decay
# and we need to serialize each part separately.
#
# For combined lineshape, we need to serialize the three parts separately.
# For other lineshapes, we use the trivial lineshape parser.

"""
    trivial_lineshape_parser(Xlineshape; k)

Trivial lineshape parser leaves all fields empty.
"""
function trivial_lineshape_parser(Xlineshape; k)
    scattering = ""
    FF_production = ""
    FF_decay = ""
    appendix = Dict()
    (; scattering, FF_production, FF_decay), appendix
end

struct CombinedFPF{T1,T2,T3} <: HadronicLineshapes.AbstractFlexFunc
    ff_production::T1
    ff_decay::T2
    propagator::T3
    ms::ThreeBodyDecays.MassTuple{Float64}
    k::Int
end
CombinedFPF(ff_production, ff_decay, propagator; ms, k) =
    CombinedFPF(ff_production, ff_decay, propagator, ms, k)
function (d::CombinedFPF)(σ)
    @unpack ms, k = d
    i, j = ij_from_k(k)
    p_production = HadronicLineshapes.breakup(ms[4], sqrt(σ), ms[k])
    p_decay = HadronicLineshapes.breakup(sqrt(σ), ms[i], ms[j])
    d.ff_production(p_production) * d.ff_decay(p_decay) * d.propagator(σ)
end


function lineshape_parser(Xlineshape; k)
    @warn "Using trivial lineshape parser for lineshape of type $(typeof(Xlineshape)). No custom parser is defined for type $(typeof(Xlineshape)). To ensure correct parsing, please implement a custom lineshape_parser method for this type. See the method for CombinedFPF as an example."
    trivial_lineshape_parser(Xlineshape; k) # fall back
end


"""
    lineshape_parser(Xlineshape::CombinedFPF; k)

A proper lineshape parser for lineshape is defined on a special type `CombinedFPF`, where
the propagator, FF_production, and FF_decay are separated.
The method generates a random name for the propagator, FF_production, and FF_decay with a common base name.
The index k is used to provide correct name to the dependence variable.
"""
function lineshape_parser(Xlineshape::CombinedFPF; k)
    appendix = Dict()

    @unpack ff_decay, ff_production, propagator = Xlineshape
    dependence_variable = "σ$k"
    #
    common_part = randstring(10)
    propagator_name = common_part * "_propagator"
    FF_decay_name = common_part * "(decay)"
    FF_production_name = common_part * "(production)"

    a = Dict(
        propagator_name =>
            NamedArgFunc(propagator, [dependence_variable]) |> serializeToDict |> first,
    )
    merge!(appendix, a)

    a = Dict(
        FF_decay_name => ff_decay |> serializeToDict |> first,
        FF_production_name => ff_production |> serializeToDict |> first,
    )
    merge!(appendix, a)

    (;
        scattering = propagator_name,
        FF_production = FF_production_name,
        FF_decay = FF_decay_name,
    ),
    appendix
end
