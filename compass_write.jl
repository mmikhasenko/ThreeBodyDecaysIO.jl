using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.ThreeBodyDecays.PartialWaveFunctions
using Parameters
using DelimitedFiles
using DataFrames




# ## Lineshapes
using ThreeBodyDecaysIO.HadronicLineshapes

@with_kw struct BreitWignerRhoNoSqrt <: HadronicLineshapes.AbstractFlexFunc
    m::Float64
    Γ::Float64
    mπ::Float64
    d::Float64
end
function (bw::BreitWignerRhoNoSqrt)(σ)
    @unpack m, Γ, mπ, d = bw
    p, p0 = sqrt(σ / 4 - mπ^2), sqrt(m^2 / 4 - mπ^2)
    ff = BlattWeisskopf{1}(d)
    mΓ = m * Γ * p / p0 * ff(p)^2 / ff(p0)^2
    1 / (m^2 - σ - 1im * mΓ)
end

@with_kw struct KatchaevSigma <: HadronicLineshapes.AbstractFlexFunc
    sP::Vector{Float64}
    a::Vector{Float64}
    c::Vector{Float64}
end
function (bw::KatchaevSigma)(σ)
    mπ = 0.13956755
    mπ2 = mπ^2

    @unpack sP, a, c = bw
    M00 = 0.0
    for i = 1:length(sP)
        M00 += a[i] / (σ - sP[i])
    end
    # 
    mK = 0.493677
    mK0 = 0.497614
    scale = (σ / (4 * ((mK + mK0) / 2.0)^2)) - 1.0
    for (i, ci) in enumerate(c)
        M00 += scale^(i - 1) * ci
    end

    rho00 = sqrt(Kallen(σ, mπ2, mπ2)) / σ
    return 1.0 / (M00 - 1im * rho00)
end


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# pipi resonances


_fρ = let
    bw = BreitWignerRhoNoSqrt(; m=0.7685, Γ=0.1507, mπ=0.13956755, d=4.94)
    ff = BlattWeisskopf{1}(bw.d)(σ -> sqrt(σ / 4 - bw.mπ^2))
    Xlineshape = bw * ff * (1 / ff(bw.m^2))
    Xlineshape
end
# 
_ff2 = let
    bw = BreitWigner(; m=1.274, Γ=0.185, ma=0.13956755, mb=0.13956755, l=2, d=4.94)
    p(σ) = HadronicLineshapes.breakup(sqrt(σ), bw.ma, bw.mb)
    ff = BlattWeisskopf{bw.l}(bw.d)(p)
    Xlineshape = bw * ff * (bw.m * bw.Γ)
    Xlineshape
end
# 
_fρ3 = let
    bw = BreitWigner(; m=1.690, Γ=0.190, ma=0, mb=0, l=0, d=4.94)
    mπ = 0.13956755
    p(σ) = HadronicLineshapes.breakup(sqrt(σ), mπ, mπ)
    ff = BlattWeisskopf{3}(bw.d)(p)
    extra(σ) = sqrt(bw.m * bw.Γ * sqrt(σ))
    Xlineshape = bw * extra * ff
    Xlineshape
end
# 
_ff0_1500 = let
    bw = BreitWigner(; m=1.507, Γ=0.109, ma=0, mb=0, l=0, d=4.94)
    Xlineshape = bw * (bw.m * bw.Γ)
    Xlineshape
end
# 
_ff0_980 = MultichannelBreitWigner(; m=0.965,
    channels=[
        (gsq=1 / 0.965 * 0.165, ma=0.13956755, mb=0.13956755, l=0, d=1.0),
        (gsq=1 / 0.965 * 0.165 * 4.21, ma=0.493677, mb=0.493677, l=0, d=1.0)])
# 
_fσ = KatchaevSigma(sP=[-0.0074, 0.9828] .+ 1e-7, a=[0.1131, 0], c=[0.0337, -0.3185, -0.0942, -0.5927])

_fρ(1.1) ≈ -1.8472929896027317 + 0.6744244890043742im
_ff2(1.1) ≈ 0.30305342103185806 + 0.11181942166047641im
_fρ3(1.1) ≈ 0.1589433409235323 + 0.02906252876860443im
_ff0_1500(1.1) ≈ 0.13756331374474612 + 0.019296000940740514im
_ff0_980(1.1) ≈ -0.9212576583634419 + 2.1470398931994494im
_ff0_980(0.4 + 1im * nextfloat(0.0)) ≈ 0.7246075861113888 + 0.07865591658208379im
_ff0_980(0.13 + 1im * nextfloat(0.0)) ≈ 0.3947260628478383 + 0.01674732323566498im
_fσ(1.1) ≈ 0.10172436035046228 + 1.0273440332286132im
_fσ(0.135) ≈ (0.5840174783557608 + 0.26875840850408017im)

const isobarsV = [_fσ, _fρ, _ff2, _fρ3]
const isobarsS = [_fσ, _ff0_980, _ff0_1500]



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Gottfried-Jackson

function ϵWignerD(j0, ϵP, M, λ, (ϕ, cosθ, χ))
    M < 0 && return 0.0im
    W⁺ = wignerD(j0, M, λ, ϕ, cosθ, χ)
    W⁻ = wignerD(j0, -M, λ, ϕ, cosθ, χ)
    #
    n = (M == 0) ? 1 / 2 : 1 / sqrt(2)
    return (W⁺ - ϵP * ThreeBodyDecays.x"-1"^(j0 - M) * W⁻) * n
end

function gj_amplitude(process, σs, angles; ϵP, M)
    j = div(spins(process).two_h0, 2)
    sum(-j:j) do λ
        conj(ϵWignerD(j, ϵP, M, λ, angles)) *
        amplitude(process, σs, ThreeBodySpins(0, 0, 0; two_h0=x2(λ)))
    end
end


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


function build_compass_model(wave_description, m0)
    mπ = 0.13956755
    ms = ThreeBodyMasses(mπ, mπ, mπ; m0)
    # 
    @unpack J = wave_description
    tbs = ThreeBodySystem(ms, ThreeBodySpins(0, 0, 0; two_h0=J |> x2))

    @unpack S = wave_description
    l = wave_description.L
    @unpack name = wave_description

    bw_ff = (S ≥ 0) ? isobarsV[S+1] : isobarsS[1-S]
    j = S ≥ 0 ? S : 0
    # 
    q(σ) = HadronicLineshapes.breakup(ms.m0, sqrt(σ), ms.m1)
    ff_Rj = BlattWeisskopf{l}(1 / 0.2024)(q)
    Xlineshape = bw_ff * ff_Rj
    iϵ = 1im * nextfloat(0.0)

    dc1 = DecayChain(;
        k=1,
        two_j=x2(j),
        Xlineshape=σ -> Xlineshape(σ + iϵ),
        HRk=RecouplingLS((l, j) .|> x2),
        Hij=RecouplingLS((j, 0) .|> x2),
        tbs)

    # isobarsV[S+1]
    dc3 = DecayChain(dc1; k=3)
    c3 = sqrt(tbs.two_js.two_h0 + 1)
    c1 = c3 * ThreeBodyDecays.x"-1"^(j)

    wave = ThreeBodyDecay(name .=> zip([c1, c3], [dc1, dc3]))
    return wave
end



# comparison

phsp = (σ1=0.13569322768095665, cosθ1=0.5832472308560757, ϕ1=0.5079864049912346, cosθ23=-0.12538287914286417, ϕ23=-0.39836956124095346, s=2.3201214385414826)

σs_test = let
    @unpack σ1 = phsp
    σ2 = σ2of1(phsp.cosθ23, σ1, masses(wave2)^2)
    Invariants(masses(wave2); σ1, σ2)
end

angles_test = (ϕ=phsp.ϕ1, cosθ=phsp.cosθ1, χ=phsp.ϕ23)






wave_description = @NamedTuple{wn, name, J, P, M, ϵ, S, L}((2, "1-(1++)0+rhopiS", 1, "+", 0, "+", 1, 0))
value = 3.5036258938478007 - 0.6239732117186556im

wave2 = build_compass_model(wave_description, sqrt(phsp.s))


wave2.chains[1].Xlineshape(σs_test[1]) ≈ 1.288120896761017 + 0.03786584582224358im
wave2.chains[2].Xlineshape(σs_test[3]) ≈ -2.064664920993486 + 0.8309945337099337im

cal_test = gj_amplitude(wave2[1], σs_test, angles_test;
    wave_description.M, ϵP=(wave_description.ϵ == wave_description.P)) ≈
           1.8203662058242676 + 0.05351182972272878im

gj_amplitude(wave2, σs_test, angles_test;
    wave_description.M, ϵP=2 * (wave_description.ϵ == wave_description.P) - 1) ≈
value


unpolarized_intensity(wave2, σs_test)

check_points = map(eachslice(readdlm("amplitudes_point1.txt"); dims=1)) do x
    complex(x...)
end[2:end]

wavelist_df = DataFrame(readdlm("waves.txt"), [:wn, :name, :J, :P, :M, :ϵ, :S, :L])
wavelist_df = wavelist_df[2:end, :]


computed_values = map(eachrow(wavelist_df)) do wave_description
    wave = build_compass_model(wave_description, sqrt(phsp.s))
    gj_amplitude(wave, σs_test, angles_test;
        wave_description.M, ϵP=2 * (wave_description.ϵ == wave_description.P) - 1)
end

df_comp = DataFrame(; wavelist_df.name, computed_values, check_points,
    diff=computed_values - check_points, ratio=computed_values ./ check_points)

sort!(transform!(df_comp, :diff => ByRow(abs) => :absdiff), :absdiff)
