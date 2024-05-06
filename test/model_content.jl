using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.JSON
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.DataFrames
using Test


# extra likeshapes for testing
@with_kw struct BreitWignerWidthExpLikeBugg <: HadronicLineshapes.AbstractFlexFunc
    m::Float64
    Γ::Float64
    γ::Float64
end
function (BW::BreitWignerWidthExpLikeBugg)(σ)
    mK = 0.493677
    mπ = 0.13957018
    σA = mK^2 - mπ^2 / 2
    @unpack m, Γ, γ = BW
    Γt = (σ - σA) / (m^2 - σA) * Γ * exp(-γ * σ)
    1 / (m^2 - σ - 1im * m * Γt)
end
function ThreeBodyDecaysIO.dict2instance(::Type{BreitWignerWidthExpLikeBugg}, dict)
    @unpack mass, width, slope = dict
    return BreitWignerWidthExpLikeBugg(mass, width, slope)
end


list_of_models = readdir(joinpath(@__DIR__, "..", "models"))

map(list_of_models) do file_name
    @info "⭐ Reading model from $file_name ⭐"
    # @testset "Test of $file_name" begin
    let
        @show "Test of $file_name"
        # get the JSON content
        input = open(joinpath(@__DIR__, "..", "models", file_name)) do io
            JSON.parse(io)
        end
        model_descrition = first(input["distributions"])
        @test model_descrition["type"] == "HadronicUnpolarizedIntensity"

        @test haskey(model_descrition, "decay_description")
        @unpack decay_description = model_descrition

        @test haskey(decay_description, "reference_topology")
        @unpack reference_topology = decay_description

        @test flatten_topology(reference_topology) |> sort == [1, 2, 3]

        @info "🐰 Parcing kinamatics 🐰"
        @unpack kinematics = decay_description
        tbs = dict2instance(ThreeBodySystem, kinematics)


        @test haskey(input, "functions")
        @unpack functions = input

        @info "🔥 Building lineshapes functions 🔥"
        workspace = Dict{String,Any}()
        for fn in functions
            @unpack name, type = fn
            instance_type = eval(Symbol(type))
            workspace[name] = dict2instance(instance_type, fn)
        end
        # @test dict2instance(ThreeBodyDecay, decay_description; workspace) isa ThreeBodyDecay

        @info "🦊 Establishing decay chains 🦊"
        df = dict2instance.(DecayChain, decay_description["chains"]; tbs, workspace) |> DataFrame
        model = ThreeBodyDecay(Vector{Pair{String,Tuple{Complex,AbstractDecayChain}}}(df.name .=> zip(df.coupling, df.chain)))

        @test model isa ThreeBodyDecay
    end
end



