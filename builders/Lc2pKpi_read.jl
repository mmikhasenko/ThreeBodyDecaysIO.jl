# Script builds the model from JSON file

using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using HadronicLineshapes
using JSON
using Parameters
using DataFrames
using Plots
using Test



# functions
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


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


# code
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# get the JSON content
input = open(joinpath(@__DIR__, "..", "models", "Lc2ppiK.json")) do io
    JSON.parse(io)
end

# build functions from JSON array,
@unpack functions = input
# built functions will be stored in workspace
workspace = Dict{String,Any}()
# add into dictionary with the name as key
for fn in functions
    @unpack name, type = fn
    instance_type = eval(Symbol(type))
    workspace[name] = dict2instance(instance_type, fn)
end


@unpack distributions = input

map(distributions) do dist
    @unpack name, type = dist
    instance_type = eval(Symbol(type))
    workspace[name] = dict2instance(instance_type, distributions[1]; workspace)
end


let
    model = workspace["my_model_for_reaction_intensity"].model
    plot(masses(model), Base.Fix1(unpolarized_intensity, model))
end

let
    model = workspace["my_model_for_reaction_intensity"].model

    # get a random point in the phase space
    σs0 = Invariants(masses(model);
        σ1=0.7980703453578917,
        σ2=3.6486261122281745)

    # call intensity
    _I = unpolarized_intensity(model, σs0)

    # call the amplitude
    _A = amplitude(model, σs0, [1, 0, 0, 1])  # pars: model, mandelstam variables, helicity values

    @testset "Tests from the original package" begin
        # @test 
        @test _I ≈ 9345.853380852352
        # # @test 
        @test _A ≈ -45.1323269502508 + 54.85942516648639im
        # 
        @test model.chains[2].Xlineshape(σs0.σ2) ≈
              model.chains[2].Xlineshape(σs0.σ2) ≈
              -0.5636481410171861 + 0.13763637759224928im
        # 
        @test model.chains[21].Xlineshape(σs0.σ1) ≈
              model.chains[22].Xlineshape(σs0.σ1) ≈
              model.chains[23].Xlineshape(σs0.σ1) ≈
              model.chains[24].Xlineshape(σs0.σ1) ≈ 2.1687201455088894 + 23.58225917009096im
    end
end


@unpack misc = input
@unpack amplitude_model_checksums = misc
@unpack parameter_points = input

# map(amplitude_model_checksums) do check_point_info
let check_point_info = amplitude_model_checksums[1]
    @unpack name, value, distribution = check_point_info
    # 
    # pull distribution
    dist = workspace[distribution]

    # pull correct parameter point
    parameter_points_dict = array2dict(parameter_points, "name")
    parameter_point = parameter_points_dict[name]
    @unpack parameters = parameter_point
    # 
    # compute, compare
    _parameters = array2dict(parameters, "name"; apply=v -> v["value"])
    @assert value ≈ dist(_parameters) "Check-point validation failed with $distribution 🥕"
    return "🟢"
end

