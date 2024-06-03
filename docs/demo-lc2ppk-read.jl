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
    @unpack mass, width, slope, x = dict
    bw = BreitWignerWidthExpLikeBugg(mass, width, slope)
    parameters = String[]
    variables = [x]
    return NamedArgFunc(bw, variables, parameters)
end


# code
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# get the JSON content
input = open(joinpath(@__DIR__, "..", "test", "lc2ppi-lhcb-test.json")) do io
    JSON.parse(io)
end

# built functions will be stored in workspace
workspace = Dict{String,Any}()

# build functions from JSON array,
@unpack functions = input
for fn in functions
    @unpack name, type = fn
    instance_type = eval(Symbol(type))
    workspace[name] = dict2instance(instance_type, fn)
end

# build distributions from JSON array,
@unpack distributions = input
for dist in distributions
    @unpack name, type = dist
    instance_type = eval(Symbol(type))
    workspace[name] = dict2instance(instance_type, distributions[1]; workspace)
end

# perform validation

@unpack misc, parameter_points = input
@unpack amplitude_model_checksums = misc

# map(amplitude_model_checksums) do check_point_info
let check_point_info = amplitude_model_checksums[1]
    @unpack name, value, distribution = check_point_info
    # 
    # pull distribution
    dist = workspace[distribution]

    # pull correct parameter point
    parameter_points_dict = array2dict(parameter_points; key="name")
    # find the point in the list of points
    parameter_point = parameter_points_dict[name]
    # compute, compare
    _parameters = array2dict(parameter_point["parameters"];
        key="name", apply=v -> v["value"])
    @assert value ≈ dist(_parameters) "Check-point validation failed with $distribution 🥕"
    return "🟢"
end


# plot the model
let
    model = workspace["my_model_for_reaction_intensity"].model
    plot(masses(model), Base.Fix1(unpolarized_intensity, model))
end
