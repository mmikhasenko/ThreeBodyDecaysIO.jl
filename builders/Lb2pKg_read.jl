# Script builds the model from JSON file

using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using HadronicLineshapes
using JSON
using Parameters
using DataFrames
using Plots



# get the JSON content
json_content = open(joinpath(@__DIR__, "..", "models", "Lb2pKg.json")) do io
    JSON.parse(io)
end

input = copy(json_content)

# pull model description from JSON content
model_descrition = first(input["distributions"])

# make sure that we deal with three body hadronic decay
@assert model_descrition["type"] == "hadronic_cross_section_unpolarized_dist"
@unpack decay_description = model_descrition
@unpack reference_topology = decay_description

# updated_input = update2values(input, decay_description["appendix"])

# recursively vcat(x)... to flatten the topology structure
# Ex: [[1, 2], 3] -> [1, 2, 3]
flatten_topology(topology) =
    topology isa Array ? vcat(flatten_topology.(topology)...) : topology

# [TEST] the particles are labeled 1,2,3
@assert flatten_topology(reference_topology) |> sort == [1, 2, 3]

@unpack functions = input

# build functions from JSON array,
# add into dictionary with the name as key
workspace = Dict{String,Any}()
for fn in functions
    workspace[fn["name"]] = dict2lineshape(fn)
end
workspace



# 
@unpack kinematics = decay_description
tbs = dict2kinematics(kinematics)

df = dict2chain.(decay_description["chains"]; tbs, workspace) |> DataFrame
model = ThreeBodyDecay(Vector{Pair{String,Tuple{Complex,AbstractDecayChain}}}(df.name .=> zip(df.coupling, df.chain)))

dp = randomPoint(tbs)
unpolarized_intensity(model, dp.σs)

plot(masses(model), Base.Fix1(unpolarized_intensity, model); iσx=3, iσy=2)

# do the same manually
let
    ms = masses(model)
    # 
    σ3v = range(lims3(ms)[1], 2.5^2, 100)[2:end-1]
    σ2v = range(lims2(ms)..., 30)[2:end-1]
    # 
    _model = model
    f(σs) = Kibble(σs, ms^2) > 0 ? 0 : unpolarized_intensity(_model, σs)

    calv = [
        f(Invariants(ms; σ3, σ2)) for σ2 in σ2v, σ3 in σ3v]
    heatmap(σ3v, σ2v, calv)
end



