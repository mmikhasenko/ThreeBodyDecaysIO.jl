# Script builds the model from JSON file

using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using HadronicLineshapes
using JSON
using Parameters
using DataFrames
using Plots



# get the JSON content
json_content = open(joinpath(@__DIR__, "..", "models", "Lc2ppiK.json")) do io
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

plot(masses(model), Base.Fix1(unpolarized_intensity, model))

# get a random point in the phase space
σs0 = Invariants(masses(model);
    σ1=0.7980703453578917,
    σ2=3.6486261122281745)

# call intensity
_I = unpolarized_intensity(model, σs0)

# call the amplitude
_A = amplitude(model, σs0, [1, 0, 0, 1])  # pars: model, mandelstam variables, helicity values

let
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

