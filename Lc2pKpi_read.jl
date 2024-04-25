# Script builds the model from JSON file

using ThreeBodyDecaysIO
using JSON
using Graphs
using GraphPlot
using Parameters
using HadronicLineshapes


struct DecayTopology{G}
    topology::G
    names::Dict{Any,Int}
end

node2string(n) = replace(string(n), "Any" => "")

ordered_node_names(dt::DecayTopology) =
    first.(sort(collect(dt.names), by=x -> x[2]))

# build graph from bracket representation,
# Ex: [[1,2],3]
# 
#       0
#       |
#   [[1,2],3]
#      /  \
#     /    3
#   [1,2]
#   /  \
#  /    2
# 1
function parse_topology!(g::DiGraph, node_map::Dict{Any,Int}, topology, parent::Int=0)
    # This function recursively processes the topology
    # If it's an array, it needs to create a new vertex for this node
    if topology isa Array
        # Create a new node for the current structure
        add_vertex!(g)
        current_node = nv(g)
        node_map[topology] = current_node

        # Link it to its parent if there is one
        if parent != 0
            add_edge!(g, parent, current_node)
        end

        # Process each child
        for item in topology
            parse_topology!(g, node_map, item, current_node)
        end
    else
        # It's a leaf node, check if it exists, if not, add it
        if !haskey(node_map, topology)
            add_vertex!(g)
            node_map[topology] = nv(g)
        end
        # Link leaf to the current parent node
        if parent != 0
            add_edge!(g, parent, node_map[topology])
        end
    end
end

function DecayTopology(topology::Array{Any,1})
    g = DiGraph()
    node_map = Dict{Any,Int}()

    parse_topology!(g, node_map, topology)
    DecayTopology(g, node_map)
end

using Test
# test
dt = DecayTopology([[1, 2], 3])
@test Set(node2string.(keys(dt.names))) == Set(["[1, 2]", "1", "2", "3", "[[1, 2], 3]"])
# compute adjacency matrix
@test adjacency_matrix(g) ==
      [0 1 0 0 1
    0 0 1 1 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0]

dt = DecayTopology([[2, 1], 3])
@test Set(node2string.(keys(dt.names))) == Set(["[2, 1]", "1", "2", "3", "[[2, 1], 3]"])
@test adjacency_matrix(g) ==
      [0 1 0 0 1
    0 0 1 1 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0]

let
    names = ordered_node_names(dt)
    nodelabel = node2string.(names)
    gplot(dt.topology; nodelabel)
end

dt = DecayTopology([[[1, 2], 3], 4])
gplot(g; nodelabel=labels(dt))

g = dt.topology
neighbors(g, 3)

# # get children of a node
# function decay_products(dt::DecayTopology, node_name)
#     index = dt.names[node_name]
#     nei = neighbors(g, index)
#     invnames = Dict(v => k for (k, v) in dt.names)
#     getindex.(Ref(invnames), collect(nei))
# end

# decay_products(dt::DecayTopology, [2, 1])
# decay_products(dt::DecayTopology, [[2, 1], 3]) == [2, 1], 3

# get the JSON content
json_content = open("Lc2ppiK.json") do io
    JSON.parse(io)
end

input = copy(json_content)
updated_input = update2values(input, decay_description["appendix"])

# pull model description from JSON content
model_descrition = first(updated_input["distributions"])

# make sure that we deal with three body hadronic decay
@assert model_descrition["type"] == "hadronic_cross_section_unpolarized_dist"
@unpack decay_description = model_descrition
@unpack reference_topology = decay_description

# recursively vcat(x)... to flatten the topology structure
# Ex: [[1, 2], 3] -> [1, 2, 3]
flatten_topology(topology) =
    topology isa Array ? vcat(flatten_topology.(topology)...) : topology

# [TEST] the particles are labeled 1,2,3
@assert flatten_topology(reference_topology) |> sort == [1, 2, 3]




@unpack reference_topology = decay_description


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

collect(workspace)[1][2]
using ThreeBodyDecays
using DataFrames
using Plots

df = dict2chain.(decay_description["chains"]; tbs, workspace) |> DataFrame
model = ThreeBodyDecay(df.name .=> zip(df.coupling, df.chain))

dp = randomPoint(tbs)
unpolarized_intensity(model, dp.σs)


L_BW = [v for (k, v) in workspace if contains(k, r"L.*BW")]
map(L_BW) do f
    f(dp.σs[2])
end

K_BW = [v for (k, v) in workspace if contains(k, r"K.*BW")]
map(K_BW) do f
    f(dp.σs[1])
end

D_BW = [v for (k, v) in workspace if contains(k, r"D.*BW")]
map(D_BW) do f
    f(dp.σs[3])
end


plot(masses(model), Base.Fix1(unpolarized_intensity, model))
