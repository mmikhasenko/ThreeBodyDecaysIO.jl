using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.JSON
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.DataFrames
using Test

list_of_models = readdir(joinpath(@__DIR__, "..", "models"))

map(list_of_models) do file_name
    @info "⭐ Reading model from $file_name ⭐"
    @testset "Test of $file_name" begin
        # get the JSON content
        input = open(joinpath(@__DIR__, "..", "models", file_name)) do io
            JSON.parse(io)
        end
        model_descrition = first(input["distributions"])
        @test model_descrition["type"] == "hadronic_cross_section_unpolarized_dist"

        @test haskey(model_descrition, "decay_description")
        @unpack decay_description = model_descrition

        @test haskey(decay_description, "reference_topology")
        @unpack reference_topology = decay_description

        @test flatten_topology(reference_topology) |> sort == [1, 2, 3]

        @info "🐰 Parcing kinamatics 🐰"
        @unpack kinematics = decay_description
        tbs = dict2kinematics(kinematics)


        @test haskey(input, "functions")
        @unpack functions = input

        @info "🔥 Building lineshapes functions 🔥"
        workspace = Dict{String,Any}()
        for fn in functions
            workspace[fn["name"]] = dict2lineshape(fn)
        end
        workspace

        @info "🦊 Establishing decay chains 🦊"
        df = dict2chain.(decay_description["chains"]; tbs, workspace) |> DataFrame
        model = ThreeBodyDecay(Vector{Pair{String,Tuple{Complex,AbstractDecayChain}}}(df.name .=> zip(df.coupling, df.chain)))

        @test model isa ThreeBodyDecay
    end
end



