using Test

using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.JSON


model = let
    tbs = ThreeBodySystem(
        ms=ThreeBodyMasses(0.141, 0.142, 0.143; m0=3.09),
        two_js=ThreeBodySpins(0, 0, 0; two_h0=2))
    dpp = randomPoint(tbs)
    #
    two_j = 2
    ch1 = DecayChain(;
        k=1,
        two_j,
        Xlineshape=BW(4.1, 0.1),
        Hij=RecouplingLS(two_j, (two_j, 0), 0, 0),
        HRk=RecouplingLS(tbs.two_js[4], (two_j, two_j), two_j, 0),
        tbs)
    ch2 = DecayChain(ch1; k=2)
    ch3 = DecayChain(ch1; k=3)
    # 
    ThreeBodyDecay(
        "K(892)" .=> [(4.0, ch1), (2.0, ch2), (3.0, ch3)])
end

function lineshape_parser(Xlineshape)
    appendix = Dict()
    scattering, a = "X(4100)", Dict(
        "X(4100)" => Dict(
            :type => typeof(Xlineshape),
            :mass => Xlineshape.m0,
            :width => Xlineshape.Γ0))
    merge!(appendix, a)
    FF_production = "BlattWeisskopf(resonance)"
    FF_decay = "BlattWeisskopf(b-decay)"
    a = Dict(
        "BlattWeisskopf(resonance)" => Dict(
            :type => "BlattWeisskopf",
            :scale => 1.5
        ),
        "BlattWeisskopf(b-decay)" => Dict(
            :type => "BlattWeisskopf",
            :scale => 5.0
        )
    )
    merge!(appendix, a)
    (; scattering, FF_production, FF_decay), appendix
end

dict, appendix = wrap2dict(model, lineshape_parser)
dict[:appendix] = appendix
dict[:validation] = validation_section(
    model, randomPoint(model.chains[1].tbs), dict[:reference_topology])
# 
open("test.json", "w") do io
    JSON.print(io, dict, 4)
end

json_content = open("test.json") do io
    JSON.parse(io)
end

@testset "JSON has all key sections" begin
    @test haskey(json_content, "kinematics")
    @test haskey(json_content, "reference_topology")
    @test haskey(json_content, "chains")
    # 
    @test haskey(json_content, "appendix")
    @test haskey(json_content, "validation")
    # 
    @unpack chains = json_content
    chain = chains[1]
    @test haskey(chain, "vertices")
    @test haskey(chain, "propagators")
end

@testset "Specific to the test model" begin
    @test length(json_content["chains"]) == 3

    @unpack chains = json_content
    map(chains) do chain
        @test length(chain["vertices"]) == 2
        @test length(chain["propagators"]) == 1
    end
end



@testset "Parse kinematics" begin
    test_kinematics_dict = Dict{String,Any}(
        "names" => Any["p", "K", "g", "Lb"],
        "spins" => Any["1/2", "0", "1", "1/2"],
        "indices" => Any[1, 2, 3, 0],
        "masses" => Any[0.141, 0.49, 0.0, 5.2])
    # 
    test_tbs = dict2kinematics(test_kinematics_dict)
    # 
    test_tbs.ms == ThreeBodyMasses(0.141, 0.49, 0.0; m0=5.2)
    test_tbs.ms == ThreeBodySpins(1, 0, 2; two_h0=1)
end

input = copy(json_content)

@testset "Parse chain" begin
    updated_input = update2values(input, json_content["appendix"])
    tbs = dict2kinematics(updated_input["kinematics"])
    cdn = dict2chain(updated_input["chains"][1], tbs)
    @test cdn.chain isa DecayChain
    @test cdn.name isa AbstractString
    @test cdn.coupling isa Number
end

@testset "Parse model" begin
    @test dict2model(input) isa ThreeBodyDecay
end

