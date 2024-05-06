using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.HadronicLineshapes
using Test


@testset "BreitWigner from plane Dict" begin
    d = Dict("type" => "BreitWigner", "name" => "L1520_BW", "mass" => 1.0, "width" => 0.1, "ma" => 0.0, "mb" => 0.0, "l" => 0, "d" => 1.5)
    bw1 = dict2instance(BreitWigner, d)
    @test bw1 isa HadronicLineshapes.AbstractFlexFunc
end

@testset "BlattWeisskopf from plane Dict" begin
    d = Dict("type" => "BlattWeisskopf", "l" => 0, "radius" => 1.5)
    bw1 = dict2instance(BlattWeisskopf, d)
    @test bw1 isa HadronicLineshapes.BlattWeisskopf
end

@testset "MultichannelBreitWigner from a nasted Dict" begin
    d = Dict("type" => "MultichannelBreitWigner", "name" => "L1520_BW", "mass" => 1.0, "channels" => [Dict("gsq" => 0.1, "ma" => 0.0, "mb" => 0.0, "l" => 0, "d" => 1.5)])
    bw1 = dict2instance(MultichannelBreitWigner, d)
    @test bw1 isa HadronicLineshapes.AbstractFlexFunc
end
