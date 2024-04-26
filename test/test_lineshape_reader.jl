using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.HadronicLineshapes
using Test


@testset "HS3InputWrapper from nasted Dict" begin
    d = Dict("type" => "BreitWigner", "name" => "L1520_BW", "mass" => 1.0, "channels" => [Dict("gsq" => 0.1, "ma" => 0.0, "mb" => 0.0, "l" => 0, "d" => 1.5)])
    wrapper = HS3InputWrapper(d)
    @test wrapper.nt == (channels=[Dict("gsq" => 0.1, "ma" => 0.0, "mb" => 0.0, "l" => 0, "d" => 1.5)], mass=1.0)
    @test wrapper isa HS3InputWrapper{(:channels, :mass)}
    bw1 = BreitWigner(wrapper)
    @test bw1 isa HadronicLineshapes.AbstractFlexFunc
    bw2 = dict2lineshape(d)
    @test bw1.channels == bw2.channels
    @test bw1.m == bw2.m
end

@testset "HS3InputWrapper from plane Dict" begin
    d = Dict("type" => "BreitWigner", "name" => "L1520_BW", "mass" => 1.0, "width" => 0.1, "ma" => 0.0, "mb" => 0.0, "l" => 0, "d" => 1.5)
    wrapper = HS3InputWrapper(d)
    @test wrapper.nt == (width=0.1, mass=1.0, ma=0.0, mb=0.0, d=1.5, l=0)
    @test wrapper isa HS3InputWrapper{(:width, :mass, :ma, :mb, :d, :l)}
    bw1 = BreitWigner(wrapper)
    @test bw1 isa HadronicLineshapes.AbstractFlexFunc
    bw2 = dict2lineshape(d)
    @test bw1.channels == bw2.channels
    @test bw1.m == bw2.m
end

