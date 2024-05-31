using ThreeBodyDecaysIO
using MacroTools
using AssertionError
using Parameters


e(x) = 1 / (0.77^2 - x - 1im * 0.77 * 0.15 * sin(x))
# 
expression = "1/(0.77^2 - x - i*0.77*0.15*sin(x))"

@testset "Parse into a function works" begin
    f = ThreeBodyDecaysIO.parse_into_function(expression)
    @test f(1.1) == e(1.1)
    # 
    f = ThreeBodyDecaysIO.parse_into_function(expression * "*8")
    @test f(1.1) == e(1.1) * 8
    # 
    @test_throws "There are multiple symbols" ThreeBodyDecaysIO.parse_into_function(expression * "*y")
end

let
    dict = Dict(
        "type" => "GenericFunction",
        "expression" => "1/(0.77^2 - x - i*0.77*0.15*sin(x))"
    )
    @unpack type = dict
    instance_type = eval(Symbol(type))
    custom_function = dict2instance(instance_type, dict)
    @testset "GenericFunction works" begin
        @test custom_function(1.1) == e(1.1)
    end
end

