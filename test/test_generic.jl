using ThreeBodyDecaysIO: parse_into_function
using ThreeBodyDecaysIO
using Parameters
using Test

e(x) = 1 / (0.77^2 - x - 1im * 0.77 * 0.15 * sin(x))
# 
expression = "1/(0.77^2 - x - i*0.77*0.15*sin(x))"

@test expression_argument(Meta.parse(expression)) == :x
@test expression_argument(Meta.parse(replace(expression, "x" => "m12sq"))) == :m12sq

@testset "Parse into a function works" begin
    f = parse_into_function(expression)
    @test f(1.1) == e(1.1)
    # 
    f = parse_into_function(expression * "*8")
    @test f(1.1) == e(1.1) * 8
    # 
    @test_throws "There are multiple symbols" parse_into_function(expression * "*y")
end

@testset "BreitWignerWidthExpLikeBugg as a string" begin

    K700_BuggBW_lineshape = """
    1/(0.824^2 - σ - i * 0.824 * 
        (σ - 0.23397706275638377) / (0.824^2 - 0.23397706275638377) * 0.478 * exp(-0.941060 * σ))
    """

    X_K700_BuggBW = parse_into_function(K700_BuggBW_lineshape)
    @test X_K700_BuggBW(1.1) ≈ -1.6748633697157216 + 1.083006788846748im


    K1430_BuggBW_lineshape = """
    1/(1.375^2 - σ - i * 1.375 *
        (σ - 0.23397706275638377) / (1.375^2 - 0.23397706275638377) * 0.190 * exp(-0.020981 * σ))
    """
    X_K1430_BuggBW = parse_into_function(K1430_BuggBW_lineshape)
    @test X_K1430_BuggBW(1.1) ≈ 1.2297831004988296 + 0.20758230110948936im
end


let
    dict = Dict(
        "type" => "generic_function",
        "expression" => "1/(0.77^2 - x - i*0.77*0.15*sin(x))"
    )
    @unpack type = dict
    instance_type = eval(Symbol(type))
    custom_function = dict2instance(instance_type, dict)
    @testset "generic_function works" begin
        @test custom_function(Dict("x" => 1.1)) == e(1.1)
    end
end

