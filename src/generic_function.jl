
make_julia_complex(expr) = replace(string(expr), r"(?<![a-zA-Z])i(?![a-zA-Z])" => "im")

const operator_symbols = Set([
    :+,
    :-,
    :*,
    :/,
    :%,
    :^,
    :im,
    :i,
    :sin,
    :log,
    :cos,
    :tan,
    :asin,
    :acos,
    :atan,
    :sinh,
    :cosh,
    :tanh,
    :asinh,
    :acosh,
    :atanh,
    :exp,
    :sqrt,
    :abs,
    :real,
    :imag,
    :conj,
    :pow,
    :log10,
    :angle,
    :abs2,
    :sign,
    :round,
    :trunc,
    :floor,
    :ceil,
])

# Step 2: Function to traverse the syntax tree and collect symbols using MacroTools
function collect_nonregistered_symbols(expr)
    symbols = Set{Symbol}()
    MacroTools.prewalk(expr) do node
        if isa(node, Symbol) && !(node in operator_symbols)
            push!(symbols, node)
        end
        node
    end
    return symbols
end

expression_arguments(expr) = collect_nonregistered_symbols(expr)
function expression_argument(expr)
    args = expression_arguments(expr)
    length(args) != 1 && error("The expression has multiple arguments: $args")
    return first(args)
end

function print_function(body, arg, pars)
    function_def_str = "$arg -> begin\n"
    map(Tuple(keys(pars))) do k
        function_def_str *= "$k = $(pars[k])\n"
    end
    function_def_str *= "$(body)\n end"
    #
    return function_def_str
end

function parse_into_function(body, pars = Dict("i" => 1im))
    all_symbols = expression_arguments(Meta.parse(body))
    variables = [s for s in all_symbols if !haskey(pars, string(s))]
    @assert length(variables) == 1 "Number of undefined variables more than one: variables=$variables"
    arg = variables[1]
    #
    function_def_str = print_function(body, arg, pars)
    function_def = Meta.parse(function_def_str)

    # Use RuntimeGeneratedFunction to avoid world age issues
    f = @RuntimeGeneratedFunction function_def
    f, arg
end

struct generic_function end

function dict2instance(::Type{generic_function}, dict)
    @unpack expression = dict
    parameters_dict = Dict{String,Any}(dict)
    pop!(parameters_dict, "expression")
    haskey(parameters_dict, "name") && pop!(parameters_dict, "name")
    haskey(parameters_dict, "type") && pop!(parameters_dict, "type")
    parameters_dict["i"] = 1im
    f, arg = parse_into_function(expression, parameters_dict)
    # arg = expression_argument(parsed_expression)
    return NamedArgFunc(WrapFlexFunction(f), [string(arg)])
end
