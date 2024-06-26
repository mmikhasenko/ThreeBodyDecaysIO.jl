
make_julia_complex(expr) = replace(string(expr), r"(?<![a-zA-Z])i(?![a-zA-Z])" => "im")

const operator_symbols = Set([
    :+, :-, :*, :/, :%, :^,
    :im, :i,
    :sin, :log, :cos, :tan, :asin, :acos, :atan,
    :sinh, :cosh, :tanh, :asinh, :acosh, :atanh,
    :exp, :sqrt, :abs, :real, :imag, :conj, :pow, :log10, :angle, :abs2,
    :sign, :round, :trunc, :floor, :ceil
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

parse_into_function(expression::String) = parse_into_function(Meta.parse(make_julia_complex(expression)))

function parse_into_function(expression::Expr)
    # Collect symbols from the parsed expression
    symbols = collect_nonregistered_symbols(expression)
    # Ensure there is only one symbol
    @assert length(symbols) == 1 "There are multiple symbols in the expression: $symbols"
    symbol = first(symbols)
    f = eval(Expr(:->, Symbol(symbol), expression))
    return f
end


struct generic_function end

function dict2instance(::Type{generic_function}, dict)
    @unpack expression = dict
    d = copy(dict)
    # 
    pop!(d, "type")
    pop!(d, "expression")
    updated_expression = replace(expression, (k => "($v)" for (k, v) in d)...)
    # 
    parsed_expression = Meta.parse(make_julia_complex(updated_expression))
    f = parse_into_function(parsed_expression)
    arg = expression_argument(parsed_expression)
    NamedArgFunc(WrapFlexFunction(f), [string(arg)])
end
