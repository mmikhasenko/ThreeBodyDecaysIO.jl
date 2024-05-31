
make_julia_complex(expr) = replace(string(expr), r"(?<![a-zA-Z])i(?![a-zA-Z])" => "im")

const operator_symbols = Set([
    :+, :-, :*, :/, :%, :^, :im,
    :sin, :log, :cos, :tan, :asin, :acos, :atan,
    :sinh, :cosh, :tanh, :asinh, :acosh, :atanh,
    :exp, :sqrt, :abs, :real, :imag, :conj, :pow, :log10, :angle, :abs2,
    :sign, :round, :trunc, :floor, :ceil
])

# Step 2: Function to traverse the syntax tree and collect symbols using MacroTools
function collect_symbols(expr)
    symbols = Set{Symbol}()
    MacroTools.prewalk(expr) do node
        if isa(node, Symbol) && !(node in operator_symbols)
            push!(symbols, node)
        end
        node
    end
    return symbols
end


function parse_into_function(expression)
    # Parse the expression
    parsed_expr = Meta.parse(make_julia_complex(expression))
    # Collect symbols from the parsed expression
    symbols = collect_symbols(parsed_expr)
    # Ensure there is only one symbol
    @assert length(symbols) == 1 "There are multiple symbols in the expression: $symbols"
    symbol = first(symbols)
    f = eval(Expr(:->, Symbol(symbol), parsed_expr))
    return f
end


struct GenericFunction end

function dict2instance(::Type{GenericFunction}, dict)
    @unpack expression = dict
    f = parse_into_function(expression)
    WrapFlexFunction(f)
end
