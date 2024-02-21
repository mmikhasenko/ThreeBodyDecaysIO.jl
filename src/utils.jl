function update2values(x, ref)
	if x isa String && haskey(ref, x)
		return ref[x]
	end
	if x isa Dict
		_d = Dict{String, Any}()
		for (k, v) in x
			_d[k] = update2values(v, ref)
		end
		return _d
	end
	if x isa AbstractArray
		return update2values.(x, Ref(ref))
	end
	return x
end

string2complex(s) = eval(Meta.parse(replace(s, "i" => "im")))
