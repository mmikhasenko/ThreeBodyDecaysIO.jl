using DataFrames


import ThreeBodyDecays: x2
x2(s::String) = Int(eval(Meta.parse(s)) * 2)

function parse_kinematics(dict)
	df = DataFrame(dict)
	sort!(df, :indices)
	N = size(df, 1) - 1
	first(df.indices) != 0 && error("particle 0 should be indicating the overall system")
	!all(df.indices .== 0:N) &&
		error("Expected intices 0,1,...,$N for the reaction with N final particles")
	N != 3 && error("ThreeBodyDecays.jl only works with three0body systems only")
	# 

	ms = ThreeBodyMasses(df.masses[2:end]...; m0 = df.masses[1])
	two_spins = df.spins .|> x2
	two_js = ThreeBodySpins(two_spins[2:end]...; two_h0 = two_spins[1])
	# 
	tbs = ThreeBodySystem(ms, two_js)
	return tbs
end


struct BW
	Γ0::Float64
	m0::Float64
end
(bw::BW)(σ) = 1 / (bw.m0^2 - σ - 1im * bw.m0 * bw.Γ0)


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

function dict2chain(dict, tbs)
	coupling = string2complex(dict["weight"])
	name = dict["name"]
	# 
	@unpack vertices, propagators = dict
	@unpack topology = dict
	# 
	k = topology2k(topology)
	i, j = ij_from_k(k)
	# spin of suchannel resonance
	resonance = first(propagators)
	@assert resonance["node"] == [i, j]
	# 
	spin = resonance["spin"]
	two_j = spin |> x2
	two_js = tbs.two_js
	# build two vertices manualy
	ind_ij = findfirst(v -> v["node"] == [i, j], vertices)
	@show vertex_ij = vertices[ind_ij]
	Hij = dict2H(vertex_ij, (; two_j_fin = [two_js[i], two_js[j]], two_j_ini = two_j))
	# 
	ind_Rk = findfirst(v -> v["node"] == [[i, j], k], vertices)
	vertex_Rk = vertices[ind_Rk]
	HRk = dict2H(vertex_Rk, (; two_j_fin = [two_j, two_js[k]], two_j_ini = two_js[4]))

	# build lineshape
	lineshape = resonance["lineshape"]
	X = dict2X(lineshape, ())
	chain = DecayChain(; k, two_j, Xlineshape = X, Hij, HRk, tbs)
	(; coupling, chain, name)
end

function dict2model(dict)
	tbs = parse_kinematics(dict["kinematics"])
	df = dict2chain.(dict["chains"], Ref(tbs)) |> DataFrame
	return ThreeBodyDecay(df.name .=> zip(df.coupling, df.chain))
end


function dict2X(dict, properties)
	if dict["type"] == "BW"
		return BW(; m0 = dict["mass"], Γ0 = dict["width"])
	end
	error("Unknown type: $(dict["type"])")
end


function dict2H(dict, properties)
	if dict["type"] == "ls"
		@unpack l, s = dict
		two_ja, two_jb = properties.two_j_fin
		@unpack two_j_ini = properties
		return RecouplingLS(; two_ls = (l, s) .|> x2, two_ja, two_jb, two_j = two_j_ini)
	end
	if dict["type"] == "helicity"
		@unpack lambda_a, lambda_b = dict
		two_λa = lambda_a |> x2
		two_λb = lambda_b |> x2
		return NoRecoupling(two_λa, two_λb)
	end
	if dict["type"] == "parity"
		@unpack lambda_a, lambda_b, parity = dict
		two_λa = lambda_a |> x2
		two_λb = lambda_b |> x2
		return ParityRecoupling(two)
	end
	error("Unknown type: $(dict["type"])")
end

