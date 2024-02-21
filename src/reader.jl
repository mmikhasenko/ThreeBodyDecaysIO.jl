using DataFrames

ThreeBodyDecays.x2(s::String) = Int(eval(Meta.parse(s)) * 2)

function dict2kinematics(dict)
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
	Hij = dict2recoupling(vertex_ij, (; two_j_fin = [two_js[i], two_js[j]], two_j_ini = two_j))
	# 
	ind_Rk = findfirst(v -> v["node"] == [[i, j], k], vertices)
	vertex_Rk = vertices[ind_Rk]
	HRk = dict2recoupling(vertex_Rk, (; two_j_fin = [two_j, two_js[k]], two_j_ini = two_js[4]))

	# build lineshape
	lineshape = resonance["lineshape"]
	X = dict2lineshape(lineshape, ())
	chain = DecayChain(; k, two_j, Xlineshape = X, Hij, HRk, tbs)
	(; coupling, chain, name)
end

function dict2model(dict)
	tbs = dict2kinematics(dict["kinematics"])
	df = dict2chain.(dict["chains"], Ref(tbs)) |> DataFrame
	return ThreeBodyDecay(df.name .=> zip(df.coupling, df.chain))
end


function dict2lineshape(dict, properties)
	if dict["type"] == "BW"
		return BW(; m0 = dict["mass"], Γ0 = dict["width"])
	end
	error("Unknown type: $(dict["type"])")
end


function dict2recoupling(dict, properties)
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

