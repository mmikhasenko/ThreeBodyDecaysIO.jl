
function wrap2dict(H::RecouplingLS)
	type = "ls"
	l, s = H.two_ls .|> d2
	LittleDict{Symbol, Any}(pairs((; type, l, s)))
end

function wrap2dict(H::ParityRecoupling)
	type = "parity"
	helicities = [H.two_λa, H.two_λb] .|> d2
	parity_factor = ηηηphaseisplus ? '+' : '-'
	LittleDict{Symbol, Any}(pairs((; type, helicities, parity_factor)))
end

function wrap2dict(H::NoRecoupling)
	type = "parity"
	helicities = [H.two_λa, H.two_λb] .|> d2
	LittleDict{Symbol, Any}(pairs(; type, helicities))
end

function lineshape2dict(chain, name)
	lineshape_type = typeof(chain.Xlineshape)
	type = "$(name)_$(lineshape_type)"
	spin = chain.two_j |> d2
	LittleDict{Symbol, Any}(pairs((; type, spin)))
end

function wrap2dict(chain::AbstractDecayChain, name::AbstractString)
	k = 3
	i, j = ij_from_k(k)
	# 
	vertices = [
		push!(wrap2dict(chain.HRk), :node => [[i, j], k]),
		push!(wrap2dict(chain.Hij), :node => [i, j])]
	# 
	propagators = [
		push!(lineshape2dict(chain, name), :node => [i, j])]
	# 
	topology = [[i, j], k]
	Dict{Symbol, Any}(pairs((; vertices, propagators, topology)))
end

function wrap2dict(tbs::ThreeBodySystem)
	Dict{Symbol, Any}(
		:indices => (1, 2, 3, 0),
		:names => ("p", "K", "g", "Lb"),
		:masses => values(tbs.ms),
		:spins => values(tbs.two_js) .|> d2,
	)
end

function wrap2dict(model::ThreeBodyDecay)
	@unpack chains, names, couplings = model

	_kinematics = wrap2dict(model.chains[1].tbs)
	# 
	_chains = map(zip(chains, names, couplings)) do (chain, name, coupling)
		dict = wrap2dict(chain, name)
		push!(dict, :weight => replace(string(coupling), "im" => "i"))
	end
	#
	_lineshapes = Dict(map(zip(chains, names)) do (chain, name)
		lineshape = chain.Xlineshape
		type = typeof(lineshape)
		fields = fieldnames(type)
		dict = Dict(fields .=> getfield.(Ref(lineshape), fields))
		haskey(dict, :L) && pop!(dict, :L)
		"$(name)_$(type)" => dict
	end)
	#
	LittleDict{Symbol, Any}(
		:kinematics => _kinematics,
		:lineshapes => _lineshapes,
		:reference_topology => [[1, 2], 3],
		:chains => _chains)
end


function cascade_nodes(σs, ms; k)
	i, j = ij_from_k(k)
	σk = σs[k]
	_cosθij = cosθij(σs, ms^2; k)
	Dict(
		"kinematic_point" => [
			Dict(
				"node" => [[i, j], k],
				"phi" => 0.0,
				"theta" => 0.0,
				"mass" => ms.m0),
			Dict(
				"node" => [i, j],
				"phi" => 0.0,
				"theta" => _cosθij,
				"mass" => sqrt(σk)),
		])
end

function fully_flatten(x)
	_x = copy(x)
	i = 0
	while !(_x isa AbstractVector{<:Number})
		i += 1
		@show _x, i
		_x = collect(Iterators.flatten(_x))
		i > 10 && break
	end
	return _x
end


function topology2k(topology::AbstractArray)
	indices = fully_flatten(topology)
	length(indices) != 3 && error("Topology with more that three particles is not implemented")
	!(topology ∈ ([[1, 2], 3], [[2, 3], 1], [[3, 1], 2])) &&
		error("Only regular topologies [[1, 2], 3], [[2, 3], 1], [[3, 1], 2] are implememnted")
	# 
	return topology[2]
end



function validation_section(model, dpp, topology)
	k = topology2k(topology)
	ms = masses(model)
	@unpack σs, two_λs = dpp
	# 
	amplitudes = map(model.chains) do chain
		value = amplitude(chain, σs, two_λs)
		replace(string(value), "im" => "i")
	end
	# 
	_validation = Dict{Symbol, Any}(
		:kinematic_point => Dict(
			:masses_angles => cascade_nodes(σs, ms; k),
			:spin_projections => values(two_λs) .|> d2,
		),
		:chains_values => amplitudes,
	)
	return _validation
end

# 
1.1;
# Dict(
# 	"kinematic_point" => [
# 		Dict(
# 			"node" => [[[3, 1], 2], 4],
# 			"phi" => 0.0,
# 			"theta" => 0.0,
# 			"mass" => 5.6),
# 		Dict(
# 			"node" => [[3, 1], 2],
# 			"phi" => 0.0,
# 			"theta" => 0.3,
# 			"mass" => 4.4),
# 		Dict(
# 			"node" => [3, 1],
# 			"phi" => 0.0,
# 			"theta" => 0.3,
# 			"mass" => 1.2),
# 	])


