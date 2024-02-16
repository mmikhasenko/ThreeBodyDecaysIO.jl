
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

function valudation_section(model, dpp)
	@unpack σs, two_λs = dpp
	amplitudes = amplitude.(model.chains, Ref(σs), Ref(two_λs))
	# 
	_validation = Dict{Symbol, Any}(
		:kinematic_point => Dict(
			:mandelstam_variables => Dict(pairs(σs)),
			:spin_projections => values(two_λs) .|> d2,
		),
		:chains_values => amplitudes,
	)
	return _validation
end
