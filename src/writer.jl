
function serializeToDict(H::RecouplingLS)
    type = "ls"
    l, s = H.two_ls .|> d2
    H_dict = LittleDict{Symbol,Any}(pairs((; type, l, s)))
    (H_dict, Dict())
end

function serializeToDict(H::ParityRecoupling)
    type = "parity"
    helicities = [H.two_λa, H.two_λb] .|> d2
    parity_factor = H.ηηηphaseisplus ? '+' : '-'
    H_dict = LittleDict{Symbol,Any}(pairs((; type, helicities, parity_factor)))
    (H_dict, Dict())
end

function serializeToDict(H::NoRecoupling)
    type = "helicity"
    helicities = [H.two_λa, H.two_λb] .|> d2
    H_dict = LittleDict{Symbol,Any}(pairs((; type, helicities)))
    (H_dict, Dict())
end


function serializeToDict(chain::AbstractDecayChain, name::AbstractString, lineshape_parser)
    k = 3
    i, j = ij_from_k(k)
    # 
    appendix = Dict()
    # energy-dependence
    X, a = lineshape_parser(chain.Xlineshape)
    merge!(appendix, a)
    @unpack scattering, FF_production, FF_decay = X
    propagator = Dict(
        :spin => chain.two_j |> d2,
        :parametrization => scattering,
        :node => [i, j])
    propagators = [propagator]
    # 
    H1, a = serializeToDict(chain.HRk)
    push!(H1, :node => [[i, j], k], :formfactor => FF_production),
    merge!(appendix, a)
    # 
    H2, a = serializeToDict(chain.Hij)
    push!(H2, :node => [i, j], :formfactor => FF_decay)
    merge!(appendix, a)
    # 
    vertices = [H1, H2]
    # 
    topology = [[i, j], k]
    chain_dict = Dict{Symbol,Any}(pairs((; vertices, propagators, topology, name)))
    (chain_dict, appendix)
end

function serializeToDict(tbs::ThreeBodySystem)
    system_dict = Dict{Symbol,Any}(
        :indices => (1, 2, 3, 0),
        :names => ("A", "B", "C", "X"),
        :masses => values(tbs.ms),
        :spins => values(tbs.two_js) .|> d2,
    )
    system_dict, Dict()
end

function serializeToDict(model::ThreeBodyDecay, lineshape_parser)
    @unpack chains, names, couplings = model

    appendix = Dict()
    _kinematics, a = serializeToDict(model.chains[1].tbs)
    merge!(appendix, a)
    # 
    _chains = map(zip(chains, names, couplings)) do (chain, name, coupling)
        dict, a = serializeToDict(chain, name, lineshape_parser)
        merge!(appendix, a)
        push!(dict, :weight => replace(string(coupling), "im" => "i"))
    end
    #
    model_dict = LittleDict{Symbol,Any}(
        :kinematics => _kinematics,
        :reference_topology => [[1, 2], 3],
        :chains => _chains)
    model_dict, appendix
end


function topology2k(topology::AbstractArray)
    indices = fully_flatten(topology)
    length(indices) != 3 && error("Topology with more that three particles is not implemented")
    !(topology ∈ ([[1, 2], 3], [[2, 3], 1], [[3, 1], 2])) &&
        error("Only regular topologies [[1, 2], 3], [[2, 3], 1], [[3, 1], 2] are implememnted")
    # 
    return topology[2]
end
