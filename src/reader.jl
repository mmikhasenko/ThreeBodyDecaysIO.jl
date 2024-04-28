
ThreeBodyDecays.x2(s::String) = Int(eval(Meta.parse(s)) * 2)

function dict2kinematics(dict)
    initial_state = dict["initial_state"]
    final_states = dict["final_state"]

    # Parse initial state
    initial_index = initial_state["index"]
    initial_name = initial_state["name"]
    initial_spin = initial_state["spin"]
    initial_mass = initial_state["mass"]

    # Parse final states
    final_indices = [state["index"] for state in final_states]
    final_names = [state["name"] for state in final_states]
    final_spins = [state["spin"] for state in final_states]
    final_masses = [state["mass"] for state in final_states]

    # Create DataFrame
    df = DataFrame(
        indices=[initial_index; final_indices],
        names=[initial_name; final_names],
        spins=[initial_spin; final_spins],
        masses=[initial_mass; final_masses]
    )

    sort!(df, :indices)
    N = size(df, 1) - 1
    first(df.indices) != 0 && error("particle 0 should be indicating the overall system")
    !all(df.indices .== 0:N) &&
        error("Expected indices 0,1,...,$N for the reaction with N final particles")
    N != 3 && error("ThreeBodyDecays.jl only works with three-body systems only")

    ms = ThreeBodyMasses(df.masses[2:end]...; m0=df.masses[1])
    two_spins = df.spins .|> x2
    two_js = ThreeBodySpins(two_spins[2:end]...; two_h0=two_spins[1])

    tbs = ThreeBodySystem(ms, two_js)
    return tbs
end

function dict2chain(dict; tbs, workspace=Dict())
    coupling = string2complex(dict["weight"])
    name = dict["name"]
    # 
    @unpack vertices, propagators = dict
    @unpack topology = dict
    # 
    k = topology2k(topology)
    @show topology, k
    i, j = ij_from_k(k)
    # spin of suchannel resonance
    resonance = first(propagators)
    @show resonance["node"], [i, j]
    @assert resonance["node"] == [i, j]
    # 
    spin = resonance["spin"]
    two_j = spin |> x2
    two_js = tbs.two_js
    # build two vertices manualy
    ind_ij = findfirst(v -> v["node"] == [i, j], vertices)
    vertex_ij = vertices[ind_ij]
    Hij = dict2recoupling(vertex_ij, (; two_j_fin=[two_js[i], two_js[j]], two_j_ini=two_j))
    # 
    ind_Rk = findfirst(v -> v["node"] == [[i, j], k], vertices)
    vertex_Rk = vertices[ind_Rk]
    HRk = dict2recoupling(vertex_Rk, (; two_j_fin=[two_j, two_js[k]], two_j_ini=two_js[4]))

    # build lineshape
    scattering = resonance["parametrization"]
    @show name, scattering, k
    bw = scattering isa Dict ? dict2lineshape(scattering) : workspace[scattering]
    X = bw
    if vertex_Rk["formfactor"] != ""
        FF_Rk = build_or_fetch(vertex_Rk["formfactor"], workspace)
        # single variablre function for the form factor
        @unpack ms = tbs
        mR = :m ∈ fieldnames(typeof(bw)) ? bw.m : error("Value for the resonance mass (:m) not found in the lineshape type, $(typeof(bw))")
        # for 0->Rk decay
        p(σ) = breakup(ms[4], sqrt(σ), ms[k])
        FF_Rk_svf = FF_Rk(p) * (1 / FF_Rk(p(mR^2)))
        X *= FF_Rk_svf
    end
    if vertex_ij["formfactor"] != ""
        FF_ij = build_or_fetch(vertex_ij["formfactor"], workspace)
        # single variablre function for the form factor
        @unpack ms = tbs
        mR = :m ∈ fieldnames(typeof(bw)) ? bw.m : error("Value for the resonance mass (:m) not found in the lineshape type, $(typeof(bw))")
        # for R->ij decay
        q(σ) = breakup(sqrt(σ), ms[i], ms[j])
        FF_ij_svf = FF_ij(q) * (1 / FF_ij(q(mR^2)))
        X *= FF_ij_svf
    end
    chain = DecayChain(; k, two_j, Xlineshape=X, Hij, HRk, tbs)
    (; coupling, chain, name)
end

function build_or_fetch(key, workspace)
    key isa Dict && return dict2lineshape(key)
    return workspace[key]
end

function dict2model(input)
    model_descrition = first(input["distributions"])
    # 
    @assert model_descrition["type"] == "hadronic_cross_section_unpolarized_dist"
    @unpack decay_description = model_descrition
    # 
    @unpack functions = input
    workspace = Dict{String,Any}()
    for fn in functions
        workspace[fn["name"]] = dict2lineshape(fn)
    end
    # 
    @unpack kinematics = decay_description
    tbs = dict2kinematics(kinematics)
    # 
    df = dict2chain.(decay_description["chains"]; tbs, workspace) |> DataFrame
    model = ThreeBodyDecay(Vector{Pair{String,Tuple{Complex,AbstractDecayChain}}}(df.name .=> zip(df.coupling, df.chain)))
    return model
end

function dict2recoupling(dict, properties)
    if dict["type"] == "ls"
        @unpack l, s = dict
        two_ja, two_jb = properties.two_j_fin
        @unpack two_j_ini = properties
        return RecouplingLS(; two_ls=(l, s) .|> x2, two_ja, two_jb, two_j=two_j_ini)
    end
    if dict["type"] == "helicity"
        @unpack helicities = dict
        two_λa = helicities[1] |> x2
        two_λb = helicities[2] |> x2
        return NoRecoupling(two_λa, two_λb)
    end
    if dict["type"] == "parity"
        @unpack helicities, parity_factor = dict
        two_λa = helicities[1] |> x2
        two_λb = helicities[2] |> x2
        return ParityRecoupling(two_λa, two_λb, parity_factor == "+")
    end
    error("Unknown type: $(dict["type"])")
end

