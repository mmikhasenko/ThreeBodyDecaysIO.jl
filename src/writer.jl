
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


function trivial_lineshape_parser(Xlineshape)
    scattering = LittleDict()
    FF_production = LittleDict()
    FF_decay = LittleDict()
    appendix = LittleDict()
    (; scattering, FF_production, FF_decay), appendix
end


"""
    serializeToDict(chain::AbstractDecayChain;
        name::AbstractString="my_decay_chain",
        lineshape_parser::Function=trivial_lineshape_parser)

Writes a `DecayChain` model to a dictionary including `propagators`, and `vertices`, and `topology`.
Function `lineshape_parser` is called on `chain.Xlineshape` to split it into `scattering`, `FF_production`, and `FF_decay`.
See `trivial_lineshape_parser`.
"""
function serializeToDict(
    chain::AbstractDecayChain;
    name::AbstractString = "my_decay_chain",
    lineshape_parser::Function = trivial_lineshape_parser,
)
    @unpack k = chain
    i, j = ij_from_k(k)
    #
    appendix = Dict()
    # energy-dependence
    X, a = lineshape_parser(chain.Xlineshape)
    merge!(appendix, a)
    @unpack scattering, FF_production, FF_decay = X
    propagator =
        Dict(:spin => chain.two_j |> d2, :parametrization => scattering, :node => [i, j])
    propagators = [propagator]
    #
    H1, a = serializeToDict(chain.HRk)
    push!(H1, :node => [[i, j], k], :formfactor => FF_production), merge!(appendix, a)
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

function serializeToDict(
    tbs::ThreeBodySystem;
    particle_labels::NTuple{4,String} = ("A", "B", "C", "X"),
)
    #
    @unpack ms, two_js = tbs
    system_dict = Dict(
        :initial_state => Dict(
            :name => particle_labels[4],
            :mass => ms[4],
            :index => 0,
            :spin => d2(two_js[4]),
        ),
        :final_state => [
            Dict(
                :name => particle_labels[1],
                :mass => ms[1],
                :index => 1,
                :spin => d2(two_js[1]),
            ),
            Dict(
                :name => particle_labels[2],
                :mass => ms[2],
                :index => 2,
                :spin => d2(two_js[2]),
            ),
            Dict(
                :name => particle_labels[3],
                :mass => ms[3],
                :index => 3,
                :spin => d2(two_js[3]),
            ),
        ],
    )
    system_dict, Dict()
end


"""
serializeToDict::ThreeBodyDecay;
        lineshape_parser::Function=trivial_lineshape_parser,
        particle_labels::NTuple{4,String}=("A", "B", "C", "X"))

Writes a `ThreeBodyDecay` model to a dictionary. The argument `lineshape_parser` is passed to the chain-serialization function.
The argument `particle_labels` is passed to the kinematics serialization function.

## Arguments

- `model::ThreeBodyDecay`: The model to serialize.
- `lineshape_parser`: A function that takes a `Xlineshape` and returns
a special tuple of a type `((; scattering, FF_production, FF_decay), appendix)`. See `trivial_lineshape_parser`.
- `particle_labels`: a tuple with labels of the particles in the decay.
```
"""
function serializeToDict(
    model::ThreeBodyDecay;
    lineshape_parser::Function = trivial_lineshape_parser,
    particle_labels::NTuple{4,String} = ("A", "B", "C", "X"),
)
    #
    @unpack chains, names, couplings = model

    appendix = Dict()
    _kinematics, a = serializeToDict(model.chains[1].tbs; particle_labels)
    merge!(appendix, a)
    #
    _chains = map(zip(chains, names, couplings)) do (chain, name, coupling)
        dict, a = serializeToDict(chain; name, lineshape_parser)
        merge!(appendix, a)
        push!(dict, :weight => replace(string(coupling), "im" => "i"))
    end
    #
    i, j, k = ijk(first(chains).k)
    #
    model_dict = LittleDict{Symbol,Any}(
        :kinematics => _kinematics,
        :reference_topology => [[i, j], k],
        :chains => _chains,
    )
    model_dict, appendix
end


function add_hs3_fields(decay_description, appendix, model_name = "my_amplitude_model")

    k = topology2k(decay_description[:reference_topology])
    variable_groups = variablesToDict(k)
    variable_names = vcat(map(variable_groups) do v
        v[:mass_phi_costheta]
    end...)

    dict = OrderedDict(
        :distributions => [
            OrderedDict(
                :type => "HadronicUnpolarizedIntensity",
                :name => model_name,
                :decay_description => decay_description,
                :variables => variable_groups,
            ),
        ],
        :functions => [(v[:name] = k; v) for (k, v) in appendix],
        :domains => [
            OrderedDict(
                :name => "default",
                :type => "product_domain",
                :axes => [
                    OrderedDict(:name => name, :min => -1, :max => 1) for
                    name in variable_names
                ],
            ),
        ],
        :misc => Dict(),
        :parameter_points => Dict(),
    )
    return dict
end


function variablesToDict(k::Int)
    i, j = ij_from_k(k)
    ij_str = "$(i)$(j)"
    ij_k_str = "$(i)$(j)_$(k)"
    return [
        Dict(
            :node => [i, j],
            :mass_phi_costheta => ["m_" * ij_str, "phi_" * ij_str, "cos_theta_" * ij_str],
        ),
        Dict(
            :node => [[i, j], k],
            :mass_phi_costheta =>
                ["m_" * ij_k_str, "phi_" * ij_k_str, "cos_theta_" * ij_k_str],
        ),
    ]
end
