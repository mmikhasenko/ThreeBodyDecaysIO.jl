

@with_kw struct HadronicUnpolarizedIntensity{M,P,D}
    model::M
    reference_k::Int
    parameters::P
    mass_angles_cascade_names::D
end

function (dist::HadronicUnpolarizedIntensity)(pars)
    @unpack model, mass_angles_cascade_names, reference_k = dist
    mass_angles_values = mass_angles_replace(pars, mass_angles_cascade_names)
    angles, σs = angles_invariants(mass_angles_values, masses(model); k=reference_k)
    unpolarized_intensity(model, σs)
end

function dict2instance(::Type{HadronicUnpolarizedIntensity}, dict; workspace)
    @unpack parameters, variables, decay_description = dict
    model = dict2instance(ThreeBodyDecay, decay_description; workspace)
    # 
    @unpack reference_topology = decay_description
    reference_k = topology2k(reference_topology)
    # 
    mass_angles_cascade_names = extract_mass_angles_names(variables, reference_k)
    HadronicUnpolarizedIntensity(; model, reference_k, parameters, mass_angles_cascade_names)
end



# help functions

function extract_mass_angles_names(variables, reference_k)
    variables_dict = array2dict(variables, "node"; apply=v -> v["mass_angles"])
    i, j, k = ijk(reference_k)
    # 
    mass_angles_Rk = variables_dict[[[i, j], k]]
    mass_angles_ij = variables_dict[[i, j]]
    return (mass_angles_Rk, mass_angles_ij)
end

function mass_angles_replace(parameter_values, mass_angles_names)
    return map(mass_angles_names) do three
        map(name -> parameter_values[name], three)
    end
end