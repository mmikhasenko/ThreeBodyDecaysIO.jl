@with_kw struct BW
    Γ0::Float64
    m0::Float64
end
(bw::BW)(σ) = 1 / (bw.m0^2 - σ - 1im * bw.m0 * bw.Γ0)

function serializeToDict(x::BW)
    return Dict{String,Any}(
        "type" => "BreitWigner",
        "width" => x.Γ0,
        "mass" => x.m0,
    ), Dict()
end
