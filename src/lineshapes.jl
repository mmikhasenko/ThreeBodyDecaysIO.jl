@with_kw struct BW
    m0::Float64
    Γ0::Float64
end
(bw::BW)(σ) = 1 / (bw.m0^2 - σ - 1im * bw.m0 * bw.Γ0)

function serializeToDict(x::BW)
    return Dict{Symbol,Any}(
        :type => "BreitWigner",
        :width => x.Γ0,
        :mass => x.m0,
        :ma => 0.0,
        :mb => 0.0,
        :l => 0,
        :d => 0.0
    ), Dict()
end
