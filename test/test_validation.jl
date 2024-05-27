using Test
using ThreeBodyDecaysIO
using ThreeBodyDecays

@testset "Order in angles_invariants" begin
    ms = ThreeBodyMasses(1.0, 2.0, 1.3; m0=8.0)
    σs = x2σs([0.1, 0.2], ms; k=2)

    _cosθ12 = cosθ12(σs, ms^2)
    mass_angles_cascade = [(ms.m0, 0.0, 1.0), (sqrt(σs.σ3), 0, _cosθ12)]
    angles, _σs = angles_invariants(mass_angles_cascade, ms; k=3)
    @test all(collect(_σs) .≈ collect(σs))

    _cosθ23 = cosθ23(σs, ms^2)
    mass_angles_cascade = [(ms.m0, 0.0, 1.0), (sqrt(σs.σ1), 0, _cosθ23)]
    angles, _σs = angles_invariants(mass_angles_cascade, ms; k=1)
    @test all(collect(_σs) .≈ collect(σs))

    _cosθ31 = cosθ31(σs, ms^2)
    mass_angles_cascade = [(ms.m0, 0.0, 1.0), (sqrt(σs.σ2), 0, _cosθ31)]
    angles, _σs = angles_invariants(mass_angles_cascade, ms; k=2)
    @test all(collect(_σs) .≈ collect(σs))
end




model = let
    tbs = ThreeBodySystem(
        ms=ThreeBodyMasses(0.141, 0.142, 0.143; m0=3.09),
        two_js=ThreeBodySpins(0, 0, 0; two_h0=2))
    dpp = randomPoint(tbs)
    #
    two_j = 2
    ch1 = DecayChain(;
        k=1,
        two_j,
        Xlineshape=BW(1.1, 0.1),
        Hij=RecouplingLS((two_j, 0)),
        HRk=RecouplingLS((two_j, two_j)),
        tbs)
    ch2 = DecayChain(ch1; k=2)
    ch3 = DecayChain(ch1; k=3)
    # 
    ThreeBodyDecay(
        "K892" .=> [(4.0, ch1), (2.0, ch2), (3.0, ch3)])
end

@testset "Order and values in validation_fields" begin
    ms = masses(model)
    σs = x2σs([0.2, 0.3], ms; k=3)
    model_name, point_name = "mymodel", "mypoint"
    _fields1 = validation_fields(model, σs; k=3, point_name, model_name)
    _fields2 = validation_fields(model, (m=sqrt(σs[3]), cosθ=-0.4); k=3, point_name, model_name)
    @test _fields1[2]["value"] ≈ _fields2[2]["value"]
    @test _fields1[2]["point"] == _fields2[2]["point"] == point_name
    @test _fields1[2]["distribution"] == _fields2[2]["distribution"] == model_name
    d1 = array2dict(_fields1[1]["parameters"]; key="name", apply=x -> x["value"])
    d2 = array2dict(_fields2[1]["parameters"]; key="name", apply=x -> x["value"])
    @test collect(keys(d1)) == collect(keys(d2)) == ["m_12_3", "phi_12_3", "cos_theta_12_3", "m_12", "phi_12", "cos_theta_12"]
    @test all(values(d1) .≈ values(d2))
end