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
