using ThreeBodyDecaysIO
using Test

@testset "Reorder function" begin
    @test (1, 2, 3) |> reorder(3) == (1, 2, 3)
    @test (2, 3, 1) |> reorder(1) == (1, 2, 3)
    @test (3, 1, 2) |> reorder(2) == (1, 2, 3)
end