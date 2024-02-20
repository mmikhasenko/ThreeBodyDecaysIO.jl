using Test

using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.JSON

struct BW
	Γ0::Float64
	m0::Float64
end
(bw::BW)(σ) = 1 / (bw.m0^2 - σ - 1im * bw.m0 * bw.Γ0)

model = let
	tbs = ThreeBodySystem(
		ms = ThreeBodyMasses(0.141, 0.142, 0.143; m0 = 3.09),
		two_js = ThreeBodySpins(0, 0, 0; two_h0 = 2))
	dpp = randomPoint(tbs)
	#
	two_j = 2
	ch1 = DecayChain(;
		k = 1,
		two_j,
		Xlineshape = BW(4.1, 0.1),
		Hij = RecouplingLS(two_j, (two_j, 0), 0, 0),
		HRk = RecouplingLS(tbs.two_js[4], (two_j, two_j), two_j, 0),
		tbs)
	ch2 = DecayChain(ch1; k = 2)
	ch3 = DecayChain(ch1; k = 3)
	# 
	ThreeBodyDecay(
		"K(892)" .=> [(4.0, ch1), (2.0, ch2), (3.0, ch3)])
end

dict = wrap2dict(model)
dict[:validation] = validation_section(
	model, randomPoint(model.chains[1].tbs), dict[:reference_topology])

json_content = let
	buffer = IOBuffer()
	JSON.print(buffer, dict, 4)
	s = String(take!(buffer))
	JSON.parse(s)
end

@testset "JSON has all key sections" begin
	@test haskey(json_content, "kinematics")
	@test haskey(json_content, "lineshapes")
	@test haskey(json_content, "reference_topology")
	@test haskey(json_content, "chains")
	@test haskey(json_content, "validation")
	# 
	@unpack chains = json_content
	chain = chains[1]
	@test haskey(chain, "vertices")
	@test haskey(chain, "propagators")
end

@testset "Specific to the test model" begin
	@test length(json_content["chains"]) == 3
	@test length(json_content["lineshapes"]) == 1

	@unpack chains = json_content
	map(chains) do chain
		@test length(chain["vertices"]) == 2
		@test length(chain["propagators"]) == 1
	end
end



@testset "Parse kinematics" begin
	test_kinematics_dict = Dict{String, Any}(
		"names"   => Any["p", "K", "g", "Lb"],
		"spins"   => Any["1/2", "0", "1", "1/2"],
		"indices" => Any[1, 2, 3, 0],
		"masses"  => Any[0.141, 0.49, 0.0, 5.2])
	# 
	test_tbs = parse_kinematics(test_kinematics_dict)
	# 
	test_tbs.ms == ThreeBodyMasses(0.141, 0.49, 0.0; m0 = 5.2)
	test_tbs.ms == ThreeBodySpins(1, 0, 2; two_h0 = 1)
end