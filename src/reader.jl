using DataFrames


import ThreeBodyDecays: x2
x2(s::String) = Int(eval(Meta.parse(s)) * 2)

function parse_kinematics(dict)
	df = DataFrame(dict)
	sort!(df, :indices)
	N = size(df, 1) - 1
	first(df.indices) != 0 && error("particle 0 should be indicating the overall system")
	!all(df.indices .== 0:N) &&
		error("Expected intices 0,1,...,$N for the reaction with N final particles")
	N != 3 && error("ThreeBodyDecays.jl only works with three0body systems only")
	# 

	ms = ThreeBodyMasses(df.masses[2:end]...; m0 = df.masses[1])
	two_spins = df.spins .|> x2
	two_js = ThreeBodySpins(two_spins[2:end]...; two_h0 = two_spins[1])
	# 
	tbs = ThreeBodySystem(ms, two_js)
	return tbs
end


struct BW
	Γ0::Float64
	m0::Float64
end
(bw::BW)(σ) = 1 / (bw.m0^2 - σ - 1im * bw.m0 * bw.Γ0)

json_content = let
	buffer = IOBuffer()
	JSON.print(buffer, dict, 4)
	s = String(take!(buffer))
	JSON.parse(s)
end

