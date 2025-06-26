### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 3a3ad176-73a2-11ef-3d28-a959d3f2b6da
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
		PackageSpec(url="https://github.com/mmikhasenko/ThreeBodyDecaysIO.jl"),
		PackageSpec("Plots")
	])
	# 
	using ThreeBodyDecaysIO
	using ThreeBodyDecaysIO.ThreeBodyDecays
	using ThreeBodyDecaysIO.HadronicLineshapes
	using ThreeBodyDecaysIO.Parameters
	using ThreeBodyDecaysIO.JSON
	# 
	using Plots
end

# ╔═╡ a8c65293-d012-411b-b3e5-d8c313ab4ad0
theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    xlim=(:auto, :auto), ylim=(:auto, :auto),
    lab="")

# ╔═╡ 51129209-08f2-4409-a111-77b487f1e5e4
json_content = 
	open(joinpath(@__DIR__, "..", "..", "AmplitudeModelSerialization", "models", "lc2ppik-lhcb-2683025.json")) do io
	JSON.parse(io)
end

# ╔═╡ 82630a3d-a98a-436b-ab03-d67b7d842bf0
begin
	workspace = Dict{String,Any}()
	
	@unpack functions = json_content
	for fn in functions
	    @unpack name, type = fn
	    instance_type = eval(Symbol(type))
	    workspace[name] = dict2instance(instance_type, fn)
	end
	
	@unpack distributions = json_content
	for dist in distributions
	    @unpack name, type = dist
	    instance_type = eval(Symbol(type))
	    workspace[name] = dict2instance(instance_type, distributions[1]; workspace)
	end
end

# ╔═╡ b8d5cb23-10ec-4aa1-b560-bfa060771fd0
@unpack model = workspace["default_model"]

# ╔═╡ bc754a5a-5d3d-4f5f-9dd7-abb61ae041a4
σs_test = randomPoint(masses(model))

# ╔═╡ e7579201-69ac-4951-889c-1c39932fac23
plot(masses(model), σs->unpolarized_intensity(model, σs))

# ╔═╡ bf39fc32-2b77-487f-8400-3fa49c75e1e8
f(x,y) = unpolarized_intensity(model, x2σs([x,y], masses(model); k=1))

# ╔═╡ 63be27a1-a319-476f-9f27-7735d2c2f49c
heatmap(
	0:0.01:1 |> shift_by_half,
	0:0.01:1 |> shift_by_half,
	f)

# ╔═╡ 70b12858-8026-4958-9cf2-af125c01cbbb
begin
	struct binned2dDensity
	    grid::Matrix{Array{Float64,1}}
	    cumarr::Array{Float64}
	    density
	end
	# 
	arg1_lims(binned2dDensity) = (binned2dDensity.grid[1, 1][1], binned2dDensity.grid[end, end][1])
	arg2_lims(binned2dDensity) = (binned2dDensity.grid[1, 1][2], binned2dDensity.grid[end, end][2])
	# 
	function getbinned2dDensity(g, arg1_lims, arg2_lims, N1, N2)
	    grid = hcat([[[x1, x2]
	                  for x2 = range(arg2_lims[1], arg2_lims[2], length=N1)]
	                 for x1 = range(arg1_lims[1], arg1_lims[2], length=N2)]...)
	    #
	    weights = [g(v...) for v in (grid[2:end, 2:end] .+ grid[1:end-1, 1:end-1]) ./ 2]
	    # normalize by the cell size
	    for i = 1:size(grid, 1)-1, j = 1:size(grid, 2)-1
	        weights[i, j] *= (grid[i+1, j][2] - grid[i, j][2]) * (grid[i, j+1][1] - grid[i, j][1])
	    end
	    weights ./= sum(weights)
	    cumarr = [0; cumsum(vcat(weights...), dims=1)]
	    return binned2dDensity(grid, cumarr, g)
	end

	function Base.rand(bD::binned2dDensity, n::Int)
	    # get
	    binind = findfirst(bD.cumarr .> rand()) - 2
	    # pars back
	    Nrows, Ncols = size(bD.grid)
	    indCol = div(binind, Nrows - 1) + 1
	    indRow = mod(binind, Nrows - 1) + 1
	    indRow + 1 > Nrows && error("indRow+1 > Nrows: binind = $binind")
	    indCol + 1 > Ncols && error("indCol+1 > Ncols: binind = $binind")
	    s1, σ1 = bD.grid[indRow, indCol]
	    s2, σ2 = bD.grid[indRow+1, indCol+1]
	    #         s1 != s2 && error("Something is wrong!")
	    s = s1 + rand() * (s2 - s1)
	    σ = σ1 + rand() * (σ2 - σ1)
	    bD.density(s, σ) == 0.0 && return rand(bD)
	    return [s, σ]
	end
end

# ╔═╡ bb57df63-0f58-4f4e-83f7-cf540c1bbfd1
const f_dens_0 = getbinned2dDensity(f, (0,1), (0,1), 100, 100);

# ╔═╡ 08a21b31-de90-46e4-ac08-f3fc692c3f55
heatmap(reshape(diff(f_dens_0.cumarr), size(f_dens_0.grid) .- 1))

# ╔═╡ 8f9edfde-8352-4e3e-908d-e669e0d6c6cb
function f_disc(f_dens, x,y)
	Δx = f_dens.grid[1,2][1] - f_dens.grid[1,1][1]
	Δy = f_dens.grid[2,1][2] - f_dens.grid[1,1][2]
	nx, ny = div(x, Δx) |> Int, div(y, Δy) |> Int
	# 
	function_values = reshape(diff(f_dens.cumarr), size(f_dens.grid) .- 1)
	function_values[ny+1,nx+1]
end

# ╔═╡ 200efdee-817f-492a-8817-8d272aa95e34
heatmap(
	0:0.01:1 |> shift_by_half,
	0:0.01:1 |> shift_by_half,	
	(x,y) -> f_disc(f_dens_0, x,y))

# ╔═╡ c347c83a-996e-4918-ba4a-714a0a351f2e
function efficiency(weights)
	ϵ = sum(weights)^2 / sum(weights.^2) / length(weights)
end

# ╔═╡ 032393a8-68d4-4c40-bf0a-7760a502faef
function efficiency(p,q; sampler, n = 1000)
	data = [rand(sampler, 2) for N in 1:1000]
	true_weights = map(x->p(x...), data)
	q_weights = map(x->q(x...), data)
	ws = true_weights./ q_weights
	efficiency(ws)
end

# ╔═╡ 640809c5-95dd-4176-b25b-ce482861382e
q(x,y) = f_disc(f_dens_0, x,y)

# ╔═╡ bce511e2-4eba-493e-ba00-e4b3ae916fd9
efficiency(f, q; sampler=f_dens_0)

# ╔═╡ b53eceea-c2a3-4d8e-8476-d5563ba2b2fb
efficiency(f, q; sampler=Float64)

# ╔═╡ b9ac9bb6-dd17-462c-9136-479d3e2ec5db
experiment = let
	Ns = 6:5:70
	ϵs = map(Ns) do n
		N = round(Int, n)
		f_dens = getbinned2dDensity(f, (0,1), (0,1), N, N)
		q(x,y) = f_disc(f_dens, x, y)
		efficiency(f, q; sampler=Float64, n=200)
	end
	(; Ns, ϵs)
end

# ╔═╡ 3099feec-e86f-42ed-a423-ae7d2b5842b3
plot(experiment.Ns.^2, experiment.ϵs, m=(4,:o),
	xlab="N", ylab="eff", ylims=(0.75, 1.001), xscale=:log10,
	lab="2d binning")

# ╔═╡ Cell order:
# ╠═3a3ad176-73a2-11ef-3d28-a959d3f2b6da
# ╠═a8c65293-d012-411b-b3e5-d8c313ab4ad0
# ╠═51129209-08f2-4409-a111-77b487f1e5e4
# ╠═82630a3d-a98a-436b-ab03-d67b7d842bf0
# ╠═b8d5cb23-10ec-4aa1-b560-bfa060771fd0
# ╠═bc754a5a-5d3d-4f5f-9dd7-abb61ae041a4
# ╠═e7579201-69ac-4951-889c-1c39932fac23
# ╠═bf39fc32-2b77-487f-8400-3fa49c75e1e8
# ╠═63be27a1-a319-476f-9f27-7735d2c2f49c
# ╠═70b12858-8026-4958-9cf2-af125c01cbbb
# ╠═bb57df63-0f58-4f4e-83f7-cf540c1bbfd1
# ╠═08a21b31-de90-46e4-ac08-f3fc692c3f55
# ╠═8f9edfde-8352-4e3e-908d-e669e0d6c6cb
# ╠═200efdee-817f-492a-8817-8d272aa95e34
# ╠═c347c83a-996e-4918-ba4a-714a0a351f2e
# ╠═032393a8-68d4-4c40-bf0a-7760a502faef
# ╠═640809c5-95dd-4176-b25b-ce482861382e
# ╠═bce511e2-4eba-493e-ba00-e4b3ae916fd9
# ╠═b53eceea-c2a3-4d8e-8476-d5563ba2b2fb
# ╠═b9ac9bb6-dd17-462c-9136-479d3e2ec5db
# ╠═3099feec-e86f-42ed-a423-ae7d2b5842b3
