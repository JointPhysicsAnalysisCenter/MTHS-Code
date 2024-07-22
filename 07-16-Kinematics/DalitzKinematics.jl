### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 023422ce-3f85-11ef-0185-19a1e34afd9e
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate()
	
	using ThreeBodyDecays
	using LinearAlgebra
	using Plots
	using Plots.PlotMeasures: mm
	using PlutoUI
	using LaTeXStrings
	# 
	theme(:wong2, frame=:box, grid=false, minorticks=true,
	    guidefontvalign=:top, guidefonthalign=:right,
	    foreground_color_legend=nothing, lab="")
end

# ╔═╡ 4cd08b25-b024-49be-8919-36f2dd2af2ed
md"""
# Dalitz plot coordinates (16.07 / Tue)
"""

# ╔═╡ 5b45ab02-54de-41c4-98f9-917f4ac91373
ms = ThreeBodyMasses(1,4,2; m0=18.0);

# ╔═╡ 958e46bc-b801-4548-80a4-45b316c48434
md"""
## CMS kinematics
"""

# ╔═╡ c1eb13f3-e827-49cd-9b1c-bcfc27ef8920
function four_momenta(σs, ms)
	ps = [breakup_Rk(σs[k], ms; k) for k in 1:3]
	p1, p2, p3 = ps
	
	c12 = cosζ12_for0(σs, ms^2)
	c31 = cosζ31_for0(σs, ms^2)

	E1, E2, E3 = sqrt.((ms^2)[1:3] .+ ps.^2)
	# 
	p1v = [0,0,-p1, E1]
	p2v = [p2*sqrt(1-c12^2),0, -p2*c12, E2]
	p3v = [-p3*sqrt(1-c31^2),0, -p3*c31, E3]
	# 
	@assert all((p1v+p2v+p3v)[1:3] .+ 1 .≈ 1.0) "p1v+p2v+p3v is not 0"
	return p1v, p2v, p3v
end

# ╔═╡ d4e66a78-26cb-4613-ab05-6d6e258dd05e
function plot_momenta!(σs, ms; mc)
	ps = four_momenta(σs, ms);
	#
	for (i,(p,m)) in enumerate(zip(ps, collect(ms)[1:3]))
		plot!(sp=1, [0im,  dot(p, [1im, 0, 1, 0])], arrow=true,
			lab=latexstring("\$p_$i\\,(m=$(round(Int, m)))\$"))
	end
	scatter!(sp=1, [0im]; ms=10, mc)
	plot!(sp=1, xlab="", ylab="")
	# 
	annotate!(sp=1,
		[(reim(dot(p,[1im, 0, 1, 0])/2)..., 
				string(round(norm(p[1:3]), digits=1)), :bottom) for p in ps])
end

# ╔═╡ be341733-eaaa-4950-ad1a-36bd4e20b353
function plot_two(σs, ms; mc=2)
	# 
	plot(layout=grid(1,2), size=(700,350), aspect_ratio=1)
	plot_momenta!(σs, ms; mc)
	# 
	plot!(sp=2, border13(ms), lw=2)
	scatter!(sp=2, [σs.σ1], [σs.σ3], ms=10; mc)
	# 
	plot!(sp=2, xlab=L"m_{23}^2", ylab=L"m_{13}^2")
end

# ╔═╡ 75c498ee-9345-48dd-a389-73fec3edc5df
md"""
## Interactive exploration
"""

# ╔═╡ 05750d26-626b-4ba2-b56d-996bbb0608b5
@bind x Slider((0:0.01:1)[2:end-1], show_value=true, default=0.3)

# ╔═╡ 121d754a-182d-4783-a90d-d7aa25d488e2
@bind y Slider((0:0.01:1)[2:end-1], show_value=true, default=0.6)

# ╔═╡ faada165-0e32-4922-9fcd-d851e7feb09f
σ_test = x2σs([x,y], ms; k=1)

# ╔═╡ 5caaa628-ee75-41f7-af6c-d58465622be2
plot_two(σ_test, ms)

# ╔═╡ 173591df-abb8-4a96-be6a-b61f273a0339
md"""
## Examples
"""

# ╔═╡ 751c031f-a03e-4b8e-bd71-f0f32033f441
xy_array = [
	[0.9,0.9], [0.1,0.5],
	[0.3,0.1], [0.9,0.3],
	[0.7,0.7], [0.4,0.9],
	[0.3,0.9], [0.3,0.5]
	];

# ╔═╡ 06cc4774-d4c6-4295-89d3-80206c0a234a
is_saved = false

# ╔═╡ 6382f4cc-1781-4a7b-8f75-42f98584681b
for (mc, xy) in enumerate(xy_array)
	plot(aspect_ratio=1, size=(400,350))
	σs = x2σs(xy, ms; k=2)
	plot_momenta!(σs, ms; mc)
	is_saved && savefig("~/Downloads/four_momenta_$mc.pdf")
end

# ╔═╡ cd553fc8-3a72-4d55-b3a0-7f0a923cf4fa
let
	plot(border13(ms), lw=2, lab="", aspect_ratio=1)
	for (mc, xy) in enumerate(xy_array)
		σs = x2σs(xy, ms; k=2)
		scatter!([σs.σ1], [σs.σ3], ms=10; mc, lab="point $mc")
	end
	plot!()
	savefig("~/Downloads/four_momenta_solution.pdf")
	plot!(xlab=L"m_{23}^2", ylab=L"m_{13}^2")
end

# ╔═╡ Cell order:
# ╟─4cd08b25-b024-49be-8919-36f2dd2af2ed
# ╠═023422ce-3f85-11ef-0185-19a1e34afd9e
# ╠═5b45ab02-54de-41c4-98f9-917f4ac91373
# ╟─958e46bc-b801-4548-80a4-45b316c48434
# ╠═c1eb13f3-e827-49cd-9b1c-bcfc27ef8920
# ╠═d4e66a78-26cb-4613-ab05-6d6e258dd05e
# ╠═be341733-eaaa-4950-ad1a-36bd4e20b353
# ╟─75c498ee-9345-48dd-a389-73fec3edc5df
# ╠═5caaa628-ee75-41f7-af6c-d58465622be2
# ╠═05750d26-626b-4ba2-b56d-996bbb0608b5
# ╠═121d754a-182d-4783-a90d-d7aa25d488e2
# ╠═faada165-0e32-4922-9fcd-d851e7feb09f
# ╟─173591df-abb8-4a96-be6a-b61f273a0339
# ╠═751c031f-a03e-4b8e-bd71-f0f32033f441
# ╠═06cc4774-d4c6-4295-89d3-80206c0a234a
# ╠═6382f4cc-1781-4a7b-8f75-42f98584681b
# ╠═cd553fc8-3a72-4d55-b3a0-7f0a923cf4fa
