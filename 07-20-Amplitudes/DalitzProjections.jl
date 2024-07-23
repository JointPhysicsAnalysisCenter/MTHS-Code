### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ f787599e-45c8-11ef-2dd4-61b67a34bcbd
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(".")
    #
    using Plots
    using Plots.PlotMeasures: mm
    using Parameters
    using JSON
    using ThreeBodyDecaysIO
    using ThreeBodyDecaysIO.ThreeBodyDecays
    using ThreeBodyDecaysIO.HadronicLineshapes
    #
    using QuadGK
    using ThreeBodyDecaysIO
    #
    theme(
        :wong2,
        frame = :box,
        grid = false,
        minorticks = true,
        guidefontvalign = :top,
        guidefonthalign = :right,
        foreground_color_legend = nothing,
        lab = "",
    )
end

# ╔═╡ 5da141a7-c68c-470f-9b5e-8083a8367b9d
md"""
# Projections of Dalitz
"""

# ╔═╡ f7967859-2aa2-44b0-8d21-ff4f7b53ab1c
md"""
We use the expression of the phase space in term of $\cos\theta_{ij}$ and $m_{ij}^2$,

$d \Phi_3 = \frac{\lambda^{1/2}(\sigma_k, m_i^2, m_j^2) \lambda^{1/2}(m_0^2, m_k^2, \sigma_k)}{\sigma_k} \frac{d \cos\theta_{ij}}{2} d \sigma_k$


The differential cross section reads,

$
\frac{\Gamma}{\cos\theta_{ij}} = \int \left| \mathcal{M} \right|^2 \frac{\lambda^{1/2}(\sigma_k, m_i^2, m_j^2) \lambda^{1/2}(m_0^2, m_k^2, \sigma_k)}{\sigma_k}\,d \sigma_k$
"""

# ╔═╡ 15da60bc-9dc6-48c4-95e0-15079822104f
function project_cosθij_integrant(fs, ms, z; k)
    #
    i, j = ij_from_k(k)
    mi, mj, mk, m0 = ms[i], ms[j], ms[k], ms[4]
    #
    function integrand(xσk)
        σs = x2σs([xσk, (z + 1) / 2], ms; k)
        σk = σs[k]
        jac = sqrt(Kallen(σk, mk^2, m0^2) * Kallen(σk, mi^2, mj^2)) / σk
        return fs(σs) / 2 * jac * ((m0 - mk)^2 - (mi + mj)^2)
    end
    return integrand
end

# ╔═╡ 7a49d79e-73ba-4b61-8f13-00fa8b15f765
md"""
## What do I integrate
"""

# ╔═╡ 7dd4f654-1a7d-4a41-8b14-72422ebdf51c
md"""
## Load the model
"""

# ╔═╡ 857440e6-f0e4-4e67-919d-60717acf512b
begin
    input = open(joinpath(@__DIR__, "x2pipipi-compass-1391643.json")) do io
        JSON.parse(io)
    end

    workspace = Dict{String,Any}()

    @unpack functions = input
    for fn in functions
        @unpack name, type = fn
        instance_type = eval(Symbol(type))
        workspace[name] = dict2instance(instance_type, fn)
    end

    @unpack distributions = input
    for dist in distributions
        @unpack name, type = dist
        instance_type = eval(Symbol(type))
        workspace[name] = dict2instance(instance_type, dist; workspace)
    end
end

# ╔═╡ 71f12722-c866-4954-9f1c-ce93f1f28e83
model_names = filter(x -> occursin("compass", x), keys(workspace)) |> collect

# ╔═╡ e845a971-315c-4f03-9e3f-d9e6cdfc1bb0
model_name = model_names[4]

# ╔═╡ 399664e3-fbb6-40a4-ac7d-479a98ec3e48
const model = workspace[model_name].model;

# ╔═╡ f840f8bc-4606-4240-911d-bb1894c8f6d0
md"""
## Kinematics
"""

# ╔═╡ 06b54e37-8859-4d6c-a509-ebb690628ae5
ms = model.chains[1].tbs.ms;

# ╔═╡ 6f1f3684-e7a4-49ae-b144-0235f3855af5
plot(
    layout = grid(1, 3),
    size = (800, 250),
    plot(σs -> σs[1], ms, c = cgrad(:roma, 10, categorical = true), aspect_ratio = 1),
    plot(σs -> σs[2], ms, c = cgrad(:roma, 10, categorical = true), aspect_ratio = 1),
    plot(σs -> σs[3], ms, c = cgrad(:roma, 10, categorical = true), aspect_ratio = 1),
)

# ╔═╡ b60e2811-75f8-43e4-8c81-8b119aa2bddc
plot(
    layout = grid(1, 3),
    size = (800, 250),
    plot(
        σs -> cosθ23(σs, ms^2),
        ms,
        c = cgrad(:roma, 10, categorical = true),
        aspect_ratio = 1,
    ),
    plot(
        σs -> cosθ12(σs, ms^2),
        ms,
        c = cgrad(:roma, 10, categorical = true),
        aspect_ratio = 1,
    ),
    plot(
        σs -> cosθ31(σs, ms^2),
        ms,
        c = cgrad(:roma, 10, categorical = true),
        aspect_ratio = 1,
    ),
)

# ╔═╡ 4bf7f8a2-bf99-4a32-b101-05393a01c6bd
I(σs) = unpolarized_intensity(model, σs)

# ╔═╡ 7f3ba720-f43b-412c-87d0-a0f0fb5d162c
begin
    plot(ms, I; iσx = 3, iσy = 1, aspect_ratio = 1)
    plot!(border12(ms), l = (3, :black), aspect_ratio = 1)
end

# ╔═╡ 25f5a348-6344-4527-b704-1363bb6d305f
begin
    p = plot(
        layout = grid(2, 3),
        size = (800, 500),
        yaxis = nothing,
        map(1:3) do k
            i, j = ij_from_k(k)
            plot(0:0.03:2, xlab = "m($i$j)") do σk
                integrand = projection_integrand(I, ms, σk; k)
                quadgk(integrand, 0, 1)[1]
            end
        end...,
        map(1:3) do k
            i, j = ij_from_k(k)
            plot(shift_by_half(-1:0.03:1), xlab = "cosθ($i$j)") do z
                integrand = project_cosθij_integrant(I, ms, z; k)
                quadgk(integrand, 0, 1)[1]
            end
        end...,
    )
    plot!(
        layout = grid(2, 1, heights = (0.01, 0.99)),
        size = (800, 500 / 0.9),
        plot(title = model_name, xaxis = false, yaxis = false),
        p,
    )
    savefig("~/Downloads/dalitz_projections_$(model_name).pdf")
    plot!()
end

# ╔═╡ 0252d2fa-cf26-4e8e-9176-805b392c3b52
begin
    plot(layout = (@layout [A B{0.8w,0.8h}; D C]), size = (700, 700), link = :both)
    plot!(sp = 1, 0:0.01:2, permute = (:x, :y), yaxis = nothing, xlab = "m²(12)") do σk
        integrand = projection_integrand(I, ms, σk; k = 1)
        quadgk(integrand, 0, 1)[1]
    end
    plot!(sp = 2, I, ms, iσx = 1, iσy = 3, inset = (2, bbox(0.65, 0.0, 0.35, 0.35)))
    plot!(sp = 3, xaxis = false, yaxis = false, top_margin = -10mm)
    plot!(sp = 4, 0:0.01:2, yaxis = nothing, xlab = "m²(23)", top_margin = -10mm) do σk
        integrand = projection_integrand(I, ms, σk; k = 3)
        quadgk(integrand, 0, 1)[1]
    end
    plot!(sp = 2, xaxis = nothing, yaxis = nothing)
    plot!(sp = 5, 0:0.01:2, yaxis = nothing, xlab = "m²(31)") do σk
        integrand = projection_integrand(I, ms, σk; k = 2)
        quadgk(integrand, 0, 1)[1]
    end
end

# ╔═╡ Cell order:
# ╟─5da141a7-c68c-470f-9b5e-8083a8367b9d
# ╠═f787599e-45c8-11ef-2dd4-61b67a34bcbd
# ╟─f7967859-2aa2-44b0-8d21-ff4f7b53ab1c
# ╠═15da60bc-9dc6-48c4-95e0-15079822104f
# ╟─7a49d79e-73ba-4b61-8f13-00fa8b15f765
# ╟─6f1f3684-e7a4-49ae-b144-0235f3855af5
# ╟─b60e2811-75f8-43e4-8c81-8b119aa2bddc
# ╟─7dd4f654-1a7d-4a41-8b14-72422ebdf51c
# ╠═857440e6-f0e4-4e67-919d-60717acf512b
# ╠═71f12722-c866-4954-9f1c-ce93f1f28e83
# ╠═e845a971-315c-4f03-9e3f-d9e6cdfc1bb0
# ╠═399664e3-fbb6-40a4-ac7d-479a98ec3e48
# ╟─f840f8bc-4606-4240-911d-bb1894c8f6d0
# ╠═06b54e37-8859-4d6c-a509-ebb690628ae5
# ╠═4bf7f8a2-bf99-4a32-b101-05393a01c6bd
# ╠═7f3ba720-f43b-412c-87d0-a0f0fb5d162c
# ╠═25f5a348-6344-4527-b704-1363bb6d305f
# ╠═0252d2fa-cf26-4e8e-9176-805b392c3b52
