### A Pluto.jl notebook ###
# v0.18.1

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

# ╔═╡ 011c70e0-83e8-4ec7-9dc7-19c9dd69142e
begin
	import Pkg
	Pkg.activate(".")
	using EpitaxialDeposition
	using DataFrames
	using CairoMakie
	using JFVM
	using PlutoUI
end

# ╔═╡ ceac564a-9ff6-11ec-0ba7-13f6890e2021
md"""
# Coupled Kinetics and Transport
"""

# ╔═╡ 4411f086-c12c-45c4-a5b0-bd6d9ac2bcd9
# set up initial conditions and static variables

# ╔═╡ 54e026fa-12c1-4906-b96e-88b11eabfb38
begin
	y_SiCl₄ = 0.5
	Psys = 1.
	temperature = 1100.
	
	u₀ = [
	        :SiCl₄   => (y_SiCl₄ * Psys) / EpitaxialDeposition.PARAMS[:R][:L_Torr_mol_K] / temperature,
	        :H₂      => Psys * (1 - y_SiCl₄) / EpitaxialDeposition.PARAMS[:R][:L_Torr_mol_K] / temperature,
	        :Si_dep  => 0.0,
	        :Si_etch => 0.0,
	        :HCl     => 0.0, 
	        :SiCl₂   => 0.0,
	        :k       => EpitaxialDeposition.PARAMS[:kKp][:k],
	        :T       => temperature
	    ]
	
	t_max = 10.  

	timesteps = 2

	Δt = 1

	Lx = 30
	Ly = 30
	Nx = Ny = 150
	x = [1:Nx...] * Lx / Nx
	y = [1:Ny...] * Ly / Ny

	wafer_diameter = 10

	D = Dict(
		:SiCl₄ => DAB(:SiCl₄, :H₂),
		:HCl => DAB(:HCl, :H₂),
		:H₂ => DAB(:H₂, :H₂),
		:SiCl₂ => DAB(:SiCl₂, :H₂)
	)

	# initial concentration of each species
	c₀ = Dict([
		:H₂ => 1.,
		:SiCl₄ => 1.,
		:SiCl₂ => 0.,
		:HCl => 0.
	])

	mesh = createMesh2D(Nx + 2, Ny + 2, Lx, Ly)
end

# ╔═╡ c1d044da-26df-4eb4-96fa-a16349a2d396
# run a Catalyst sim with call-outs to JFVM

# ╔═╡ 2c89c6df-1e94-43b2-9624-280463aaf2b9
function f()
	kinetics = []
	transport = []
	u = u₀

	colnames = Symbol.(String.([u[1] for u in u₀]) .* "(t)")
	sol_to_u = Dict([name => Symbol(chop(String(name), tail=3)) for name in colnames])

	# loop over time points
	for t in 1:timesteps
		# run the kinetic model
		k_sol = solve(ODEProblem(deposition_rxn_network, u, ((t - 1) * Δt, t * Δt)), Tsit5(), saveat=EpitaxialDeposition.PARAMS[:Δt][:Catalyst], maxiters=1e8)

		# update the kinetic model's state for next iteration
		u = [sol_to_u[col] => DataFrame(k_sol)[end, col] for col in colnames]

		# append the kinetics data
		push!(kinetics, k_sol)

		# run the transport model
		t_sol = Dict([key => trans_diff_Neumann(mesh, value, c₀[key], Nx) for (key, value) in D])

		# append the transport data
		push!(transport, t_sol)
	end

	# merge timestep outputs
	kinetics = reduce(append!, DataFrame.(kinetics))

	# calculate film deposition
	δ = film_thickness.(kinetics[:, "Si_dep(t)"], kinetics[:, "Si_etch(t)"])

	return kinetics, transport, δ
end

# ╔═╡ 0c83de9c-9f4e-4f3c-97e1-efe01fecc4a3
kinetics, transport, δ = f();

# ╔═╡ 48063896-50eb-4e95-96a8-c563e8baaf8c
# visualize the results

# ╔═╡ 5296cf87-ee7f-483b-ba52-996ea974aaec
begin
    local fig = Figure()
    local ax = Axis(
        fig[1,1],
        title = "Gas-Phase Species",
        ylabel="Concentration [mmol/L]",
        xlabel="Time [s]"
    )

    plotted_species = setdiff(
        String.(Symbol.(species(deposition_rxn_network))),
        ["T(t)", "k(t)", "Si_dep(t)", "Si_etch(t)"]
    )
    
    for name in plotted_species
        lines!(kinetics[:, :timestamp], 1000 .* kinetics[:, name], label=chop(name, tail=3)) #mol/L --> mmol/L
    end

    fig[:,2] = Legend(fig, ax)

    Axis(
        fig[2,1],
        title="Si Film Growth on a $wafer_diameter cm Wafer, T = $temperature [K]",
        xlabel="Time [s]",
        ylabel="δ [μm]"
    )
    
    lines!(kinetics[:, :timestamp], δ)

    fig
end

# ╔═╡ 7da4da48-be95-4df4-ba69-8927d02da60b
## THIS IS PICKING TIMESTEP, NOT ACTUAL TIME

md"""
Profile snapshot time: $(@bind t_i PlutoUI.Slider(0:Δt:t_max, default=0, show_value=true)) s
"""

# ╔═╡ 600ebffc-d978-4813-9d06-67208b11a61d
begin
	local i = round(Int, floor(t_i)) + 1
	plot_species_profiles(transport[i], x, y)
end

# ╔═╡ Cell order:
# ╟─ceac564a-9ff6-11ec-0ba7-13f6890e2021
# ╠═011c70e0-83e8-4ec7-9dc7-19c9dd69142e
# ╠═4411f086-c12c-45c4-a5b0-bd6d9ac2bcd9
# ╠═54e026fa-12c1-4906-b96e-88b11eabfb38
# ╠═c1d044da-26df-4eb4-96fa-a16349a2d396
# ╠═2c89c6df-1e94-43b2-9624-280463aaf2b9
# ╠═0c83de9c-9f4e-4f3c-97e1-efe01fecc4a3
# ╠═48063896-50eb-4e95-96a8-c563e8baaf8c
# ╠═5296cf87-ee7f-483b-ba52-996ea974aaec
# ╠═600ebffc-d978-4813-9d06-67208b11a61d
# ╠═7da4da48-be95-4df4-ba69-8927d02da60b
