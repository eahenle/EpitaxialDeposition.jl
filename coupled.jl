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
	using StatsBase
end

# ╔═╡ ceac564a-9ff6-11ec-0ba7-13f6890e2021
md"""
# Coupled Kinetics and Transport
"""

# ╔═╡ 632a0cc6-6c8f-4a93-b72c-d6802fa3fef4
md"""
!!! note \" Linguistic Side-Note \"
	In German, compound words are often used to express novel ideas. Example: when one is simultaneously proud and embarrassed, the feeling is called "Schamstolz" (portmanteau of the words for "shame" and "pride").

	This code fills me with Schamstolz.
"""

# ╔═╡ 4411f086-c12c-45c4-a5b0-bd6d9ac2bcd9
md"""
## Set up initial conditions and static variables
"""

# ╔═╡ 25f27d9b-4acb-4803-a2c1-845801de7ba7
md"""
!!! warning \"Caveat Quaesitor\"
	The runtime is approximately 10-20 seconds per simulated second.  
	
	Yep.  It's literally faster to do the experiment!
"""

# ╔═╡ 54e026fa-12c1-4906-b96e-88b11eabfb38
begin
	y_SiCl₄ = 0.5
	Psys = 760.
	temperature = 1000.
	
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

	timesteps = 5

	Δt = 1

	Lx = 30
	Ly = 30
	Nx = Ny = 150
	x  = [1:Nx...] * Lx / Nx
	y  = [1:Ny...] * Ly / Ny

	wafer_diameter = 10

	D = Dict(
		:SiCl₄ 	=> 10 * DAB(:SiCl₄, :H₂),
		:HCl 	=> 10 * DAB(:HCl, :H₂),
		:H₂ 	=> 10 * DAB(:H₂, :H₂),
		:SiCl₂ 	=> 10 * DAB(:SiCl₂, :H₂)
	)

	# initial concentration of each species
	local udict = Dict(u₀)
	c₀ = Dict([
		:H₂ 	=> udict[:H₂],
		:SiCl₄ 	=> udict[:SiCl₄],
		:SiCl₂ 	=> udict[:SiCl₂],
		:HCl 	=> udict[:HCl]
	])

	mesh = createMesh2D(Nx + 2, Ny + 2, Lx, Ly)

	wafer_bounds = EpitaxialDeposition.wafer_boundary(Nx, Lx, wafer_diameter)

	# book-keeping variables
	colnames = Symbol.(String.([u[1] for u in u₀]) .* "(t)")
	sol_to_u = Dict([name => Symbol(chop(String(name), tail=3)) for name in colnames])
	u_to_sol = Dict([value => key for (key, value) in sol_to_u])
end;

# ╔═╡ c1d044da-26df-4eb4-96fa-a16349a2d396
md"""
## Run a coupled simulation with Catalyst/JFVM
"""

# ╔═╡ 2492cf1b-f352-46e1-8e2a-1e2877568b41
try
	LocalResource("dfd.png")
catch
	try
		LocalResource("../dfd.png")
	catch e
		@warn e
	end
end

# ╔═╡ f00a16d4-9fbd-4ef8-aa0f-93476443c372
function transport_step(m::MeshStructure, D::Float64, c₀::Matrix{Float64}, Nx::Int, 
		ϕ::Float64; N_steps=5, Lx=30., Ly=30.)
    # Define boundary conditions
    BC = createBC(m)

    BC.left.a[:] .= 1.0 
    BC.right.a[:] .= 1.0

    BC.left.b[:] .= 0.0
    BC.right.b[:] .= 0.0

    BC.left.c[:] .= 0.0 
    BC.right.c[:] .= 0.0

    # top and bottom boundary conditions for each coefficient
    BC.top.a[1,:] .= 1.0 
    BC.bottom.a[1,:] .= 1.0

    BC.top.b[1,:] .= 0.0 
    BC.bottom.b[1,:] .= 0.0

    BC.top.c[1,:] .= 0.0 
    BC.bottom.c[1,:] .= 0.0

    BC.bottom.c[wafer_bounds, 1] .= -ϕ

    # Give a value for the diffusion coefficient based on current system
    D_cell = createCellVariable(m, D) # assign the diffusion coefficient as a variable to each cell of the mesh
    D_face = geometricMean(D_cell) # choose an averaging scheme, for how the diffusion coefficient on a cell face is calculated
    
    c = createCellVariable(m, c₀, BC) # assign the "initial" species concentration to cells [mol/L]
    
    # Discretize the problem and build our solution
    M_diff = diffusionTerm(D_face) # matrix of diffusion term coefficients
    (M_bc, RHS_bc) = boundaryConditionTerm(BC) # matrix composed of coefficients and the right hand side for the BC

	Δt = EpitaxialDeposition.PARAMS[:Δt][:JFVM]
    
    for i = 1:N_steps
        (M_t, RHS_t) = transientTerm(c, Δt, 1.0)
        M = M_t - M_diff + M_bc # adding together all sparse matrices of coefficients
        RHS = RHS_bc + RHS_t # add all RHS's to each other
        c = solveLinearPDE(m, M, RHS)
    end

	return c
end

# ╔═╡ 059e8e10-d1d1-4499-b884-8ff8cec42bad
function coupled_simulation()
	# arrays to store time points
	kinetics = []
	transport = []
	# variables to store states
	u = u₀ # surface concentration
	c = Dict([key => ones(Nx + 2, Ny + 2) * value for (key, value) in c₀]) # cells
	ϕ = Dict([key => 1e-32 for (key, value) in c₀]) # surface Δc (pseudo-zero)

	# loop over time points
	for t in 1:timesteps
		# run the kinetic model
		k_sol = solve(ODEProblem(deposition_rxn_network, u, ((t - 1) * Δt, t * Δt)), Tsit5(), saveat=EpitaxialDeposition.PARAMS[:Δt][:Catalyst], maxiters=1e8)

		# run the transport models
		t_sol = Dict([key => transport_step(mesh, value, c[key], Nx, ϕ[key], N_steps=20) for (key, value) in D])

		# update the transport models' cell states for next iteration
		c = Dict([species => cell.value[2:end-1, 2:end-1] for (species, cell) in t_sol])

		u_sol = Dict([sol_to_u[col] => DataFrame(k_sol)[end, col] for col in colnames])

		# change in concentrations due to chemical reactions
		old_u = Dict(u)
		ϕ = Dict([species => new_u - old_u[species] for (species, new_u) in u_sol])
		
		for (species, ci) in c
			c[species][wafer_bounds, 1] .+= ϕ[species] #- Δu[species]
		end
		
		# update the kinetic model's state for next iteration
		u′ = [
	        :SiCl₄   => mean(c[:SiCl₄][wafer_bounds, 1]),
	        :H₂      => mean(c[:H₂][wafer_bounds, 1]),
	        :Si_dep  => u_sol[:Si_dep],
	        :Si_etch => u_sol[:Si_etch],
	        :HCl     => mean(c[:HCl][wafer_bounds, 1]), 
	        :SiCl₂   => mean(c[:SiCl₂][wafer_bounds, 1]),
	        :k       => u_sol[:k],
	        :T       => u_sol[:T]
	    ]
		u = u′

		# append the output data
		push!(kinetics, k_sol)
		push!(transport, t_sol)
	end

	# merge timestep outputs
	kinetics = reduce(append!, DataFrame.(kinetics))

	# calculate film deposition
	δ = film_thickness.(kinetics[:, "Si_dep(t)"], kinetics[:, "Si_etch(t)"])

	return kinetics, transport, δ
end

# ╔═╡ 0c83de9c-9f4e-4f3c-97e1-efe01fecc4a3
kinetics, transport, δ = coupled_simulation();

# ╔═╡ 48063896-50eb-4e95-96a8-c563e8baaf8c
md"""
## Results
"""

# ╔═╡ 5296cf87-ee7f-483b-ba52-996ea974aaec
begin
    local fig = Figure()
    local ax = Axis(
        fig[1,1],
        title = "Surface Gas-Film Species",
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
Profile snapshot: $(@bind t_i PlutoUI.NumberField(0:length(transport)-1, default=0))
"""

# ╔═╡ 600ebffc-d978-4813-9d06-67208b11a61d
begin
	local i = round(Int, floor(t_i)) + 1
	plot_species_profiles(transport[i], x, y)
end

# ╔═╡ Cell order:
# ╟─ceac564a-9ff6-11ec-0ba7-13f6890e2021
# ╟─632a0cc6-6c8f-4a93-b72c-d6802fa3fef4
# ╠═011c70e0-83e8-4ec7-9dc7-19c9dd69142e
# ╟─4411f086-c12c-45c4-a5b0-bd6d9ac2bcd9
# ╟─25f27d9b-4acb-4803-a2c1-845801de7ba7
# ╠═54e026fa-12c1-4906-b96e-88b11eabfb38
# ╟─c1d044da-26df-4eb4-96fa-a16349a2d396
# ╟─2492cf1b-f352-46e1-8e2a-1e2877568b41
# ╠═f00a16d4-9fbd-4ef8-aa0f-93476443c372
# ╠═059e8e10-d1d1-4499-b884-8ff8cec42bad
# ╠═0c83de9c-9f4e-4f3c-97e1-efe01fecc4a3
# ╟─48063896-50eb-4e95-96a8-c563e8baaf8c
# ╟─5296cf87-ee7f-483b-ba52-996ea974aaec
# ╟─600ebffc-d978-4813-9d06-67208b11a61d
# ╠═7da4da48-be95-4df4-ba69-8927d02da60b
