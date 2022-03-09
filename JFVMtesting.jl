### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ d6b3a9fa-3c6d-422d-afff-30e49502e145
begin
	import Pkg; using Pkg
	Pkg.add(url="https://github.com/simulkade/JFVMvis.jl")
	Pkg.add(url="https://github.com/simulkade/JFVM.jl")
	# Pkg.add(["Dierckx", "CairoMakie", "PyPlot", "LsqFit"])
	
end

# ╔═╡ d2eb8940-9b16-11ec-1d74-3bc00f5f950f
using JFVM, CairoMakie, LsqFit

# ╔═╡ 6e4302dc-b2ce-4d8f-809b-078f2c39c2dc
md"""
### Diffusion Coefficients ``\mathscr{D}``
Here we'll start with gas phase binary diffusion coefficients for each species (in hydrogen), and later consider mixtures.

The first correlation proposed by Hirschfelder, Bird, and Spotz is for dilute, **nonpolar**, spherical, **nonreactive** gases (Fundamentals of Momentum, Heat, and Mass Transfer 6th Ed, Welty et al.). The second correlation is a modification by Fuller, Schettler, and Giddings, which allow for evaluating diffusivity when Lennard-Jones parameteres aren't available:

1. `` D_{AB}=\frac{0.001858T^{3/2}(\frac{1}{M_A}+\frac{1}{M_B})^{1/2}}{P \sigma^2_{AB} \Omega_D} `` with:

- D = diffusion of A through B [cm²/s]
- P = absolute pressure [atm]
- MW = molecular weight [g/gmol]
- T = absolute temperature [K]
- σ = collision diameter of the two species [Å]
- Ω = collision integral for molecular diffusion

Values for σ and a parameter ϵ are tabulated in the book referenced, Appendix K

2. ``D_{AB} = \frac{0.001T^{1.75}(\frac{1}{M_A}+\frac{1}{M_B})^{1/2}}{P [(\Sigma \nu_i)_A ^{1/3}+(\Sigma \nu_i)_B ^{1/3}]^2}  ``
- ``\nu_i`` is an incremental quantity dependent on a molecule's functional groups, with values tabulated in Table 24.3 of the textbook

We will use method 1, using critical parameters (Tc, Pc) to estimate the collision diameter and energy ϵ for species which we can gather critical parameter information.

"""

# ╔═╡ 84ecce92-8d75-4a67-8f71-ca63a6da1137
#= TO-DO:

2) Add something more to the c₀ for each species (most likely values exported from Catalyst)

=#

# ╔═╡ e137888b-fc7e-4ddb-b25e-3c53bb8ee23a
begin
	# random system parameters and necessary constants
	const κ = 1.38e-16 # [ergs/K]
	const T = 1200 # [K]
	const P = 1 # [atm]
end

# ╔═╡ 457a34ad-8514-4c6a-b16b-c94e6515ba40
begin

	# note that we're pulling critical parameters for SiCl₂ by taking them to be equal to SiH₂Cl₂
	MW_dict = Dict(:SiCl₄ => 169.9, :H₂ => 2.0, :HCl => 36.5, :SiCl₂ => 99) # [g/gmol]
	σ_dict = Dict(:SiCl₄ => 5.08, :H₂ => 2.968, :HCl => 3.305) # [Å]
	ϵ_dict = Dict(:SiCl₄ => 358*κ, :H₂ => 33.3*κ, :HCl => 360*κ) # [ergs]
	ν_dict = Dict(:HCl => [19.5, 1.98], :H₂ => 7.07)
	Tc_dict = Dict(:HCl => 324.68, :SiCl₄ => 508.1, :H₂ => 33.18, :SiCl₂ => 178.9 + 273) # [K]
	Pc_dict = Dict(:HCl => 82.56 / 1.013, :SiCl₄ => 35.93 / 1.013, :H₂ => 13.00 / 1.013, :SiCl₂ => 44) # [atm]

	Params = Dict(:MW => MW_dict, 
				  :σ => σ_dict, 
				  :ϵ => ϵ_dict, 
				  :ν => ν_dict,
				  :Pc => Pc_dict,
				  :Tc => Tc_dict) # store all data in a single dictionary
end

# ╔═╡ 1333c9f1-6202-4df9-9b09-4d19f2b91d65
begin
	# load in Lennard-Jones Constants (App. K, Welty et al.)
	a = [0.30 2.785 2.662 1.80 1.221 1.116
	0.35 2.628 2.476 1.85 1.209 1.105
	0.40 2.492 2.318 1.90 1.197 1.094
	0.45 2.368 2.184 1.95 1.186 1.084
	0.50 2.257 2.066 2.00 1.175 1.075
	0.55 2.156 1.966 2.10 1.156 1.057
	0.60 2.065 1.877 2.20 1.138 1.041
	0.65 1.982 1.798 2.30 1.122 1.026
	0.70 1.908 1.729 2.40 1.107 1.012
	0.75 1.841 1.667 2.50 1.093 0.9996
	0.80 1.780 1.612 2.60 1.081 0.9878
	0.85 1.725 1.562 2.70 1.069 0.9770
	0.90 1.675 1.517 2.80 1.058 0.9672
	0.95 1.629 1.476 2.90 1.048 0.9576
	1.00 1.587 1.439 3.00 1.039 0.9490
	1.05 1.549 1.406 3.10 1.030 0.9406
	1.10 1.514 1.375 3.20 1.022 0.9328
	1.15 1.482 1.346 3.30 1.014 0.9256
	1.20 1.452 1.320 3.40 1.007 0.9186
	1.25 1.424 1.296 3.50 0.9999 0.9120
	1.30 1.399 1.273 3.60 0.9932 0.9058
	1.35 1.375 1.253 3.70 0.9870 0.8998
	1.40 1.353 1.233 3.80 0.9811 0.8942
	1.45 1.333 1.215 3.90 0.9755 0.8888
	1.50 1.314 1.198 4.00 0.9700 0.8836
	1.55 1.296 1.182 4.10 0.9649 0.8788
	1.60 1.279 1.167 4.20 0.9600 0.8740
	1.65 1.264 1.153 4.30 0.9553 0.8694
	1.70 1.248 1.140 4.40 0.9507 0.8652
	4.50 0.9464 0.8610 10.0 0.8242 0.7424
	4.60 0.9422 0.8568 20.0 0.7432 0.6640
	4.70 0.9382 0.8530 30.0 0.7005 0.6232
	4.80 0.9343 0.8492 40.0 0.6718 0.5960
	4.90 0.9305 0.8456 50.0 0.6504 0.5756
	5.0 0.9269 0.8422 60.0 0.6335 0.5596
	6.0 0.8963 0.8124 70.0 0.6194 0.5464
	7.0 0.8727 0.7896 80.0 0.6076 0.5352
	8.0 0.8538 0.7712 90.0 0.5973 0.5256];

	
	b = vcat(a[:, [1, 3]], a[:, [4,6]]);

	# order the x/y data via some mask which orders the array in increasing order
	sp = sortperm(b[:, 1])
	
	xdata = b[:, 1][sp][60:end, 1] # values of κT/ϵ
	
	ydata = b[:, 2][sp][60:end, 1] # values of Ω
	
end

# ╔═╡ fb36836c-bc8c-4a5d-812b-5e8eceaed845
begin
	
	@. model(x,p) = p[1] / (x*p[2] + p[3]) + p[4] # inverse x model
	
	p0 = [0.25, 0.75, 0.5, 0.5]
	
	fit = curve_fit(model, xdata, ydata, p0);

	collision_interp(x) = model(x, fit.param)

end

# ╔═╡ 6e4c3904-6543-41e9-a862-22cb749ed178
begin

	xs = range(xdata[1], xdata[end], length = 150)
	
	fig = Figure()
	ax = Axis(fig[1,1],
		title = "Least Squares Fit for Linear Interpolation",
		ylabel = "Collision Integral Ω",
		xlabel = "κT/ϵ")

	scatter!(xdata,ydata)
	lines!(xs, model(xs, fit.param), color = :blue)

	fig

end

# ╔═╡ c26da7bb-b415-4d74-9dd2-89c675724874
function DAB(species1::Symbol, species2::Symbol, T, P)

	# calculate the collision diameter from critical parameters
	σ1 = 2.44*(Params[:Tc][species1] / Params[:Pc][species1])^(1/3)
	σ2 = 2.44*(Params[:Tc][species2] / Params[:Pc][species2])^(1/3)

	# calculate the ϵ's
	ϵ1 = 0.77 * κ * Params[:Tc][species1]
	ϵ2 = 0.77 * κ * Params[:Tc][species2]

	# combined parameters for calculating the diffusion coefficient
	σ12 = (σ1 + σ2) / 2
	ϵ12 = sqrt(ϵ1 * ϵ2)

	# calculate the value of κT/ϵ that we need to interpolate between
	Ω = collision_interp(κ*T/ϵ12)
	
	
	return 0.001858 * T^(3/2) * sqrt(1 / Params[:MW][species1] + 1 / Params[:MW][species2]) / P / σ12^2 / Ω # [cm²/s]
	
end

# ╔═╡ 00f98dbf-4df2-4972-8f3e-15425b735077
D_SiCl₄= DAB(:SiCl₄, :H₂, T, P)

# ╔═╡ 9912d690-5dd6-45d9-94cf-8e0ca301fd29
D_HCl = DAB(:HCl, :H₂, T, P)

# ╔═╡ ac05f8e1-e284-48d0-8c81-c9b7677d8d52
D_H₂ = DAB(:H₂, :H₂, T, P)

# ╔═╡ 86c58901-b84d-45f7-a59e-521337b9ec36
D_SiCl₂ = DAB(:SiCl₂, :H₂, T, P) # (dichlorosilane DAB)

# ╔═╡ 680313b8-d39c-428b-adb5-a10dee051860
md"""
## Finite Volume Method for Showing Concentration Gradients
"""

# ╔═╡ 7c259d4c-122c-4270-ad64-185b09be4a2f
begin
	
	# Start by generating a uniform mesh for our system, defined on the size of the domain and number of cells in each direction

	# Non-uniform meshes can be generated by defining an array of "cell positions" and passing it into `createMesh1D` instead of Nx and Lx
	
	Nx = 100 # number of cells in x direction
	Ny = 60 # numbers of cells in y direction
	Lx = 30.0 # [cm] - domain length in x direction
	Ly = 3.0 # [cm] - domain length in y direction

	N_steps = 100
	
	m = createMesh2D(Nx, Ny, Lx, Ly); # generate a uniform mesh (2D)

	# encode each cell location to its corresponding location in true space
	x_cell_reshape = [1:Nx...] * (Lx/Nx); 
	y_cell_reshape = [1:Ny...] * (Ly/Ny);
	

end

# ╔═╡ d4171c29-94b5-4135-aaac-3c9bcca63cc4
wafer_diameter = 20 # [cm]

# ╔═╡ ea81c3e9-abe2-4158-b25d-f5f6812778bb
function wafer_boundary(Nx, Lx, wafer_diameter)  
	@assert wafer_diameter < Lx "Wafer too big"
	
	half = ceil(Nx / 2 + 1)
	dticks = wafer_diameter * Nx / Lx 
	left, right = half - dticks/2, half + dticks/2

	@assert left ≥ 0.5 && right < Nx + 2
	
	return round(Int, left):round(Int, right)
end

# ╔═╡ 6a944590-be5a-498f-aaf0-8966f610dc61
function trans_diff_Neumann(m::MeshStructure, D::Float64, c₀::Float64; N=N_steps)

	# Define our boundary conditions. Package/function defaults yield Neumann boundary conditions, though Dirichlet and Robin conditions are also supported, and each cell face can have a unique boundary condition
	BC = createBC(m)

	# left and right boundary conditions for each coefficient in the boundary condition expression 
	
	#link: https://nbviewer.org/github/)simulkade/JFVM.jl/blob/master/examples/jfvm_tutorial.ipynb
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

	BC.bottom.c[wafer_boundary(Nx, Lx, wafer_diameter), 1] .= -1.0
		

	# Give a value for the diffusion coefficient based on current system
	D = D # [cm²/s]
	D_cell = createCellVariable(m, D) # assign the diffusion coefficient as a variable to each cell of the mesh
	D_face = geometricMean(D_cell) # choose an averaging scheme, for how the diffusion coefficient on a cell face is calculated
	
	c = createCellVariable(m, c₀, BC) # assign the "initial" species concentration to cells [mol/L]
	
	# Discretize the problem and build our solution
	M_diff = diffusionTerm(D_face) # matrix of diffusion term coefficients
	(M_bc, RHS_bc) = boundaryConditionTerm(BC) # matrix composed of coefficients and the right hand side for the BC



	### modify to be 1/5th of the timestep from Catalyst
	
	if Lx < Ly # calculate a time step based on the smaller length (for more resolution)
	Δt = sqrt(Lx^2 / D) / N # recommended time step calculation

	else
	Δt = sqrt(Ly^2 / D) / N

	end

	# this code seems to be for taking discrete time steps forward, with size of Δt

	### parameterize the # of steps?
	
	for i = 1:10
		(M_t, RHS_t) = transientTerm(c, Δt, 1.0)
		M = M_t - M_diff + M_bc # adding together all sparse matrices of coefficients
		RHS = RHS_bc + RHS_t # add all RHS's to each other
		c = solveLinearPDE(m, M, RHS)
	end

return c

end

# ╔═╡ c3e3ebfd-0722-4cc0-b647-15e10ef93a82
function plot_species_profile!(species::String, D::Float64, m::MeshStructure, c₀)

	sol = trans_diff_Neumann(m, D, c₀) #set up and completely solve discretization

	fig = CairoMakie.Figure()
	ax = Axis(fig[1,1],
				xlabel = "x-direction",
				ylabel = "y-direction",
				title = "Concentration Profile of $species")
	
	heatmap!(x_cell_reshape, y_cell_reshape, sol.value)

	Colorbar(fig[1,2], limits = (0, maximum(sol.value)), colormap = :viridis,
			label = "$species Concentration [mol/L]")

	fig

	return sol, fig
end

# ╔═╡ 20896c36-caaf-4dea-b5f2-0b71af22af79
sol1, fig1 = plot_species_profile!("SiCl₄", D_SiCl₄, m, 0.0);

# ╔═╡ 59f80b8d-2e48-4b71-bef0-1327888158b6
fig1

# ╔═╡ 21a4c01a-a92c-456f-86f7-663d576845e2
sol2, fig2 = plot_species_profile!("HCl", D_HCl, m, 0.0);

# ╔═╡ b34d9e43-357f-45d0-8406-42a55b116bee
fig2

# ╔═╡ f73c1d88-1605-4445-9c9c-ddc16937cc48
function Dab_nonpolar(species1::Symbol, species2::Symbol, T, P)

	# Binary diffusion coefficient of species A (1) through species B (2)
	
	σ_AB = (Params[:σ][species1] + Params[:σ][species2]) / 2
	
	ϵ_AB = sqrt(Params[:ϵ][species1] * Params[:ϵ][species2])


	Ω = 0.73 # estimated from Appendix K, didn't want to write an interpolation calculator at the time


	return 0.001858 * T^(3/2) * sqrt(1 / Params[:MW][species1] + 1 / Params[:MW][species2]) / P / σ_AB^2 / Ω # [cm²/s]
	
end

# ╔═╡ d0374368-0283-40f1-a55c-133eb1cff7e8
function Dab_polar(species1::Symbol, species2::Symbol, T, P)

	ν_species1 = sum(Params[:ν][species1])
	ν_species2 = sum(Params[:ν][species2])

	return 0.001 * T^(1.75) * sqrt(1 / Params[:MW][species1] + 1 / Params[:MW][species2]) / P / (ν_species1^(1/3) + ν_species2^(1/3))^2 # [cm²/s]
end

# ╔═╡ Cell order:
# ╠═d6b3a9fa-3c6d-422d-afff-30e49502e145
# ╠═d2eb8940-9b16-11ec-1d74-3bc00f5f950f
# ╟─6e4302dc-b2ce-4d8f-809b-078f2c39c2dc
# ╠═84ecce92-8d75-4a67-8f71-ca63a6da1137
# ╠═e137888b-fc7e-4ddb-b25e-3c53bb8ee23a
# ╠═457a34ad-8514-4c6a-b16b-c94e6515ba40
# ╠═1333c9f1-6202-4df9-9b09-4d19f2b91d65
# ╟─6e4c3904-6543-41e9-a862-22cb749ed178
# ╟─fb36836c-bc8c-4a5d-812b-5e8eceaed845
# ╟─c26da7bb-b415-4d74-9dd2-89c675724874
# ╠═00f98dbf-4df2-4972-8f3e-15425b735077
# ╠═9912d690-5dd6-45d9-94cf-8e0ca301fd29
# ╠═ac05f8e1-e284-48d0-8c81-c9b7677d8d52
# ╠═86c58901-b84d-45f7-a59e-521337b9ec36
# ╟─680313b8-d39c-428b-adb5-a10dee051860
# ╠═7c259d4c-122c-4270-ad64-185b09be4a2f
# ╠═20896c36-caaf-4dea-b5f2-0b71af22af79
# ╠═59f80b8d-2e48-4b71-bef0-1327888158b6
# ╠═21a4c01a-a92c-456f-86f7-663d576845e2
# ╠═b34d9e43-357f-45d0-8406-42a55b116bee
# ╠═d4171c29-94b5-4135-aaac-3c9bcca63cc4
# ╟─ea81c3e9-abe2-4158-b25d-f5f6812778bb
# ╟─c3e3ebfd-0722-4cc0-b647-15e10ef93a82
# ╠═6a944590-be5a-498f-aaf0-8966f610dc61
# ╟─f73c1d88-1605-4445-9c9c-ddc16937cc48
# ╟─d0374368-0283-40f1-a55c-133eb1cff7e8
