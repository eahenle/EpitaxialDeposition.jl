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

# ╔═╡ d2eb8940-9b16-11ec-1d74-3bc00f5f950f
begin
    import Pkg
    Pkg.activate(".")
    using EpitaxialDeposition
    using CairoMakie
    using StatsBase
    using PlutoUI
    using JFVM
end

# ╔═╡ 6584c779-4b89-40df-98c9-0dbf00ff6081
md"""
# Mass Transport
"""

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

# ╔═╡ 457a34ad-8514-4c6a-b16b-c94e6515ba40
md"""
!!! note 
	We use critical parameters for SiH₂Cl₂ in place of those for SiCl₂ (which are unknown)
"""

# ╔═╡ eb71d9ae-33d7-45dd-80ea-950a5a18ea78
md"""
To allow for calculating the collision integral at any temperature, we interpolate the data by fitting a model of the form:

$Ω = \frac{a}{b\kappa T/\epsilon + c} + d$
"""

# ╔═╡ 6e4c3904-6543-41e9-a862-22cb749ed178
begin

    xs = range(
        EpitaxialDeposition.WELTY_X[1], EpitaxialDeposition.WELTY_X[end], length=150
    )

    rsquared = round(sum((Ω(EpitaxialDeposition.WELTY_X) .-
        mean(EpitaxialDeposition.WELTY_Y)) .^ 2) / sum(((EpitaxialDeposition.WELTY_Y .- mean(EpitaxialDeposition.WELTY_Y)) .^ 2)), digits=3)
    
    local fig = Figure()
    ax = Axis(
        fig[1,1],
        title = "Least Squares Fit, R² = $rsquared",
        ylabel = "Collision Integral Ω",
        xlabel = "κT/ϵ"
    )

    scatter!(
        EpitaxialDeposition.WELTY_X, EpitaxialDeposition.WELTY_Y, label="Welty Data"
    )
    lines!(xs, Ω(xs), color=:red, label="Interpolation")

    fig[1, 2] = Legend(fig, ax)
    
    fig
end

# ╔═╡ 87620c55-7fbc-4004-928a-55c900a9a144
md"""
Calculated diffusion coefficients:
"""

# ╔═╡ 00f98dbf-4df2-4972-8f3e-15425b735077
D_SiCl₄ = DAB(:SiCl₄, :H₂)

# ╔═╡ 9912d690-5dd6-45d9-94cf-8e0ca301fd29
D_HCl = DAB(:HCl, :H₂)

# ╔═╡ ac05f8e1-e284-48d0-8c81-c9b7677d8d52
D_H₂ = DAB(:H₂, :H₂)

# ╔═╡ 86c58901-b84d-45f7-a59e-521337b9ec36
D_SiCl₂ = DAB(:SiCl₂, :H₂) # (dichlorosilane DAB)

# ╔═╡ 99be578d-673b-4f82-b0e4-2cd3fed8b847
D = Dict(
    :SiCl₄ => D_SiCl₄,
    :HCl => D_HCl,
    :H₂ => D_H₂,
    :SiCl₂ => D_SiCl₂
);

# ╔═╡ 680313b8-d39c-428b-adb5-a10dee051860
md"""
## Finite Volume Method for Showing Concentration Gradients
"""

# ╔═╡ 8075782d-3a72-4a43-833a-6f7a652f2cdf
md"""
## Simulation Parameters
"""

# ╔═╡ 7abb2d83-c82e-4ad0-a4f6-0c9b67633f1d
# initial concentration of each species
c₀ = Dict([
    :H₂ => 1.,
    :SiCl₄ => 1.,
    :SiCl₂ => 0.,
    :HCl => 0.
])

# ╔═╡ 77168336-8588-4f53-a236-cc517978455c
md"""
Voxel Edge $(@bind vox PlutoUI.Slider(0.:1e-3:1., default=0.2, show_value=true)) cm

Reactor Diameter $(@bind Lx PlutoUI.Slider(1.:100., default=30., show_value=true)) cm

Reactor Height $(@bind Ly PlutoUI.Slider(1.:100., default=30., show_value=true)) cm

Wafer Diameter $(@bind wafer_diam PlutoUI.Slider(1.:40., default=15., show_value=true)) cm
"""

# ╔═╡ 7c259d4c-122c-4270-ad64-185b09be4a2f
begin
    # generate uniform mesh, defined w/ size of domain and # cells in each direction
    Nx = round(Int, Lx / vox) # number of cells in x direction
    Ny = round(Int, Ly / vox) # numbers of cells in y direction    
    mesh = createMesh2D(Nx + 2, Ny + 2, Lx, Ly); # mesh

    # encode each cell location to its corresponding location in true space
    x_cell_reshape = [1:Nx...] * Lx / Nx
    y_cell_reshape = [1:Ny...] * Ly / Ny
end;

# ╔═╡ 29f239fd-446b-4e84-acc7-4075914c5383
begin
    # set up and evaluate transport model
    EpitaxialDeposition.PARAMS[:reactor][:wafer_diameter] = wafer_diam
    sol = Dict([key => trans_diff_Neumann(mesh, value, c₀[key], Nx) for (key, value) in D])
end;

# ╔═╡ 9e4a46ec-b4c0-4cd5-8b76-a7ce3bea88db
plot_species_profiles(sol, x_cell_reshape, y_cell_reshape)

# ╔═╡ Cell order:
# ╟─6584c779-4b89-40df-98c9-0dbf00ff6081
# ╠═d2eb8940-9b16-11ec-1d74-3bc00f5f950f
# ╟─6e4302dc-b2ce-4d8f-809b-078f2c39c2dc
# ╟─457a34ad-8514-4c6a-b16b-c94e6515ba40
# ╟─eb71d9ae-33d7-45dd-80ea-950a5a18ea78
# ╟─6e4c3904-6543-41e9-a862-22cb749ed178
# ╟─87620c55-7fbc-4004-928a-55c900a9a144
# ╟─00f98dbf-4df2-4972-8f3e-15425b735077
# ╟─9912d690-5dd6-45d9-94cf-8e0ca301fd29
# ╟─ac05f8e1-e284-48d0-8c81-c9b7677d8d52
# ╟─86c58901-b84d-45f7-a59e-521337b9ec36
# ╟─99be578d-673b-4f82-b0e4-2cd3fed8b847
# ╟─680313b8-d39c-428b-adb5-a10dee051860
# ╠═7c259d4c-122c-4270-ad64-185b09be4a2f
# ╠═29f239fd-446b-4e84-acc7-4075914c5383
# ╟─8075782d-3a72-4a43-833a-6f7a652f2cdf
# ╠═7abb2d83-c82e-4ad0-a4f6-0c9b67633f1d
# ╟─77168336-8588-4f53-a236-cc517978455c
# ╟─9e4a46ec-b4c0-4cd5-8b76-a7ce3bea88db
