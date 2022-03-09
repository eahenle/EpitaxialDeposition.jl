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

# ╔═╡ 7c79cb6a-00cc-49d0-8110-bf849ae02c8c
begin
    # load project functions
    import Pkg
    Pkg.activate(".")
    using EpitaxialDeposition
    # load notebook utilities
    using CairoMakie, ColorSchemes, PlutoUI
end

# ╔═╡ 2ec89215-0634-498c-94ee-389b9b8d7034
md"""
# Epitaxial Silicon Film Growth
### Ian Harreschou, Adrian Henle
"""

# ╔═╡ 036cb0e3-e198-49f7-8221-3593fbeebb95
md"""
## Chemical/Physical Parameters
"""

# ╔═╡ 4842f48e-3d89-4104-b015-c383e0c63072
EpitaxialDeposition.PARAMS # from imported project code

# ╔═╡ d7f3703a-dee2-46ce-a052-92d44938c965
md"""
## Chemical Reaction Network
"""

# ╔═╡ 1b12f164-be6e-444c-bba8-c27b1cccb730
deposition_rxn_network

# ╔═╡ abb094af-f9c7-42a5-bdc7-d033ac9a5d72
md"""
## Equilibrium/Rate Expressions
"""

# ╔═╡ bdee4229-9e84-4bdd-a836-985a2f4b96b0
md"""
$K(T) = K_0 * \exp(-E/RT)$

$K_2(T) = K_{20} * \exp(-E_{a2}/RT)$

$K_3(T)=K_{30} * \exp(-E_{a3}/RT)$

$kK_p(T) = k * \exp(a-b/T)$
"""

# ╔═╡ 63d22599-3509-4027-b6bd-898dca165f65
md"""
## System of Differential Equations
"""

# ╔═╡ 05db1f52-1c87-465f-8b1b-25e45aed6d0b
# view the system of equations
EpitaxialDeposition.odesys

# ╔═╡ b7222499-de8e-452f-ab52-e1a443109147
md"""
## Simulation
"""

# ╔═╡ b9e889d9-5edb-498e-b560-fd699fe2de73
md"""
## Physical System

Wafer Diameter: $(@bind wafer_diameter PlutoUI.Slider(10.:45., default=30., show_value = true)) cm

Temperature: $(@bind temperature PlutoUI.Slider(900.:1400., default=1000., show_value = true)) K

Mole Fraction of SiCl₄: $(@bind y_SiCl₄ PlutoUI.Slider(0.:0.05:1., default=0.5, show_value = true)) 

Total System Pressure: $(@bind Psys PlutoUI.Slider(1.:1:760., default=760., show_value = true)) mmHg

Maximum Run Time: $(@bind t_max_min PlutoUI.Slider(1.:1:9600., default=10., show_value = true)) min
"""

# ╔═╡ 4016352d-ab26-4385-ba05-e57416df55c3
t_max = 60 * t_max_min;

# ╔═╡ 44189ce2-83df-4b77-8a5d-44383838d2ee
begin
    sol = run_simulation(y_SiCl₄, Psys, temperature, t_max)
    # calculate the film thickness from the simulation data
    δ = film_thickness.(sol[:, "Si_dep(t)"], sol[:, "Si_etch(t)"])
end;

# ╔═╡ 35d612fb-b859-42b1-9ea2-ae90ad01dfb8
md"""
## Result
"""

# ╔═╡ abccd9a2-0d1a-4562-89d7-3b0a0e5b2267
begin
    max_idx = argmax([isnan(d) ? 0 : d for d in δ])

md"""
**Maximum Film Thickness**: $(round(δ[max_idx], digits=2)) μm

**Time to Maximum Thickness**: $(sol[max_idx, :timestamp] / 60) min
"""
end

# ╔═╡ fd7dc5d2-28d1-4828-b53c-7a7f54fb4460
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
        lines!(sol[:, :timestamp], 1000 .* sol[:, name], label=chop(name, tail=3)) #mol/L --> mmol/L
    end

    fig[:,2] = Legend(fig, ax)

    Axis(
        fig[2,1],
        title="Si Film Growth on a $wafer_diameter cm Wafer, T = $temperature [K]",
        xlabel="Time [s]",
        ylabel="δ [μm]"
    )
    
    lines!(sol[:, :timestamp], δ)

    fig
end

# ╔═╡ e702092a-0b2e-475e-8dc1-51308be65c64
begin
    dδ = estimate_derivative(δ)
    
    local fig = Figure()
    local ax = Axis(
        fig[1,1],
        title = "Si Film Growth Rate over Time",
        ylabel = "δ Growth Rate [μm/s]",
        xlabel = "Time [s]"
    )

    lines!(ax, sol[:, :timestamp], dδ)

    idx = findfirst(dδ .<= 0.0) # find position, if any, where dδ = 0
    if ! isnothing(idx)
        vlines!(ax, sol[idx, :timestamp], linestyle = :dash, color = (:red,0.5))
    end
    
    fig
end

# ╔═╡ 23d028df-dbd7-4188-ae74-7da3c071e549
if ! isnothing(idx)
    
md"""
**Switched from Deposition to Etching at**: $(round(sol[idx, :timestamp], digits=2)) s
"""

end

# ╔═╡ 7b58b65f-a35a-461e-a4b6-4cb16646642e
md"
## Parametric Study #1: 
### Influence of Temperature on Film Thickness"

# ╔═╡ a347066c-c48e-42e9-a051-ac530d9bb7eb
md"Temperature Range: $(@bind temp_span PlutoUI.RangeSlider(900.:50:1400.)) K"

# ╔═╡ 96528e86-2a83-4b59-a412-de66b4475d52
begin
    local fig = Figure()
    Axis(
        fig[1,1],
        title="Parametric Study of Temperature on Film Thickness",
        xlabel="Time [s]",
        ylabel="δ [μm]"
    )

    local max_δ = zeros(length(temp_span))
    
    for (i, temp) in enumerate(temp_span)
        # run simulation w/ T=temp
        local sol = run_simulation(y_SiCl₄, Psys, temp, t_max)
        local δ = film_thickness.(sol[:, "Si_dep(t)"], sol[:, "Si_etch(t)"])

        non_negative = δ .> 0.0 # find all values in δ that aren't negative
        if length(δ[non_negative]) > 0 # avoid iteration over empty collection
            max_δ[i] = maximum(δ[non_negative]) # record max
        end

        lines!(
            sol[non_negative, :timestamp], 
            δ[non_negative], 
            color=ColorSchemes.viridis[round(Int, 256 * i / length(temp_span))]
        )
    end

    Colorbar(
        fig[1,2], 
        limits=(temp_span[1], temp_span[end]), 
        colormap=:viridis,
        label="Reactor Temperature [K]"
    )

    local ax2=Axis(
        fig[2,:],
        title="",
        xlabel="Temperature [K]",
        ylabel="Maximum δ [μm]"
    )

    scatter!(ax2, temp_span, max_δ)

    fig
end

# ╔═╡ 084b9228-5deb-4cef-bf6b-d962ca884dbe
md"
## Parametric Study #2: 
### Influence of Temperature on Reaction Rates"

# ╔═╡ 4e829cbd-3d8c-4ef6-8898-1af9443f01ac
begin
    local fig = Figure()
    local ax = Axis(
        fig[1,1],
        title="Forward Reaction Rate Constants",
        ylabel="Rate Constant Value [cm/s]",
        xlabel="Temperature [K]"
    )
    
    T = collect(range(temp_span[1], step=EpitaxialDeposition.PARAMS[:Δt][:Catalyst], temp_span[end]))
    
    lines!(ax, T, K.(T), label="SiCl₄ deposition")
    lines!(ax, T, K2.(T), label="HCl etching")

    fig[1,2] = Legend(fig, ax)

    local ax = Axis(
        fig[2,1],
        ylabel="Equilibrium Constant [∅]",
        xlabel="Temperature [K]"
    )

    lines!(ax, T, kKp.(T), label="SiCl₄ etching")

    fig[2,2] = Legend(fig, ax)
    
    fig
end

# ╔═╡ df6821d6-9c47-4d12-a62f-1b458c5ac6ba
md"
## Parametric Study #3: 
### Growth Rate as a Function of ``y_{SiCl_4}``"

# ╔═╡ ef49613a-22c3-4775-b903-665544c77bdf
begin
    local fig = Figure()
    local ax = Axis(
        fig[1,1],
        title="Growth Rate versus SiCl₄ Mole Fraction, T = $temperature [K]",
        ylabel="Film Growth Rate [μm/min]",
        xlabel="Mole Fraction of SiCl₄",
        xreversed=true
    )

    gaseous_species = Symbol.(["SiCl₄(t)", "SiCl₂(t)", "H₂(t)", "HCl(t)"])
    all_gases = map(r -> sum([r[x] for x in gaseous_species]), eachrow(sol))
    
    molfrac_SiCl4 = sol[:, gaseous_species[1]] ./ all_gases

    lines!(ax, molfrac_SiCl4, dδ .* 60)

    local idx = findfirst(dδ .< 0)
    if ! isnothing(idx)
        vlines!(ax, molfrac_SiCl4[idx], linestyle = :dash, color = (:red,0.5))
    end
    
    fig
end

# ╔═╡ Cell order:
# ╟─2ec89215-0634-498c-94ee-389b9b8d7034
# ╠═7c79cb6a-00cc-49d0-8110-bf849ae02c8c
# ╟─036cb0e3-e198-49f7-8221-3593fbeebb95
# ╠═4842f48e-3d89-4104-b015-c383e0c63072
# ╟─d7f3703a-dee2-46ce-a052-92d44938c965
# ╟─1b12f164-be6e-444c-bba8-c27b1cccb730
# ╟─abb094af-f9c7-42a5-bdc7-d033ac9a5d72
# ╟─bdee4229-9e84-4bdd-a836-985a2f4b96b0
# ╟─63d22599-3509-4027-b6bd-898dca165f65
# ╟─05db1f52-1c87-465f-8b1b-25e45aed6d0b
# ╟─b7222499-de8e-452f-ab52-e1a443109147
# ╠═44189ce2-83df-4b77-8a5d-44383838d2ee
# ╟─b9e889d9-5edb-498e-b560-fd699fe2de73
# ╟─4016352d-ab26-4385-ba05-e57416df55c3
# ╟─35d612fb-b859-42b1-9ea2-ae90ad01dfb8
# ╟─abccd9a2-0d1a-4562-89d7-3b0a0e5b2267
# ╟─fd7dc5d2-28d1-4828-b53c-7a7f54fb4460
# ╟─23d028df-dbd7-4188-ae74-7da3c071e549
# ╟─e702092a-0b2e-475e-8dc1-51308be65c64
# ╟─7b58b65f-a35a-461e-a4b6-4cb16646642e
# ╟─a347066c-c48e-42e9-a051-ac530d9bb7eb
# ╟─96528e86-2a83-4b59-a412-de66b4475d52
# ╟─084b9228-5deb-4cef-bf6b-d962ca884dbe
# ╟─4e829cbd-3d8c-4ef6-8898-1af9443f01ac
# ╟─df6821d6-9c47-4d12-a62f-1b458c5ac6ba
# ╟─ef49613a-22c3-4775-b903-665544c77bdf
