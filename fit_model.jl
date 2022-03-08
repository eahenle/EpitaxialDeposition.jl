### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 09961cea-99f6-11ec-279b-79c147b6e3b4
begin
    import Pkg
    Pkg.activate(".")
    using CSV, DataFrames, EpitaxialDeposition, LsqFit
end

# ╔═╡ 95fbb48d-6339-46d3-9270-6bf88fc00c20
df = CSV.read("Steinmaier_1423K.csv", DataFrame)

# ╔═╡ 3687dea3-3ad6-42fe-9aaa-6ef9626647a5
begin
    P = 100.
    T = 1423.
    t_max = 3600.

    gaseous_species = Symbol.(["SiCl₄(t)", "SiCl₂(t)", "H₂(t)", "HCl(t)"])
end;

# ╔═╡ a5afed75-96dd-42c5-aaea-186404161489
function model(y, p)
    # set model parameters
    PARAMS[:kKp][:k] = p[1]
    PARAMS[:kKp][:a] = p[2]
    PARAMS[:kKp][:b] = p[3]
	
    # solve ODE system and find film growth rates
    sol = [run_simulation(y, P, T, t_max) for y in y]
    δ = [film_thickness.(sol[:, "Si_dep(t)"], sol[:, "Si_etch(t)"]) for sol in sol]
    dδ = [estimate_derivative(δ) for δ in δ]

    # calculate mole fractions
    all_gases = 
		[map(r -> sum([r[x] for x in gaseous_species]), eachrow(sol)) for sol in sol]
    molfrac_SiCl4 = 
		[sol[:, gaseous_species[1]] ./ all_gases[i] for (i, sol) in enumerate(sol)]
    # find index of timepoint w/ mole fraction ~y
    idx = 
		[argmin(abs.(molfrac_SiCl4 .- y[i])) for (i, molfrac_SiCl4) in enumerate(molfrac_SiCl4)]
    # return the film growth rate for the given mole fraction
    return [dδ[i][idx] for (i, idx) in enumerate(idx)]
end;

# ╔═╡ 50b5383c-1f2d-4d95-a33f-760664412bfb
fit = curve_fit(model, df.y, df.r, 
	[
		PARAMS[:kKp][:k],
		PARAMS[:kKp][:Ea],
		PARAMS[:kKp][:b]
	]
)

# ╔═╡ 95b818ba-dfdf-40b6-86e0-4afd77f20db5
fit.param

# ╔═╡ Cell order:
# ╠═09961cea-99f6-11ec-279b-79c147b6e3b4
# ╠═95fbb48d-6339-46d3-9270-6bf88fc00c20
# ╠═3687dea3-3ad6-42fe-9aaa-6ef9626647a5
# ╠═a5afed75-96dd-42c5-aaea-186404161489
# ╠═50b5383c-1f2d-4d95-a33f-760664412bfb
# ╠═95b818ba-dfdf-40b6-86e0-4afd77f20db5
