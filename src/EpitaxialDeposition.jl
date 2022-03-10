module EpitaxialDeposition

using CairoMakie, DataFrames, JFVM, LsqFit, Memoize, Reexport
@reexport using Catalyst, DifferentialEquations

include("globals.jl")
include("Welty.jl")
include("diffusion.jl")
include("reaction_model.jl")
include("misc.jl")
include("simulation.jl")

export deposition_rxn_network, film_thickness, estimate_derivative, K, K2, kKp, run_simulation, Ω, DAB, plot_species_profiles, trans_diff_Neumann

end
