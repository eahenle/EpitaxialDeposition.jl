module EpitaxialDeposition

using DataFrames, Reexport
@reexport using Catalyst, DifferentialEquations

include("params.jl")
include("reaction_model.jl")
include("misc.jl")
include("simulation.jl")

export deposition_rxn_network, film_thickness, estimate_derivative, K, K2, kKp, run_simulation, PARAMS

end
