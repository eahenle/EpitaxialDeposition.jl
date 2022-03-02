module EpitaxialDeposition

using Reexport
@reexport using Catalyst, DifferentialEquations

include("params.jl")
include("reaction_model.jl")
include("misc.jl")

export PARAMS, deposition_rxn_network, odesys, film_thickness, estimate_derivative, K, K2, kKp

end
