"""
    sol_df = run_simulation(y_SiCl₄, Psys, temperature, t_max)

Runs an epitaxial growth simulation with the specified parameters.
"""
function run_simulation(y_SiCl₄::Float64, Psys::Float64, temperature::Float64, t_max::Float64)::DataFrame
    # initial values
    local u₀ = [
        :SiCl₄   => (y_SiCl₄ * Psys) / PARAMS[:R][:L_Torr_mol_K] / temperature,
        :H₂      => Psys * (1 - y_SiCl₄) / PARAMS[:R][:L_Torr_mol_K] / temperature,
        :Si_dep  => 0.0,
        :Si_etch => 0.0,
        :HCl     => 0.0, 
        :SiCl₂   => 0.0,
        :k       => PARAMS[:kKp][:k],
        :T       => temperature
    ]

    # define the ODEs and solve
    return DataFrame(solve(ODEProblem(deposition_rxn_network, u₀, (0.0, t_max)), Tsit5(), saveat=PARAMS[:Δt], maxiters=1e8))
end
