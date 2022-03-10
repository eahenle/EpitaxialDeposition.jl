# Load Lennard-Jones parameters from Welty, et al.
include("Welty.jl")


# Global constant for diffusion calculations
const κ = 1.38e-16 # ergs/K


"""
Global parameters for simulations.  Includes:
    - rate constants
    - activation energies
    - ideal gas constants (different units)
    - densities
    - molar masses
    - equilibrium expression parameters
    - physical system parameters
    - timesteps (different submodules)
    - critical point data 
"""
PARAMS = Dict{Symbol, Dict{Symbol, Float64}}(
    :rate => Dict{Symbol, Float64}( # constants for rate expressions
        :K0     => 2950.94,
        :K10    => 1.899e-16,
        :K20    => 39.3886,
        :K30    => 952.058
    ),
    :energy => Dict{Symbol, Float64}( # activation energies
        :E      => 18.910,
        :Ea1    => 10.620,
        :Ea2    => 17.490,
        :Ea3    => 21.721
    ),
    :R => Dict{Symbol, Float64}( # ideal gas constant
        :kcal_mol_K     => 1.985e-3,
        :L_Torr_mol_K   => 62.3636
    ),
    :ρ => Dict{Symbol, Float64}( # densities, g/cm³
        :Si => 2.33,
    ),
    :MW => Dict{Symbol, Float64}( # molar masses, g/mol
        :Si     => 28.086,
        :SiCl₄  => 169.9, 
        :H₂     => 2.016, 
        :HCl    => 36.46, 
        :SiCl₂  => 98.9915
    ),
    :kKp => Dict{Symbol, Float64}( # equilibrium expression parameters
        :k  => 1.,
        :Ea => 10.38,
        :b  => 16770.
    ),
    :reactor => Dict{Symbol, Float64}( # reactor parameters
        :volume         => 10.,     # L
        :wafer_diameter => 30.,     # cm
        :pressure       => 1.,      # atm
        :temperature    => 1100.    # K
    ),
    :Δt => Dict{Symbol, Float64}( # simulation time steps
        :Catalyst   => 0.05, # s
        # JFVM should always be 1/5th of Catalyst, b/c we do 5 transport steps per kinetic step
        :JFVM       => 0.01  # s
    ),
    :Tc => Dict{Symbol, Float64}( # critical temperature, K
        :HCl    => 324.68,
        :SiCl₄  => 508.1, 
        :H₂     => 33.18, 
        :SiCl₂  => 452.05 # actually for SiH₂Cl₂, b/c unknown for SiCl₂
    ),
    :Pc => Dict{Symbol, Float64}( # critical pressure, atm
        :HCl    => 81.5, 
        :SiCl₄  => 35.469, 
        :H₂     => 12.833, 
        :SiCl₂  => 44. # actually for SiH₂Cl₂, b/c unknown for SiCl₂
    )
)
