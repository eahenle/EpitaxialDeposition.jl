PARAMS = Dict{Symbol, Union{Dict{Symbol, Float64}, Float64}}(
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
    :Si => Dict{Symbol, Float64}( # silicon physical properties
        :ρ => 2.33,
        :m => 28.086
    ),
    :kKp => Dict{Symbol, Float64}( # equilibrium expression parameters
        :k  => 1.,
        :Ea => 10.38,
        :b  => 16770.
    ),
    :reactor => Dict{Symbol, Float64}( # reactor parameters
        :volume         => 10.,
        :wafer_diameter => 30.
    ),
    :Δt => 0.05 # simulation time step, also used in calculating time-derivative of film thickness
)
