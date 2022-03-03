# define the temperature-dependent equilibrium constant expressions
K(T::Float64)::Float64   = π * (PARAMS[:reactor][:wafer_diameter] / 2)^2 * PARAMS[:rate][:K0]  * exp(-PARAMS[:energy][:E]   / (PARAMS[:R][:kcal_mol_K] * T)) # [cm/s]
K2(T::Float64)::Float64  = π * (PARAMS[:reactor][:wafer_diameter] / 2)^2 * PARAMS[:rate][:K20] * exp(-PARAMS[:energy][:Ea2] / (PARAMS[:R][:kcal_mol_K] * T)) # [cm/s]
K3(T::Float64)::Float64  = π * (PARAMS[:reactor][:wafer_diameter] / 2)^2 * PARAMS[:rate][:K30] * exp(-PARAMS[:energy][:Ea3] / (PARAMS[:R][:kcal_mol_K] * T)) # [cm/s]
kKp(T::Float64)::Float64 = π * (PARAMS[:reactor][:wafer_diameter] / 2)^2 * PARAMS[:kKp][:k] * 10^(PARAMS[:kKp][:Ea] - PARAMS[:kKp][:b] / T) # [∅]

# register the equilibrium constant expressions with Catalyst
@register_symbolic kKp(T)
@register_symbolic K(T)
@register_symbolic K2(T) 
@register_symbolic K3(T)

# define the reaction network for the deposition process with competitive etching
deposition_rxn_network = @reaction_network begin
    K(T), SiCl₄ + 2H₂ --> Si_dep + 4HCl     # Si deposition from source gas
    K2(T), 2HCl --> SiCl₂ + H₂ + Si_etch    # Si etching by HCl
    kKp(T), SiCl₄ --> 2SiCl₂ + Si_etch      # Si etching by SiCl₄
    k, 2SiCl₂ --> SiCl₄ + Si_dep            # Si deposition by SiCl₂ decomposition
end

# generate the system of ODEs
odesys = convert(ODESystem, deposition_rxn_network, combinatoric_ratelaws=false)
