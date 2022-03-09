"""
    δ = film_thickness(dep, etch)

Calacultes the thickness of the depositied film given the pseudo-concentrations of the deposited and etched material.
Takes optional parameters `V`, `d`, `ρ`, and `MW`; by default, these are pulled from `PARAMS` for the default Si deposition system.
"""
function film_thickness(dep::Float64, etch::Float64; 
        V=PARAMS[:reactor][:volume]::Float64, d=PARAMS[:reactor][:wafer_diameter]::Float64, ρ::Float64=PARAMS[:ρ][:Si], MW::Float64=PARAMS[:MW][:Si])::Float64
    net      = dep - etch   # net pseudo-concentration of deposited species
    net_mol  = net * V      # net moles deposited species
    dep_mass = net_mol * MW # net mass deposited species (g)
    dep_vol  = dep_mass / ρ # net volume deposited species (cm³)
    return dep_vol / (π * d^2 / 4) * 1e4 # film thickness (μm)
end

"""
    dδ_dt = estimate_derivative(δ)

Computes an estimate of the derivative at each point in the vector `x`.
The final point in the vector is assumed to have the same derivative as the point immediately prior to it.
"""
function estimate_derivative(x::Vector{Float64})::Vector{Float64}
    dδ = zeros(length(x))
    @inbounds for i = 2:length(x)
        dδ[i-1] = (x[i] - x[i-1]) ./ PARAMS[:Δt][:Catalyst]
    end
    dδ[end] = dδ[end-1]
    return dδ
end


@memoize function wafer_boundary(Nx, Lx, wafer_diameter)  
    @assert wafer_diameter ≤ Lx "Wafer too big! d = $wafer_diameter, Lx = $Lx"
    half = ceil(Nx / 2 + 1)
    dticks = wafer_diameter * Nx / Lx 
    left, right = half - dticks/2, half + dticks/2
    return round(Int, left):round(Int, right)
end


function plot_species_profiles(sol, x, y)
    fig = CairoMakie.Figure()

    Axis(fig[1,1], title="H₂")
    heatmap!(x, y, sol[:H₂].value)

    Axis(fig[1,2], title="SiCl₄")
    heatmap!(x, y, sol[:SiCl₄].value)

    Axis(fig[2,1], title="HCl")
    heatmap!(x, y, sol[:HCl].value)

    Axis(fig[2,2], title="SiCl₂")
    heatmap!(x, y, sol[:SiCl₂].value)

    Colorbar(fig[:,3], limits = (0, maximum(sol[:H₂].value)), colormap = :viridis,
            label = "$species Concentration [mmol/L]")

    return fig
end
