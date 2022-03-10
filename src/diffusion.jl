"""
Function to generate the interpolator for the collision integral Ω from the Welty data.
"""
function fit_Ω()
    # reformat Welty data
    welty = vcat(WELTY_DATA_RAW[:, [1, 3]], WELTY_DATA_RAW[:, [4,6]])
    # sort by x values
    sp = sortperm(welty[:, 1])
    # extract sorted values from relevant range
    xdata = welty[:, 1][sp][60:end, 1] # values of κT/ϵ
    ydata = welty[:, 2][sp][60:end, 1] # values of Ω
    # define model
    # @. makes all operations broadcasts (necessary for model fitting)
    @. model(x,p) = p[1] / (x * p[2] + p[3]) + p[4] # inverse x model
    # initial parameters
    p0 = 0.5 * ones(4)
    # fit the curve
    fit = curve_fit(model, xdata, ydata, p0);
    # generate the interpolator for the collision integral
    @memoize collision_interp(x) = model(x, fit.param)

    return collision_interp, xdata, ydata
end

# define the exported version of the collision integral interpolator
Ω, WELTY_X, WELTY_Y = fit_Ω()


"""
Calculates the diffusion coefficient for species A dissolved in species B
"""
@memoize function DAB(species1::Symbol, species2::Symbol; T=PARAMS[:reactor][:temperature], P=PARAMS[:reactor][:pressure])
    # calculate the collision diameter from critical parameters
    σ1 = 2.44 * (PARAMS[:Tc][species1] / PARAMS[:Pc][species1])^(1/3)
    σ2 = 2.44 * (PARAMS[:Tc][species2] / PARAMS[:Pc][species2])^(1/3)
    # calculate the ϵ's
    ϵ1 = 0.77 * κ * PARAMS[:Tc][species1]
    ϵ2 = 0.77 * κ * PARAMS[:Tc][species2]
    # combined parameters for calculating the diffusion coefficient
    σ12 = (σ1 + σ2) / 2
    ϵ12 = sqrt(ϵ1 * ϵ2)
    # calculate D
    return 0.001858 * T^(3/2) * sqrt(1 / PARAMS[:MW][species1] + 1 / PARAMS[:MW][species2]) / P / σ12^2 / Ω(κ * T / ϵ12) # [cm²/s]
end
