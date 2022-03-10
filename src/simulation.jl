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
    return DataFrame(solve(ODEProblem(deposition_rxn_network, u₀, (0.0, t_max)), Tsit5(), saveat=PARAMS[:Δt][:Catalyst], maxiters=1e8)) 
end



function trans_diff_Neumann(m::MeshStructure, D::Float64, c₀::Float64, Nx::Int; N_steps=5, Lx=30., Ly=30.)
    # Define our boundary conditions. Package/function defaults yield Neumann boundary conditions, though Dirichlet and Robin conditions are also supported, and each cell face can have a unique boundary condition
    BC = createBC(m)

    # left and right boundary conditions for each coefficient in the boundary condition expression 
    
    #link: https://nbviewer.org/github/)simulkade/JFVM.jl/blob/master/examples/jfvm_tutorial.ipynb
    BC.left.a[:] .= 1.0 
    BC.right.a[:] .= 1.0

    BC.left.b[:] .= 0.0
    BC.right.b[:] .= 0.0

    BC.left.c[:] .= 0.0 
    BC.right.c[:] .= 0.0

    # top and bottom boundary conditions for each coefficient
    BC.top.a[1,:] .= 1.0 
    BC.bottom.a[1,:] .= 1.0

    BC.top.b[1,:] .= 0.0 
    BC.bottom.b[1,:] .= 0.0

    BC.top.c[1,:] .= 0.0 
    BC.bottom.c[1,:] .= 0.0

    BC.bottom.c[EpitaxialDeposition.wafer_boundary(Nx, Lx, EpitaxialDeposition.PARAMS[:reactor][:wafer_diameter]), 1] .= -1.0
        

    # Give a value for the diffusion coefficient based on current system
    D = D # [cm²/s]
    D_cell = createCellVariable(m, D) # assign the diffusion coefficient as a variable to each cell of the mesh
    D_face = geometricMean(D_cell) # choose an averaging scheme, for how the diffusion coefficient on a cell face is calculated
    
    c = createCellVariable(m, c₀, BC) # assign the "initial" species concentration to cells [mol/L]
    
    # Discretize the problem and build our solution
    M_diff = diffusionTerm(D_face) # matrix of diffusion term coefficients
    (M_bc, RHS_bc) = boundaryConditionTerm(BC) # matrix composed of coefficients and the right hand side for the BC



    ### modify to be 1/5th of the timestep from Catalyst
    
    if Lx < Ly # calculate a time step based on the smaller length (for more resolution)
    Δt = sqrt(Lx^2 / D) / N_steps # recommended time step calculation

    else
    Δt = sqrt(Ly^2 / D) / N_steps

    end

    # this code seems to be for taking discrete time steps forward, with size of Δt

    ### parameterize the # of steps?
    
    for i = 1:N_steps
        (M_t, RHS_t) = transientTerm(c, Δt, 1.0)
        M = M_t - M_diff + M_bc # adding together all sparse matrices of coefficients
        RHS = RHS_bc + RHS_t # add all RHS's to each other
        c = solveLinearPDE(m, M, RHS)
    end

return c
end
