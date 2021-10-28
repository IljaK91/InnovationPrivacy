@with_kw struct Pars
    # Economic Parameters
    θ       = 2
    α       = 0.33
    α_tilde = (θ - 1)/θ*α 
    δ       = 1.0
    γ       = 1 # adjustment cost parameter
    
    # Shocks
    a_t = 0
    σ²ₐ = 1
    a_n   = 5  # Number of nodes for productivity
    
    # Partial Equilibrium
    Y   = 1
    R   = 1.04

    # Linearly spaced grid
    K_min  = 0.01
    K_max  = 15
    K_n    = 20 # Number of nodes for capital
    K_grid = collect(range(K_min, K_max, length=K_n))
end

@with_kw struct BPars
    τ
    w_a # Bayesian weight on the private signal
    w_p # Bayesian weight on the prior
    σ²  # Posterior uncertainty
end
"""
    Parameter set for the replication of 'Firm Learning and Growth' (2018, RED)
"""
@with_kw struct ParsGrowth
    #! Deep Parameters
    σ     = 2.   # Elasticity of Substitution
    δ     = 0.1   # Probability of a firm dying
    f     = 0.1   # Fixed Cost of Operating
    β     = 0.96  # discount rate
    J     = 1.    # Mass of entering firms
    ξ     = 13.09 # Curvature exp(z) prod shock
    
    z_min = 1.    # z_min is not really important!
    #! Shocks
    θ_bar = 0     # mean demand shock
    σ²_θ  = 0.997 # variance of demand shocks

    σ²_ϵ  = 1.46    # variance of transitory demand shocks with zero mean
    ϵ_bar = -σ²_ϵ/(σ*2) # mean transient demand shock. Chosen to get a mean-preserving change in the variance.
    
    #! Steady State Parameters
    P = 1         # aggregate price level
    Y = 1         # aggregate expenditure
    w = 1         # wage rate

    #! Parameters for solving the model
    N            = 15 # Number of grid points
    a_min        = -4
    a_max        = 4
    a_grid       = range(a_min, a_max, length = N)
    a_thresholds = [(a_grid[i+1] - a_grid[i])/2 + a_grid[i] for i in 1:length(a_grid)-1]

    N_z    = 30 # Number of grid points
    z_max  = 3  # Check whether few firms become so large
    z_grid = range(z_min, z_max, length = N_z) # A grid for z, which we will not really use.

    u      = 0.5 # some random value
    u_grid = (exp.(z_grid).^(σ-1).*P.^(σ-1).*Y./(w.^(σ-1))).^(1/(σ-1)) # Some initial value, but will be replaced later on
    u_star = 0 # Needs to be solved for!
end

"""
    Bayesian parameter set for the replication of 'Firm Learning and Growth' (2018, RED)
"""
@with_kw struct BParsGrowth
    n
    V_n
    w_an
    w_pn
end