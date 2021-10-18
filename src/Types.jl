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
    δ     = 0.1  # Probability of a firm dying
    f     = 0.1  # Fixed Cost of Operating
    β     = 0.96 # discount rate

    #! Shocks
    θ_bar = 0 # mean demand shock
    σ²_θ  = 1 # variance of demand shocks

    σ²_ϵ  = 1 # variance of transitory demand shocks with zero mean
    
    z_bar = 0 # mean of firm productivity shocks
    σ²_z  = 1 # variance of firm productivity shocks

    #! Steady State Parameters
    P = 1 # aggregate price level
    Y = 1 # aggregate expenditure
    w = 1 # wage rate
end

"""
    Bayesian parameter set for the replication of 'Firm Learning and Growth' (2018, RED)
"""
@with_kw struct BParsGrowth
    V_n
    w_an
    w_pn
end