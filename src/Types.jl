@with_kw struct Pars
    # Economic Parameters
    θ = 2
    α = 0.33
    α_tilde = (θ - 1) / θ * α
    δ = 1.0
    γ = 1 # adjustment cost parameter

    # Shocks
    a_t = 0
    σ²ₐ = 1
    a_n = 5  # Number of nodes for productivity

    # Partial Equilibrium
    Y = 1
    R = 1.04

    # Linearly spaced grid
    K_min = 0.01
    K_max = 15
    K_n = 20 # Number of nodes for capital
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
    σ = 2.0   # Elasticity of Substitution
    δ = 0.1   # Probability of a firm dying
    f = 0.1   # Fixed Cost of Operating
    β = 0.96  # discount rate
    J = 1.0    # Mass of entering firms
    ξ = 13.09 # Curvature exp(z) prod shock

    z_min = 1.0    # z_min is not really important!

    #! Shocks
    θ_bar = 0     # mean demand shock
    σ²_θ = 0.997 # variance of demand shocks
    ρ = 0.9   # persistence of persistent demand

    σ²_ϵ = 1.46    # variance of transitory demand shocks with zero mean
    ϵ_bar = -σ²_ϵ / (σ * 2) # mean transient demand shock. Chosen to get a mean-preserving change in the variance.

    σ²_η = 0.3 # Variance of the innovation in persistent demand

    #! Steady State Parameters
    P = 1         # aggregate price level
    Y = 1         # aggregate expenditure
    w = 1         # wage rate
    L = 1         # Labor supply

    #! Parameters for solving the model
    N = 20 # Number of grid points
    a_min = -8
    a_max = 8
    a_grid = range(a_min, a_max, length=N)
    a_thresholds = [(a_grid[i+1] - a_grid[i]) / 2 + a_grid[i] for i = 1:length(a_grid)-1]

    N_z = 20 # Number of grid points
    z_max = 7  # Check whether few firms become so large
    z_grid = range(z_min, z_max, length=N_z) # A grid for z, which we will not really use.

    u = 0.5 # some random value
    u_grid = (exp.(z_grid) .^ (σ - 1) .* P .^ (σ - 1) .* Y ./ (w .^ (σ - 1))) .^ (1 / (σ - 1)) # Some initial value, but will be replaced later on
    u_star = 0 # Needs to be solved for!

    #! Depreciation
    unc_max = 10
    unc_min = 0.5
    uncertainty_grid = range(unc_min, unc_max, length=N_unc)
end

@with_kw struct StSt
    M # Mass of Firms
    Y # Output
    P # Price Level
    E # Entry/Exit
end

"""
    Bayesian parameter set for the replication of 'Firm Learning and Growth' (2018, RED)
"""
@with_kw struct BParsGrowth
    n     # Number of periods
    V_n   # Posterior uncertainty after n periods
    w_an  # Weight on new signal
    w_pn  # Weight on prior
end

@with_kw struct ParsV2
    #! Deep Parameters
    σ = 2.0   # Elasticity of Substitution
    δ = 0.1   # Probability of a firm dying
    f = 0.1   # Fixed Cost of Operating
    β = 0.96  # discount rate
    J = 1.0   # Mass of entering firms
    ξ = 13.09 # Curvature exp(z) prod shock

    z_min = 1.0 # Minimum firm productivity

    #! Shocks
    #θ_bar = 0     # mean demand shock
    #σ²_θ  = 0.997 # variance of demand shocks
    ρ = 0.99  # persistence of persistent demand

    σ²_ϵ = 1.46    # variance of transitory demand shocks with zero mean
    ϵ_bar = -σ²_ϵ / (σ * 2) # mean transient demand shock. Chosen to get a mean-preserving change in the variance.

    σ²_η = 0.3 # Variance of the innovation in persistent demand

    σ²_ν = 1 # Scaling variance of the signal

    #! Steady State Parameters
    P = 1         # aggregate price level
    Y = 1         # aggregate expenditure
    w = 1         # wage rate
    L = 1         # Labor supply

    #! Parameters for solving the model
    # Prior mean
    N_μ = 20 # Number of grid points
    μ_min = -8
    μ_max = 8
    μ_grid = range(μ_min, μ_max, length=N_μ)

    # To get the cut-off points for the transition matrix
    μ_thresholds = [(μ_grid[i+1] + μ_grid[i]) / 2 for i = 1:length(μ_grid)-1]

    N_V = 20
    V_min = 0.5
    V_max = 10
    V_grid = range(V_min, V_max, length=N_V)

    N_z = 20 # Number of grid points
    z_max = 7  # Check whether few firms become so large
    z_grid = range(z_min, z_max, length=N_z) # A grid for z, which we will not really use.
    #! This is from the Learning and Firm Growth Paper.
    #! Basically, some part of firm profits is fixed and can be solved for later.
    u = 0.5 # some random value
    u_grid = (exp.(z_grid) .^ (σ - 1) .* P .^ (σ - 1) .* Y ./ (w .^ (σ - 1))) .^ (1 / (σ - 1)) # Some initial value, but will be replaced later on
    u_star = 0 # Needs to be solved for!
end


@with_kw struct Pars_Bundle
    #! Deep Parameters
    σ = 2.0  # Elasticity of Substitution
    v = 1.0   # Returns to Scale
    γ = 0.5  # Relative importance of goods
    α = 0.3  # labor share
    δ = 0.1  # intangible capital share
    τ = 0.1  # Icerberg Transportation Cost
    ϕ = 0.0   # How much capital is kept by the firm that is selling capital
    α_Y = (σ * v - σ + 1) / (σ * v)

    α_hat = (σ - 1) / σ * α
    δ_hat = (σ - 1) / σ * δ

    ϵ = 2 # elasticity of substitution between capital varieties
    ξ = 0.5 # productivity of exogenous capital

    #! Steady State Parameters
    P = 1       # aggregate price level
    Y = 0.95    # aggregate output
    w = 1       # wage rate
    L = 1       # Labor supply
    p_I_A = 1       # Price of Intangible Capital of A firms
    p_I_B = 1       # Price of Intangible Capital of B firms

    #! Solution Firm-Parameters
    l_A_ss = 0
    l_B_ss = 0
    k_I_P_A_ss = 0
    k_I_P_B_ss = 0
    k_I_b_B_ss = 0
    k_I_b_A_ss = 0
    k_I_B_ss = 0
    k_I_A_ss = 0

    # More robust to calculate it once and save it than calculating it multiple times.
    #! Cost parameters
    # Multiplicative Parameter
    ca_A = 0.1
    ca_B = 0.3

    # Curvature parameter
    cb_A = 2
    cb_B = 2


end
@with_kw struct Pars_Simple
    #! Deep Parameters
    σ = 2.0  # Elasticity of Substitution
    v = 1.0   # Returns to Scale
    γ = 0.5  # Relative importance of goods
    α = 0.3  # labor share
    δ = 0.1  # intangible capital share
    τ = 0.1  # Icerberg Transportation Cost
    ϕ = 0.0   # How much capital is kept by the firm that is selling capital
    α_C = (σ * v - σ + 1) / (σ * v)

    α_hat = (σ - 1) / σ * α
    δ_hat = (σ - 1) / σ * δ

    ϵ = 2 # elasticity of substitution between capital varieties
    ξ = 0.5 # productivity of exogenous capital

    #! Steady State Parameters
    P = 1       # aggregate price level
    C = 0.95    # aggregate output
    w = 1       # wage rate
    L = 1       # Labor supply
    p_I = 1       # Price of Intangible Capital

    #! Solution Firm-Parameters
    l_A_ss = 0
    l_B_ss = 0
    k_I_P_A_ss = 0
    k_I_P_B_ss = 0
    k_I_b_B_ss = 0
    k_I_b_A_ss = 0
    k_I_B_ss = 0
    k_I_A_ss = 0

    # More robust to calculate it once and save it than calculating it multiple times.
    #! Cost parameters
    # Multiplicative Parameter
    ca_A = 0.1
    ca_B = 0.3

    # Curvature parameter
    cb_A = 2
    cb_B = 2


end
@with_kw struct Pars_OLG
    #! Deep Parameters
    σ = 2.0  # Elasticity of Substitution
    γ = 0.5  # Share of Sophisticated Firms
    α = 0.3  # labor share
    δ = 0.1  # intangible capital share
    τ = 0.1  # Icerberg Transportation Cost
    β = 0.9 # Discount Rate
    α_Y = 1 / σ
    μ = 1.0

    α_hat = (σ - 1) / σ * α
    δ_hat = (σ - 1) / σ * δ

    f = 1.5 # fixed cost

    #! Steady State Parameters
    P = 1       # aggregate price level
    Y = 0.95    # aggregate output
    w = 1       # wage rate
    L = 1       # Labor supply
    p_I = 1     # Price of Intangible Capital

    #! Solution Firm-Parameters
    l_S_ss = 0
    l_U_ss = 0
    l_Y_ss = 0
    k_P_S_ss = 0
    k_P_U_ss = 0
    k_T_U_ss = 0
    k_T_S_ss = 0
    k_T_Y_ss = 0
    k_U_ss = 0
    k_S_ss = 0
    k_Y_ss = 0

    # More robust to calculate it once and save it than calculating it multiple times.
    #! Cost parameters
    # cost = a * k ^ b
    # Multiplicative Parameter
    ca_S = 0.3
    ca_U = 5

    # Curvature parameter
    cb_S = 2
    cb_U = 2
end
# """
#     The parameters and steady state values for the model written down by Ilja and Roxana.
#     Main features: Firms share structured which still needs to be processed internally to generate knowledge.
#     Additionally, data is partly non-rival (ν > 0). The rival-formulation is nested

# """
# @with_kw struct Pars_v3
#     #! Deep Parameters
#     σ = 2.0  # Elasticity of Substitution
#     c = ((σ - 1) / σ) # auxilary variable
#     w_S = 0.5  # Share of Sophisticated Firms
#     w_B = 0.5  # Share of Firms with large customer base.
#     α = 0.3  # labor share, knowledge share is 1-α
#     ϕ = 0.3  # diminishing returns data generation
#     τ = 0.1  # icerberg transportation cost
#     τ2 = 0.1  # icerberg transportation cost
#     ν = 0.05 # non-rivalry parameter
#     ε = 2    # elasticity of substition data bundle
#     ξ = 1.0  # relative usefulness external data
#     β = 0.94 # Discount Factor
#     ζ = 1.0  #! Does not work and needs to be equal to 1 at the moment!!!!

#     compl_prod::Symbol = :yes

#     α_Y = set_a_Y(compl_prod, σ, ζ) # auxiliary paramter


#     α_L_hat = (σ - 1) / σ * (1 - α) # adjusted labor factor share 
#     α_K_hat = (σ - 1) / σ * α # adjusted knowledge factor share

#     bundle = :yes # Switch in case I want to also implement a version without the data bundle, but making ε very large is similar

#     λ = 0.1234

#     #! Persistent Knowledge and Data
#     δ_D = 0.1 # depreciation rate for data
#     δ_K = 0.1 # depreciation rate for knowledge
#     second_period = :yes

#     #! Steady State Parameters
#     P = 1  # aggregate price level
#     #! This is the solution for the model with tau = 0, nu = 0, α_Y = 0
#     Y = 2.075300672340325  # aggregate output
#     w = 0.2367874818108683  # wage rate good production
#     w_P = 0.04669426512765732  # wage rate data processing
#     w_G = 0.06021347283561476  # wage rate data generation
#     L = 1  # Labor supply for production of the intermediate good
#     L_P = 1  # Labor supply for processing data
#     L_G = 1  # Labor supply for generating / structuring data
#     p_D = 0.11356059692613776  # Price of Structured Data
#     Ω = 0.0 # data multiplier, generated minus iceberg transportation plus non-rival
#     Ω_bundle = 0.0 # data multiplier, bundle relative to generated

#     Y2 = 2.075300672340325  # aggregate output
#     w2 = 0.2367874818108683  # wage rate good production
#     w_P2 = 0.04669426512765732  # wage rate data processing
#     w_G2 = 0.06021347283561476  # wage rate data generation
#     L2 = 1  # Labor supply for production of the intermediate good
#     L_P2 = 1  # Labor supply for processing data
#     L_G2 = 1  # Labor supply for generating / structuring data
#     p_D2 = 0.11356059692613776  # Price of Structured Data
#     Ω2 = 0.0 # data multiplier, generated minus iceberg transportation plus non-rival
#     Ω_bundle2 = 0.0 # data multiplier, bundle relative to generated

#     D_G_SS = 0.0 # Total data generated
#     D_S_SS = 0.0 # Total data shared

#     D_G_SS2 = 0.0 # Total data generated
#     D_S_SS2 = 0.0 # Total data shared

#     #! Firm Specific Parameters
#     A_G_common = 1.0
#     A_G_unsoph = 0.2
#     A_P_common = 1.0
#     A_P_unsoph = 1.0

#     A_G_common2 = 1.0
#     A_P_common2 = 1.0

#     A_G_S = A_G_unsoph * A_G_common # data generation productivity parameter, small
#     A_G_B = 1.0 * A_G_common # data generation productivity parameter, big
#     A_P_U = A_P_unsoph * A_P_common # data processing productivity parameter, unsophisticated
#     A_P_S = 1.0 * A_P_common # data processing productivity parameter, sophisticated

#     γ_S = 0.2  # labor share data processing sophisticated, data share is 1-γ
#     γ_U = 0.4  # labor share data processing unsophisticated, data share is 1-γ

#     Firm_sol = 0

#     #! Solution Firm-Parameters
#     l_LU_ss = 0 # Labor good production small customer base, unsophisticated
#     l_BU_ss = 0 # Labor good production big customer base, unsophisticated
#     l_LS_ss = 0 # Labor good production small customer base, sophisticated
#     l_BS_ss = 0 # Labor good production big customer base, sophisticated

#     l_P_LU_ss = 0 # Labor good production small customer base, unsophisticated
#     l_P_BU_ss = 0 # Labor good production big customer base, unsophisticated
#     l_P_LS_ss = 0 # Labor good production small customer base, sophisticated
#     l_P_BS_ss = 0 # Labor good production big customer base, sophisticated

#     l_G_LU_ss = 0 # Labor good production small customer base, unsophisticated
#     l_G_BU_ss = 0 # Labor good production big customer base, unsophisticated
#     l_G_LS_ss = 0 # Labor good production small customer base, sophisticated
#     l_G_BS_ss = 0 # Labor good production big customer base, sophisticated

#     D_S_LU_ss = 0 # data sharing small customer base, unsophisticated
#     D_S_BU_ss = 0 # data sharing big customer base, unsophisticated
#     D_S_LS_ss = 0 # data sharing small customer base, sophisticated
#     D_S_BS_ss = 0 # data sharing big customer base, sophisticated

#     D_G_LU_ss = 0 # data generated small customer base, unsophisticated
#     D_G_BU_ss = 0 # data generated big customer base, unsophisticated
#     D_G_LS_ss = 0 # data generated small customer base, sophisticated
#     D_G_BS_ss = 0 # data generated big customer base, sophisticated

#     D_I_LU_ss = 0 # data generated small customer base, unsophisticated
#     D_I_BU_ss = 0 # data generated big customer base, unsophisticated
#     D_I_LS_ss = 0 # data generated small customer base, sophisticated
#     D_I_BS_ss = 0 # data generated big customer base, sophisticated

#     D_E_LU_ss = 0 # data generated small customer base, unsophisticated
#     D_E_BU_ss = 0 # data generated big customer base, unsophisticated
#     D_E_LS_ss = 0 # data generated small customer base, sophisticated
#     D_E_BS_ss = 0 # data generated big customer base, sophisticated

#     D_LU_ss = 0 # data generated small customer base, unsophisticated
#     D_BU_ss = 0 # data generated big customer base, unsophisticated
#     D_LS_ss = 0 # data generated small customer base, sophisticated
#     D_BS_ss = 0 # data generated big customer base, sophisticated

#     D_last_LU = 0 # data generated small customer base, unsophisticated
#     D_last_BU = 0 # data generated big customer base, unsophisticated
#     D_last_LS = 0 # data generated small customer base, sophisticated
#     D_last_BS = 0 # data generated big customer base, sophisticated

#     K_LU_ss = 0 # knowledge employed small customer base, unsophisticated
#     K_BU_ss = 0 # knowledge employed big customer base, unsophisticated
#     K_LS_ss = 0 # knowledge employed small customer base, sophisticated
#     K_BS_ss = 0 # knowledge employed big customer base, sophisticated

#     K_last_LU = 0 # knowledge employed small customer base, unsophisticated
#     K_last_BU = 0 # knowledge employed big customer base, unsophisticated
#     K_last_LS = 0 # knowledge employed small customer base, sophisticated
#     K_last_BS = 0 # knowledge employed big customer base, sophisticated

#     prof_LU_ss = 0 # knowledge employed small customer base, unsophisticated
#     prof_BU_ss = 0 # knowledge employed big customer base, unsophisticated
#     prof_LS_ss = 0 # knowledge employed small customer base, sophisticated
#     prof_BS_ss = 0 # knowledge employed big customer base, sophisticated
# end

@with_kw struct Pars_v3
    #! Deep Parameters
    σ = 2.0  # Elasticity of Substitution
    c = (σ - 1) / σ # auxilary variable
    w_S = 0.5  # Share of Sophisticated Firms
    w_B = 0.5  # Share of Firms with large customer base.
    α = 0.3  # labor share, knowledge share is 1-α
    ϕ = 0.3  # diminishing returns data generation
    τ = 0.1  # icerberg transportation cost
    ν = 0.05 # non-rivalry parameter
    ε = 2    # elasticity of substition data bundle
    ξ = 1.0  # relative usefulness external data
    β = 0.94 # Discount Factor
    ζ = 1.0  #! Does not work and needs to be equal to 1 at the moment!!!!

    α_S = 0.8  # data share data processing sophisticated, data share is 1-γ
    α_U = 0.6  # data share data processing unsophisticated, data share is 1-γ
    
    compl_prod::Symbol = :yes

    α_Y = set_a_Y(compl_prod, σ, ζ) # auxiliary paramter

    α_S_L_hat = (σ - 1) / σ * (1 - α_S) # adjusted labor factor share sophisticated
    α_S_K_hat = (σ - 1) / σ * α_S # adjusted knowledge factor share sophisticated

    α_U_L_hat = (σ - 1) / σ * (1 - α_U) # adjusted labor factor share unsophisticated
    α_U_K_hat = (σ - 1) / σ * α_U # adjusted knowledge factor share unsophisticated

    bundle = :yes # Switch in case I want to also implement a version without the data bundle, but making ε very large is similar

    #! Rate of convergence
    λ = 0.1234

    #! Steady State Parameters
    P = 1  # aggregate price level
    #! This is the solution for the model with tau = 0, nu = 0, α_Y = 0
    Y = 2.075300672340325  # aggregate output
    w = 0.2367874818108683  # wage rate good production
    w_G = 0.06021347283561476  # wage rate data generation
    L = 1  # Labor supply for production of the intermediate good
    L_G = 1  # Labor supply for generating / structuring data
    p_D = 0.11356059692613776  # Price of Structured Data
    Ω = 0.0 # data multiplier, generated minus iceberg transportation plus non-rival
    Ω_bundle = 0.0 # data multiplier, bundle relative to generated

    D_G_SS = 0.0 # Total data generated
    D_S_SS = 0.0 # Total data shared


    #! Firm Specific Parameters
    A_G_common = 1.0
    A_G_unsoph = 0.2
    A_P_common = 1.0
    A_P_unsoph = 1.0

    A_G_S = A_G_unsoph * A_G_common # data generation productivity parameter, small
    A_G_B = 1.0 * A_G_common # data generation productivity parameter, big
    A_P_U = A_P_unsoph * A_P_common # data processing productivity parameter, unsophisticated
    A_P_S = 1.0 * A_P_common # data processing productivity parameter, sophisticated

    Firm_sol = 0

    #! Solution Firm-Parameters
    l_LU_ss = 0 # Labor good production small customer base, unsophisticated
    l_BU_ss = 0 # Labor good production big customer base, unsophisticated
    l_LS_ss = 0 # Labor good production small customer base, sophisticated
    l_BS_ss = 0 # Labor good production big customer base, sophisticated

    l_G_LU_ss = 0 # Labor good production small customer base, unsophisticated
    l_G_BU_ss = 0 # Labor good production big customer base, unsophisticated
    l_G_LS_ss = 0 # Labor good production small customer base, sophisticated
    l_G_BS_ss = 0 # Labor good production big customer base, sophisticated

    D_S_LU_ss = 0 # data sharing small customer base, unsophisticated
    D_S_BU_ss = 0 # data sharing big customer base, unsophisticated
    D_S_LS_ss = 0 # data sharing small customer base, sophisticated
    D_S_BS_ss = 0 # data sharing big customer base, sophisticated

    D_G_LU_ss = 0 # data generated small customer base, unsophisticated
    D_G_BU_ss = 0 # data generated big customer base, unsophisticated
    D_G_LS_ss = 0 # data generated small customer base, sophisticated
    D_G_BS_ss = 0 # data generated big customer base, sophisticated

    D_I_LU_ss = 0 # data generated small customer base, unsophisticated
    D_I_BU_ss = 0 # data generated big customer base, unsophisticated
    D_I_LS_ss = 0 # data generated small customer base, sophisticated
    D_I_BS_ss = 0 # data generated big customer base, sophisticated

    D_E_LU_ss = 0 # data generated small customer base, unsophisticated
    D_E_BU_ss = 0 # data generated big customer base, unsophisticated
    D_E_LS_ss = 0 # data generated small customer base, sophisticated
    D_E_BS_ss = 0 # data generated big customer base, sophisticated

    D_LU_ss = 0 # data generated small customer base, unsophisticated
    D_BU_ss = 0 # data generated big customer base, unsophisticated
    D_LS_ss = 0 # data generated small customer base, sophisticated
    D_BS_ss = 0 # data generated big customer base, sophisticated

    K_LU_ss = 0 # knowledge employed small customer base, unsophisticated
    K_BU_ss = 0 # knowledge employed big customer base, unsophisticated
    K_LS_ss = 0 # knowledge employed small customer base, sophisticated
    K_BS_ss = 0 # knowledge employed big customer base, sophisticated

    prof_LU_ss = 0 # knowledge employed small customer base, unsophisticated
    prof_BU_ss = 0 # knowledge employed big customer base, unsophisticated
    prof_LS_ss = 0 # knowledge employed small customer base, sophisticated
    prof_BS_ss = 0 # knowledge employed big customer base, sophisticated
end

@with_kw struct Firm_sol
    #! Solution for firm parameters
    l_LU_ss = 0 # Labor good production small customer base, unsophisticated
    l_BU_ss = 0 # Labor good production big customer base, unsophisticated
    l_LS_ss = 0 # Labor good production small customer base, sophisticated
    l_BS_ss = 0 # Labor good production big customer base, sophisticated

    l_G_LU_ss = 0 # Labor good production small customer base, unsophisticated
    l_G_BU_ss = 0 # Labor good production big customer base, unsophisticated
    l_G_LS_ss = 0 # Labor good production small customer base, sophisticated
    l_G_BS_ss = 0 # Labor good production big customer base, sophisticated

    D_S_LU_ss = 0 # data sharing small customer base, unsophisticated
    D_S_BU_ss = 0 # data sharing big customer base, unsophisticated
    D_S_LS_ss = 0 # data sharing small customer base, sophisticated
    D_S_BS_ss = 0 # data sharing big customer base, sophisticated

    D_G_LU_ss = 0 # data generated small customer base, unsophisticated
    D_G_BU_ss = 0 # data generated big customer base, unsophisticated
    D_G_LS_ss = 0 # data generated small customer base, sophisticated
    D_G_BS_ss = 0 # data generated big customer base, sophisticated

    D_I_LU_ss = 0 # data generated small customer base, unsophisticated
    D_I_BU_ss = 0 # data generated big customer base, unsophisticated
    D_I_LS_ss = 0 # data generated small customer base, sophisticated
    D_I_BS_ss = 0 # data generated big customer base, sophisticated

    D_E_LU_ss = 0 # data generated small customer base, unsophisticated
    D_E_BU_ss = 0 # data generated big customer base, unsophisticated
    D_E_LS_ss = 0 # data generated small customer base, sophisticated
    D_E_BS_ss = 0 # data generated big customer base, sophisticated

    D_LU_ss = 0 # data generated small customer base, unsophisticated
    D_BU_ss = 0 # data generated big customer base, unsophisticated
    D_LS_ss = 0 # data generated small customer base, sophisticated
    D_BS_ss = 0 # data generated big customer base, sophisticated

    D_last_LU = 0 # data generated small customer base, unsophisticated
    D_last_BU = 0 # data generated big customer base, unsophisticated
    D_last_LS = 0 # data generated small customer base, sophisticated
    D_last_BS = 0 # data generated big customer base, sophisticated
    
    prof_LU_ss = 0 # knowledge employed small customer base, unsophisticated
    prof_BU_ss = 0 # knowledge employed big customer base, unsophisticated
    prof_LS_ss = 0 # knowledge employed small customer base, sophisticated
    prof_BS_ss = 0 # knowledge employed big customer base, sophisticated
end