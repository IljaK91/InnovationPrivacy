"""
    Realized firm profits taking the interest rate R as given
"""
function firm_profit(a, K, par::Pars)
    @unpack_Pars par
    exp(a^((θ - 1)/θ))*K^α_tilde - R*K
end

"""
    Constructor for the Bayesian Parameters (weights on private signal and prior)
"""
function BPars(β, par::Pars)
    @unpack_Pars par
    
    σ²  = 1/(β + σ²ₐ^-1) # posterior uncertainty
    w_a = β*σ² # weight on private signal  
    w_p = σ²ₐ^-1*σ² # weight on prior

    BPars(w_a = w_a,
          w_p = w_p,
          σ²  = σ²)
end

function EA(a, par::Pars, bpar::BPars)
    @unpack_Pars par
    @unpack_BPars bpar
    exp(((θ - 1)/θ)*w_a*a + ((θ - 1)/θ)^2*σ²/2)
end

function capital_demand(a, par::Pars, bpar::BPars)
    @unpack_Pars par
    @unpack_BPars bpar
    (Y^(1/θ)*EA(a, par, bpar)/R)^(1/(1-α_tilde))
end