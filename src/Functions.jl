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
function BPar(β,par::Pars)
    @unpack_Pars par
    
    w_a = β/(β + σ²ₐ^-1)  
    w_p = σ²ₐ^-1/(β + σ²ₐ^-1)

    BPars(w_a = w_a,
          w_p = w_p)
end

function capital_demand(A)

end