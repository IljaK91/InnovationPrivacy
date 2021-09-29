@with_kw struct Pars
    # Economic Parameters
    θ       = 4
    α       = 0.5
    α_tilde = (θ - 1)/θ*α 
    δ       = 0.9
    # Shocks
    a_t = 0
    σ²ₐ = 1

    # Partial Equilibrium
    Y   = 1
    R   = 1.04
end

@with_kw struct BPars
    w_a # Bayesian weight on the private signal
    w_p # Bayesian weight on the prior
    σ²  # Posterior uncertainty
end