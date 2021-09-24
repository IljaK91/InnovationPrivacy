@with_kw struct Pars begin
    # Economic Parameters
    θ       = 4
    α       = 0.5
    α_tilde = (θ - 1)/θ*α 

    # Shocks
    a_t = 0
    σ²ₐ = 1

    # Partial Equilibrium
    Y   = 1
    R   = 1
end

@with_kw struct BPars begin
    w_a # Bayesian weight on the private signal
    w_p # Bayesian weight on the prior
end