"""
    Realized firm profits taking the interest rate R as given
"""
function firm_profit(a, K; par::Pars)
    @unpack_Pars par
    Y^(1/θ)*exp(a)^((θ - 1)/θ)*K^α_tilde - R*K
end

function Eπ(s, K; par::Pars, bpar::BPars)
    @unpack_Pars par
    ea = EA(s, par, bpar)
    Y^(1/θ)*ea*K^α_tilde - R*K
end

"""
    Constructor for the Bayesian Parameters (weights on private signal and prior)
"""
function BPars(τ, par::Pars)
    @unpack_Pars par
    τ   = τ
    σ²  = 1/(τ + σ²ₐ^-1) # posterior uncertainty
    w_a = τ*σ²           # weight on private signal  
    w_p = σ²ₐ^-1*σ²      # weight on prior

    BPars(τ   = τ,
          w_a = w_a,
          w_p = w_p,
          σ²  = σ²)
end

signal(a, eps, bpar::BPars) = a + eps/√bpar.τ
signal(a, eps, τ) = a + eps/√τ
Tau_of_K_old(K) = K

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

extractK(pos; K_grid) = K_grid[pos[2]]


function solve_model(par::Pars; V_guess = nothing, conv = 1e-8)
    @unpack_Pars par

    V_guess    === nothing ? V_guess = similar(K_grid) : nothing
    K_sol      = zeros(a_n^2,K_n)
    V_new      = zeros(K_n)
    V_next_big = repeat(V_guess, inner = a_n^2)
    nodes, w   = qnwnorm(a_n, [par.a_t, 0], [par.σ²ₐ, 1])
    a          = nodes[:,1]
    ε          = nodes[:,2]
    diff       = 1e8
    while abs(diff) > conv
        V_new   = iterate(V_guess, K_sol, V_new, a, ε, w, V_next_big, par)
        diff    = sum((V_guess - V_new).^2)
        @show diff
        V_guess = copy(V_new)
    end
    return V_new
end

function iterate(V_guess, K_sol, V_new, a, ε, w, V_next_big, par::Pars)
    @unpack_Pars par
    for i in 1:K_n
        K_old          = K_grid[i]
        τ              = Tau_of_K_old(K_old)
        s              = signal.(a, ε, τ)
        grid           = gridmake(s, K_grid)
        W_big          = Eπ.(grid[:,1], grid[:,2]; par, bpar = BPars(τ, par)) .+ (1-δ).*K_old - γ/2*grid[:,2].^2/K_old + R^-1*V_next_big
        W_big_reshaped = reshape(W_big, (a_n^2, K_n))
        W, pos         = findmax(W_big_reshaped, dims = 2)
        K_sol[:,i]     = extractK.(pos; K_grid)
        #V_new_big[:,i] = sum(w.*(firm_profit.(a, K_sol[:,i]; par)[:] + par.R^-1*extractK.(pos; K_grid = V_guess)[:]))
        V_new[i]       = sum(w.*(firm_profit.(a, K_sol[:,i]; par)[:] .+ (1-δ).*K_old - γ/2*K_sol[:,i].^2/K_old + R^-1*extractK.(pos; K_grid = V_guess)[:]))
    end
    return V_new
end


function K_solve(s, V, K_old; par::Pars)
    @unpack_Pars par
    W_big  = Eπ.(s, K_grid; par, bpar = BPars(Tau_of_K_old(K_old), par)) .+ (1-δ).*K_old - γ/2*K_grid.^2/K_old + R^-1*V
    W, pos = findmax(W_big)

    return K_grid[pos]
end