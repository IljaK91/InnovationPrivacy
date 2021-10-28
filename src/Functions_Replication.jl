"""
    Constructor for the Bayesian Parameters (weights on private signal and prior).
        Tested
"""
function BPars(n, par::ParsGrowth)
    @assert n >= 0
    @unpack_ParsGrowth par
    V_n  = 1/(σ²_θ^-1 + n*σ²_ϵ^-1)
    BParsGrowth(n    = n,
                V_n  = V_n,
                w_pn = σ²_θ^-1*V_n,
                w_an = n*σ²_ϵ^-1*V_n)
end

BPars(n; par::ParsGrowth) = BParsGrowth(n, par)

"""
    The expected mean given n signals with mean a_bar.
"""
μ_n(a_bar, par::ParsGrowth, bpar::BParsGrowth) = bpar.w_pn*par.θ_bar + bpar.w_an*a_bar 
μ_n(a_bar; par::ParsGrowth, bpar::BParsGrowth) = μ_n(a_bar, par, bpar)
# Now lets write a function to compute the weights on each node

function Transition_Probability(par::ParsGrowth, bpar::BParsGrowth; normalize::Symbol = :yes)
    @unpack_ParsGrowth par
    @unpack_BParsGrowth bpar
    Π    = zeros(N, N)
    μ_ns = μ_n.(a_grid, par = par, bpar = bpar)

    # Rows: Today's state i
    # Columns: Tomorrow's state j
    for i in 1:N, j in 1:N
        if j == 1
            d      = (n+1)*(a_grid[j+1] - a_grid[j])
            Π[i,j] = normcdf(((n+1)*a_grid[j] - n*a_grid[i] - μ_ns[i] + d/2)/√(V_n + σ²_ϵ))
        elseif j == N
            d      = (n+1)*(a_grid[j] - a_grid[j-1])
            Π[i,j] = 1 - normcdf(((n+1)*a_grid[j] - n*a_grid[i] - μ_ns[i] - d/2)/√(V_n + σ²_ϵ))
        else
            d_up      = (n+1)*(a_grid[j+1] - a_grid[j])
            d_down    = (n+1)*(a_grid[j]   - a_grid[j-1])
            Π[i,j] = normcdf(((n+1)*a_grid[j] - n*a_grid[i] - μ_ns[i] + d_up/2)/√(V_n + σ²_ϵ)) -
                     normcdf(((n+1)*a_grid[j] - n*a_grid[i] - μ_ns[i] - d_down/2)/√(V_n + σ²_ϵ))
        end
    end
    # Renormalize the transition matrix
    if normalize == :yes
        for i in 1:N
            Π[i,:] = Π[i,:]./sum(Π[i,:])
        end
    end
    return Π
end

Transition_Probability(n::Integer, par::ParsGrowth; normalize::Symbol = :yes) = Transition_Probability(par::ParsGrowth, BPars(n, par); normalize)

function find_limit(par::ParsGrowth; tol = 1e-8)
    diff = 1e8
    n    = 5
    while diff > tol
        Π_first = Transition_Probability(n, par)
        Π_next  = Transition_Probability(n+1, par)
        
        diff = sum((Π_first - Π_next).^2)
        if diff > tol
            n += 1
        end
    end
    return n
end

"""
    Expectation of next exp(a/σ) having received n signals with mean a_bar.
        Seems correct.
"""
function EA(a_bar, par::ParsGrowth, bpar::BParsGrowth)
    @unpack_ParsGrowth par
    @unpack_BParsGrowth bpar
    exp(μ_n(a_bar, par, bpar)/σ + ϵ_bar/σ + 1/(2*σ^2)*(V_n + σ²_ϵ))
end

"""
    Expected Profits given the optimal labor demand decision.
        Tested
"""
function Eπ(a_bar, par::ParsGrowth, bpar::BParsGrowth)
    @unpack_ParsGrowth par
    @unpack_BParsGrowth bpar
    ea  = EA(a_bar, par, bpar)
    con = (σ-1)^(σ-1)/(σ^σ)
    u^(σ-1)*ea^σ*con - w*f
end

Eπ(a_bar; par::ParsGrowth, bpar::BParsGrowth) = Eπ(a_bar, par, bpar)

"""
    Expected Revenue given the optimal labor demand decision.
        Tested
"""
function Er(a_bar, par::ParsGrowth, bpar::BParsGrowth)
    @unpack_ParsGrowth par
    @unpack_BParsGrowth bpar
    ea  = EA(a_bar, par, bpar)
    con = ((σ-1)/σ)^(σ-1)
    u^(σ-1)*ea^σ*con
end

Er(a_bar; par::ParsGrowth, bpar::BParsGrowth) = Er(a_bar, par, bpar)

"""
    Labor costs resulting from the optimal labor demand decision.
"""
function labor_costs(a_bar, par::ParsGrowth, bpar::BParsGrowth)
    @unpack_ParsGrowth par
    @unpack_BParsGrowth bpar
    ea  = EA(a_bar, par, bpar)
    con = ((σ-1)/σ)^σ
    u^(σ-1)*ea^σ*con
end

labor_costs(a_bar; par::ParsGrowth, bpar::BParsGrowth) = labor_costs(a_bar, par, bpar)

"""
    For some n_max, compute the steady state value function.
"""
function find_V_n_max(n_max, par::ParsGrowth; bpar::BParsGrowth = BPars(n_max, par))
    @unpack_ParsGrowth par
    Π       = Transition_Probability(n_max, par)
    V_guess = zeros(N)
    V_new   = zeros(N)
    diff = 1e8
    while diff > 1e-8
        V_new = max.(Eπ.(a_grid; par, bpar) + β*(1-δ)*Π*V_guess,0)
        diff  = sum((V_new - V_guess).^2)
        diff  > 1e-8 ? V_guess = copy(V_new) : nothing
    end
    return V_new
end

function find_V_n(n, V_n_plus_1, par::ParsGrowth; bpar::BParsGrowth = BPars(n, par))
    @unpack_ParsGrowth par
    Π       = Transition_Probability(n, par)
    max.(Eπ.(a_grid, par = par, bpar = bpar) + β*(1-δ)*Π*V_n_plus_1,0)
end

"""
    Value function of an entrant with forced continuation. Necessary to find the cutoff for u.
"""
function find_V_0(V_1, par::ParsGrowth; bpar::BParsGrowth = BPars(0, par))
    @unpack_ParsGrowth par
    Π       = Transition_Probability(0, par)
    Eπ.(a_grid, par = par, bpar = bpar) + β*(1-δ)*Π*V_1
end

function find_all_n(n_max, V_n_max, par::ParsGrowth)
    V        = zeros(par.N, n_max+1)
    V[:,end] = V_n_max
    for n in n_max:-1:1
        V[:,n] = find_V_n(n-1, V[:,n+1], par)
    end
    return V
end

function find_u_star(par::ParsGrowth)
    n_max = find_limit(par, tol = 1e-16)
    f(u) = begin
        @set! par.u = u
        V_n_max = find_V_n_max(n_max, par)
        V       = find_all_n(n_max, V_n_max, par)
        return find_V_0(V[:,2], par)[1] # This is the zeroth value function with forced continuation.
    end
    find_zero(f, [1e-8, 2])
end

function u_to_z(par::ParsGrowth)
    @unpack_ParsGrowth par
    1/(σ-1)*log(u^(σ-1)*w^(σ-1)/(P^(σ-1)*Y))
end

function z_to_u(z, par::ParsGrowth)
    @unpack_ParsGrowth par
    (exp(z)^(σ-1)*P^(σ-1)*Y/(w^(σ-1)))^(1/(σ-1))
end

# function Firm_Distribution(par::ParsGrowth; normalize::Symbol = :yes)
#     @unpack_ParsGrowth par
#     Π    = zeros(N_z)
#     d    = z_grid[2] - z_grid[1] # equidistant grid 
#     # Rows: Today's state i
#     # Columns: Tomorrow's state j
#     for i in 1:N_z
#         if i == 1
#             Π[i] = normcdf((z_grid[i] - z_bar + d/2)/√(σ²_z))
#         elseif i == N
#             Π[i] = 1 - normcdf((z_grid[i] - z_bar - d/2)/√(σ²_z))
#         else
#             Π[i] = normcdf((z_grid[i] - z_bar + d/2)/√(σ²_z)) - normcdf((z_grid[i] - z_bar - d/2)/√(σ²_z))
#         end
#     end
#     # Renormalize the transition matrix
#     if normalize == :yes
#             Π = Π./sum(Π)
#     end
#     return Π
# end

function Firm_Distribution(par::ParsGrowth; normalize::Symbol = :yes)
    @unpack_ParsGrowth par
    Π    = zeros(N_z)
    d    = (u_grid[2:end] - u_grid[1:end-1])./2
    # Rows: Today's state i
    # Columns: Tomorrow's state j
    for i in 1:N_z
        if i == 1
            Π[i] = cdf(Pareto(ξ,u_star), u_grid[i] + d[i])
        elseif i == N_z
            Π[i] = 1 - cdf(Pareto(ξ,u_star), u_grid[i] - d[i-1])
        else
            Π[i] = cdf(Pareto(ξ,u_star), u_grid[i] + d[i]) - cdf(Pareto(ξ,u_star),u_grid[i] - d[i-1])
        end
    end
    # Renormalize the transition matrix
    if normalize == :yes
            Π = Π./sum(Π)
    end
    return Π
end

"""
    How the population of firms evolves from one period to the next.
"""
function iterate_population(n_max, Π_set, V_set, par::ParsGrowth)
    # Initialize array for the distribution of firms
    pop              = zeros(par.N, n_max)
    # Which firms survive / decide not to exit each period
    V_set[V_set.>0] .= 1 
    # Set the first column of the firm distribution. All of these firms survive and see their first sales and learn from them.
    # A fraction δ dies after production in their first period.
    pop[:,1]         = (1-par.δ).*Π_set[1][1,:] # Distribution over different a_bar after the initial period

    # Then start from the second period (n=2) on and go to the last period.
    for n in 2:1:n_max
        # A fraction δ dies.
        # Then take the population from last period and take away all firms that decide to exit due to too negative demand realizations.
        # Mismatch in timing: V_set starts at period 0, pop starts in period 1. Therefore n-1 for pop and n for V_set.
        # Then study which firms go to another state.
        # Finally, the transition matrix is appropriately transposed.
        pop[:,n] = (1-par.δ).*Π_set[n]'*(pop[:,n-1].*V_set[:,n])
    end
    return pop
end

function Equilibrium_Distribution(n_max, par::ParsGrowth)
    pop   = zeros(par.N_z, par.N, n_max)
    Π_u   = Firm_Distribution(par)
    Π_set = [Transition_Probability(par,BPars(n,par)) for n in 0:n_max] # find all Transition Matrices 
    for i in 1:length(par.u_grid)
        @set! par.u = par.u_grid[i]
        V_set      = find_all_n(n_max, find_V_n_max(n_max,par), par)# find all V for a given u   
        pop[i,:,:] = Π_u[i].*iterate_population(n_max, Π_set, V_set, par)
    end
    return pop
end

function Eπ_Distribution(n_max, par::ParsGrowth)
    Eπs = zeros(par.N_z, par.N, n_max)
    for i in 1:length(par.u_grid), j in 1:length(par.a_grid), n in 1:n_max
        @set! par.u = par.u_grid[i]
        bpar        = BPars(n, par)
        Eπs[i,j,n]  = Eπ(par.a_grid[j], par, bpar)    
    end
    return Eπs
end

function Er_Distribution(n_max, par::ParsGrowth)
    Ers = zeros(par.N_z, par.N, n_max)
    for i in 1:length(par.u_grid), j in 1:length(par.a_grid), n in 1:n_max
        @set! par.u = par.u_grid[i]
        bpar        = BPars(n, par)
        Ers[i,j,n]  = Er(par.a_grid[j], par, bpar)    
    end
    return Ers
end