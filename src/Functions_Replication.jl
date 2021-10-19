"""
    Constructor for the Bayesian Parameters (weights on private signal and prior)
"""
function BPars(n, par::ParsGrowth)
    @unpack_ParsGrowth par
    BParsGrowth(n = n,
                V_n  = 1/(σ²_θ^-1 + n*σ²_ϵ^-1),
                w_pn = σ²_θ^-1/(σ²_θ^-1 + n*σ²_ϵ^-1),
                w_an = n*σ²_ϵ^-1/(σ²_θ^-1 + n*σ²_ϵ^-1))
end

BPars(n; par::ParsGrowth) = BParsGrowth(n, par)

μ_n(a_bar, par::ParsGrowth, bpar::BParsGrowth) = bpar.w_pn*par.θ_bar + bpar.w_an*a_bar 
μ_n(a_bar; par::ParsGrowth, bpar::BParsGrowth) = μ_n(a_bar, par, bpar)
# Now lets write a function to compute the weights on each node

function Transition_Probability(par::ParsGrowth, bpar::BParsGrowth; normalize::Symbol = :yes)
    @unpack_ParsGrowth par
    @unpack_BParsGrowth bpar
    @assert n >= 0
    Π    = zeros(N, N)
    μ_ns = μ_n.(a_grid, par = par, bpar = bpar)
    d    = (n+1)*(a_grid[2] - a_grid[1])
    # Rows: Today's state i
    # Columns: Tomorrow's state j
    for i in 1:N, j in 1:N
        if j == 1
            Π[i,j] = normcdf(((n+1)*a_grid[j] - n*a_grid[i] - μ_ns[i] + d/2)/√(V_n + σ²_ϵ))
        elseif j == N
            Π[i,j] = 1 - normcdf(((n+1)*a_grid[j] - n*a_grid[i] - μ_ns[i] - d/2)/√(V_n + σ²_ϵ))
        else
            Π[i,j] = normcdf(((n+1)*a_grid[j] - n*a_grid[i] - μ_ns[i] + d/2)/√(V_n + σ²_ϵ)) -
                     normcdf(((n+1)*a_grid[j] - n*a_grid[i] - μ_ns[i] - d/2)/√(V_n + σ²_ϵ))
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
Transition_Probability(n::Integer, par::ParsGrowth; normalize::Symbol) = Transition_Probability(par::ParsGrowth, BPars(n, par); normalize)