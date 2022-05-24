#! In this folder I collect functions that are for the version of the model with depreciating information.

"""
    Calculate prior uncertainty after n periods of a stable data regulation regime.
"""
function CalcPriorUncertainty(n, past_posterior; par::ParsGrowth)
    if n == 0
        par.σ²_θ
    else
        par.ρ^2 * past_posterior + par.σ²_η
    end
end

function CalcPosteriorUncertainty(past_prior, par::ParsGrowth)
    #1/(past_prior^-1 + par.σ²_ϵ^-1 + par.σ²_η^-1)
    1/(past_prior^-1 + par.σ²_ϵ^-1)
end

function CalcAllUncertainty(n_max; par::ParsGrowth)
    prior_unc = zeros(n_max + 1)
    post_unc = zeros(n_max + 1)
    for n in 0:n_max
        if n == 0
            prior_unc[n+1] = CalcPriorUncertainty(n, 0; par)
        else
            prior_unc[n+1] = CalcPriorUncertainty(n, post_unc[n]; par)
        end
        post_unc[n+1] = CalcPosteriorUncertainty(prior_unc[n+1], par::ParsGrowth)
    end

    return prior_unc, post_unc
end

function Transition_Probability_dep(par::ParsGrowth, bpar::BParsGrowth; normalize::Symbol = :yes)
    @unpack_ParsGrowth par
    @unpack_BParsGrowth bpar
    Π    = zeros(N, N)

    # Rows: Today's state i
    # Columns: Tomorrow's state j
    for i = 1:N, j = 1:N
        if j == 1
            d = a_grid[j+1] - a_grid[j]

            c_plus = (σ²_ϵ * ((V_n^-1 + σ²_ϵ^-1) * (a_grid[j] + d / 2) / ρ - V_n^-1 * a_grid[i]) - a_grid[i]) / √(V_n + σ²_ϵ)

            Π[i, j] = normcdf(c_plus)
        elseif j == N
            d = a_grid[j] - a_grid[j-1]

            c_minus = (σ²_ϵ * ((V_n^-1 + σ²_ϵ^-1) * (a_grid[j] - d / 2) / ρ - V_n^-1 * a_grid[i]) - a_grid[i]) / √(V_n + σ²_ϵ)

            Π[i, j] = 1 - normcdf(c_minus)
        else
            d_up = a_grid[j+1] - a_grid[j]
            d_down = a_grid[j] - a_grid[j-1]
        
            c_plus = (σ²_ϵ * ((V_n^-1 + σ²_ϵ^-1) * (a_grid[j] + d_up / 2) / ρ - V_n^-1 * a_grid[i]) - a_grid[i]) / √(V_n + σ²_ϵ)
        
            c_minus = (σ²_ϵ * ((V_n^-1 + σ²_ϵ^-1) * (a_grid[j] - d_down / 2) / ρ - V_n^-1 * a_grid[i]) - a_grid[i]) / √(V_n + σ²_ϵ)
        
            Π[i, j] = normcdf(c_plus) - normcdf(c_minus)
        end
    end
    # Renormalize the transition matrix
    if normalize == :yes
        for i = 1:N
            Π[i, :] = Π[i, :] ./ sum(Π[i, :])
        end
    end
    return Π
end