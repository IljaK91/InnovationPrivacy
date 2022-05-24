"""
    Labor demand by A firms.
"""
function labor_demand_A(k_I, par)
    @unpack_Pars_Simple par
    ((γ * α_hat * v * C^α_C * k_I^δ_hat) / w)^(1 / (1 - α_hat))
end

"""
    Labor demand by B firms.
"""
function labor_demand_B(k_I, par)
    @unpack_Pars_Simple par
    (((1 - γ) * α_hat * v * C^α_C * k_I^δ_hat) / w)^(1 / (1 - α_hat))
end

"""
    How much capital A firms have in place given their production and selling decision.
"""
total_capital_A(k_I_b, k_I_P, par) = (1 - par.ϕ) * k_I_b + k_I_P

"""
    How much capital B firms have in place given their production and selling decision.
"""
total_capital_B(k_I_b, k_I_P) = k_I_b + k_I_P

"""
    Market clearing condition for traded capital. For any unit shipped, only 1-τ units arrive. Needs to be scaled by sector size (γ for A and 1-γ for B).
"""
market_clearing(k_I_b_B, k_I_b_A, par) = (1 - par.γ) * k_I_b_B + (1 - par.τ) * par.γ * k_I_b_A

"""
    Marginal cost of capital production by A firms.
"""
mc_A(k_I, par) = par.ca_A * par.cb_A * k_I^(par.cb_A - 1)

"""
    Marginal cost of capital production by B firms.
"""
mc_B(k_I, par) = par.ca_B * par.cb_B * k_I^(par.cb_B - 1)

"""
    Cost of capital production by A firms.
"""
c_A(k_I, par) = par.ca_A * k_I^par.cb_A
"""
    Cost of capital production by B firms.
"""
c_B(k_I, par) = par.ca_B * k_I^par.cb_B


"""
    How much capital is produced by A firms if they are indifferent between
    producing one more unit or selling one unit less.

        AND

    The firm does not sell all of its capital! Otherwise producing one more unit has a different payoff (sell and keep the rest instead of keeping it fully).
"""
k_I_prod_A(par) = ((par.p_I / (1 - par.ϕ)) / (par.ca_A * par.cb_A))^(1 / (par.cb_A - 1))

"""
    How much capital is produced by B firms if they are indifferent between
    producing one more unit or selling one unit less.

        AND

    The firm does not sell all of its capital! Otherwise producing one more unit has a different payoff (sell and keep the rest instead of keeping it fully).
"""
k_I_prod_B(par) = ((par.p_I / (1 - par.τ)) / (par.ca_B * par.cb_B))^(1 / (par.cb_B - 1))

"""
    The price for the B-good. C_A is the quantity of A-goods.
"""
function price_A(C_A, par::Pars_Simple)
    @unpack_Pars_Simple par
    γ * v * C^α_C * C_A^(-1 / σ)
end

"""
    The price for the B-good. C_B is the quantity of B-goods.
"""
function price_B(C_B, par::Pars_Simple)
    @unpack_Pars_Simple par
    (1 - γ) * v * C^α_C * C_B^(-1 / σ)
end

"""
    The wage rate that leads to market clearing (l_A + l_B = 1)
"""
function wage_rate(k_I_A, k_I_B; par)
    @unpack_Pars_Simple par
    c = (α_hat * v * C^α_C)
    c * ((γ * k_I_A^δ_hat)^(1 / (1 - α_hat)) + ((1 - γ) * k_I_B^δ_hat)^(1 / (1 - α_hat)))^(1 - α_hat)
end


"""
    Profit for both type of firms.
"""
function profit(par)
    @unpack_Pars_Simple par
    k_I_b_A, k_I_b_B = solve_firm_problems(par)

    if par.ϕ > 0 && loss_A(k_I_b_A, par) > 1e-8
        println("")
        k_I_A, l_A = capital_and_labor_con_A(k_I_b_A, par)
    else
        k_I_A, l_A = capital_and_labor_unc_A(k_I_b_A, par)
    end
    k_I_B, l_B = capital_and_labor_B(k_I_b_B, par)

    k_I_P_A = k_I_A - (1 - ϕ) * k_I_b_A # Sold capital is only partially lost.
    k_I_P_B = k_I_B - k_I_b_B

    profit_A = γ * v * C^α_C * l_A^α_hat * k_I_A^δ_hat - w * l_A - p_I * k_I_b_A - c_A(k_I_P_A, par)
    profit_B = (1 - γ) * v * C^α_C * l_B^α_hat * k_I_B^δ_hat - w * l_B - p_I / (1 - τ) * k_I_b_B - c_B(k_I_P_B, par)

    return profit_A, profit_B
end

# function labor_demand(k_I_j, par::Pars_Simple)
#     @unpack_Pars_Simple par
#     (α_hat * γ * C^α_C * (k_I_j)^δ_hat / w)^(1 / (1 - α_hat))
# end