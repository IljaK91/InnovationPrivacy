#! Mistake I was making before: Treating each firm as if it was some kind of monopolist of a given size.
#! However, firms have actually equal size!

#! Remove strategic complementarities

"""
    Labor demand by old and sophisticated (OS) firm.
"""
function labor_demand_S(k, par)
    @unpack_Pars_OLG par
    ((α_hat * Y^α_Y * k^δ_hat) / w)^(1 / (1 - α_hat))
end

"""
    Labor demand by old and unsophisticated (OU) firm.
"""
function labor_demand_U(k, par)
    @unpack_Pars_OLG par
    ((α_hat * Y^α_Y * k^δ_hat) / w)^(1 / (1 - α_hat))
end

"""
    Labor demand by young (Y) firm.
"""
function labor_demand_Y(k, par)
    @unpack_Pars_OLG par
    ((α_hat * Y^α_Y * k^δ_hat) / w)^(1 / (1 - α_hat))
end

"""
    How much capital A firms have in place given their production and selling decision.
"""
total_capital(k_b, k_P) = k_b + k_P


"""
    Market clearing condition for traded capital. For now, unsophisticated and young firms buy data from old sophisticated firms. For any unit OU firms ship, only 1-τ units arrive. Needs to be scaled by sector size (γ for A and 1-γ for B).
"""
function market_clearing(k_S, k_U, k_Y, par)
    @unpack_Pars_OLG par
    #! Check whether k_U is positive or negative
    if k_U > 0
        μ * k_Y + μ * (1 - γ) * k_U + (1 - τ) * μ * γ * k_S
    elseif k_U < 0
        μ * k_Y + (1-τ) * μ * (1 - γ) * k_U + (1 - τ) * μ * γ * k_S
    end
end

"""
    Marginal cost of capital production by sophisticated firms.
"""
mc_S(k_I, par) = par.ca_S * par.cb_S * k_I^(par.cb_S - 1)

"""
    Marginal cost of capital production by B firms.
"""
mc_U(k_I, par) = par.ca_U * par.cb_U * k_I^(par.cb_U - 1)

"""
    Yost of capital production by A firms.
"""
c_S(k_I, par) = par.ca_S * k_I^par.cb_S
"""
    Yost of capital production by B firms.
"""
c_U(k_I, par) = par.ca_U * k_I^par.cb_U


"""
    How much capital is produced by A firms if they are indifferent between
    producing one more unit or selling one unit less.

        AND

    The firm does not sell all of its capital! Otherwise producing one more unit has a different payoff (sell and keep the rest instead of keeping it fully).
"""
k_prod_S(par) = (par.p_I / (par.ca_S * par.cb_S))^(1 / (par.cb_S - 1))

"""
    How much capital is produced by B firms if they are indifferent between
    producing one more unit or selling one unit less.

        AND

    The firm does not sell all of its capital! Otherwise producing one more unit has a different payoff (sell and keep the rest instead of keeping it fully).
"""
k_prod_U_buy(par) = ((par.p_I / (1 - par.τ)) / (par.ca_U * par.cb_U))^(1 / (par.cb_U - 1))
k_prod_U_sell(par) = (par.p_I / (par.ca_U * par.cb_U))^(1 / (par.cb_U - 1))

"""
    The price for an intermediate good is downwards sloping in the quantity produced.
"""
function price(Y_i, par::Pars_OLG)
    @unpack_Pars_OLG par
    Y^α_Y / (Y_i)^(1 / σ)
end

# """
#     The price for the B-good. Y_U is the quantity of B-goods.
# """
# function price_U(Y_U, par::Pars_OLG)
#     @unpack_Pars_OLG par
#     (Y / Y_U)^(1 / σ)
# end

"""
    The wage rate that leads to market clearing (l_S + l_U + l_Y = 1). Here, labor demand by each firm needs to be multiplied by the mass of the group of firms.
"""
function wage_rate(k_S, k_U, k_Y; par)
    @unpack_Pars_OLG par
    a = (1 / (1 - α_hat))
    α_hat * Y^α_Y * μ * ((γ * k_S^δ_hat)^a + ((1 - γ) * k_U^δ_hat)^a + (k_Y^δ_hat)^a)^(a^-1)
end


"""
    Profit for both type of firms.
"""
function profit(par)
    @unpack_Pars_OLG par
    k_T_S, k_T_U, k_T_Y = solve_firm_problems(par)

    k_S, l_S = capital_and_labor_S(k_T_S, par)
    k_U, l_U = capital_and_labor_U(k_T_U, par)
    k_Y, l_Y = capital_and_labor_Y(k_T_Y, par)

    k_P_S = k_S - k_T_S
    k_P_U = k_U - k_T_U
    
    @assert k_T_Y == k_Y

    #@assert k_T_U > 0
    
    k_P_Y = 0

    profit_S = Y^α_Y * l_S^α_hat * k_S^δ_hat - w * l_S - p_I * k_T_S - c_S(k_P_S, par)
    profit_U = Y^α_Y * l_U^α_hat * k_U^δ_hat - w * l_U - p_I / (1 - τ) * k_T_U - c_U(k_P_U, par)
    profit_Y = Y^α_Y * l_Y^α_hat * k_Y^δ_hat - w * l_Y - p_I / (1 - τ) * k_T_Y

    #@show k_T_S, k_T_U, k_T_Y, k_P_S, k_P_U
    return profit_S, profit_U, profit_Y
end

"""
    Solve the capital problem of each firm
"""
function solve_firm_problems(par)
    k_P_S = k_prod_S(par)
    min_val_S = -k_P_S + 0.00001
    f_S(x) = loss_S(x^2 + min_val_S, par)
    k_T_S = find_zero(f_S, 0.25)^2 + min_val_S

    # This part needs adjusting. I may be that the unsophisticated firm is in different regions!
    k_T_U = solve_U_problem(par)

    f_Y(x) = loss_Y(x^2, par)
    k_T_Y = find_zero(f_Y, 0.25)^2

    return k_T_S, k_T_U, k_T_Y
end

function solve_U_problem(par)
    k_P_U_buy = k_prod_U_buy(par)
    # Here I can improve and avoid setting bounds
    min_val_U_buy = -k_P_U_buy + 0.00001
    f_U_buy(x) = loss_U_buy(x^2 + min_val_U_buy, par)
    k_T_U_buy = find_zero(f_U_buy, 0.25)^2 + min_val_U_buy

    k_P_U_sell = k_prod_U_sell(par)
    min_val_U_sell = -k_P_U_sell + 0.00001
    f_U_sell(x) = loss_U_sell(x^2 + min_val_U_sell, par)
    k_T_U_sell = find_zero(f_U_sell, 0.25)^2 + min_val_U_sell

    #@show k_T_U_buy
    #@show k_T_U_sell

    #! Do some disambiguation
    if k_T_U_buy > 0 # First case: U firms actually buy at the buying price p_I / ( 1 - τ ).
        return k_T_U_buy
    elseif k_T_U_sell < 0 # Second case: U firms actually sell at the selling price p_I.
        return k_T_U_sell
    elseif k_T_U_buy < 0 && k_T_U_sell > 0 # Third case: Unsophisticated firms find it best to self-produce.
        return 0 # Unsophisticated firms do not trade
    end
end

function loss_S(k_T, par)
    @unpack_Pars_OLG par
    k, l = capital_and_labor_S(k_T, par)
    return FOC_K_S(k, l, par)
end

function loss_U_buy(k_T, par)
    @unpack_Pars_OLG par
    k, l = capital_and_labor_U_buy(k_T, par)
    return FOC_K_U_buy(k, l, par)
end

function loss_U_sell(k_T, par)
    @unpack_Pars_OLG par
    k, l = capital_and_labor_U_sell(k_T, par)
    return FOC_K_U_sell(k, l, par)
end

function loss_Y(k_T, par)
    @unpack_Pars_OLG par
    k, l = capital_and_labor_Y(k_T, par)
    return FOC_K_Y(k, l, par)
end

"""
    First order condition for capital employment sophisticated firms.

    Selling one unit less costs p_I.
"""
#! Why not solve analytically
function FOC_K_S(k, l, par)
    @unpack_Pars_OLG par
    p_I - δ_hat * Y^α_Y * l^α_hat * k^(δ_hat - 1)
end

"""
    First order condition for capital employment unsophisticated firms.

    Need to do some disambiguation in the case with three different types of firms:

    1. The unsophisticated firms buys capital from sophisticated firms.
    2. The unsophisticated firm completely self-produces data.
    3. The unsophisticated firm sells capital to young firms.

    Buying one unit more costs p_I / (1-τ).
"""
function FOC_K_U_buy(k, l, par)
    @unpack_Pars_OLG par
    p_I / (1 - τ) - δ_hat * Y^α_Y * l^α_hat * k^(δ_hat - 1)
end

function FOC_K_U_sell(k, l, par)
    @unpack_Pars_OLG par
    p_I - δ_hat * Y^α_Y * l^α_hat * k^(δ_hat - 1)
end

function FOC_K_U_prod(k, l, par)
    @unpack_Pars_OLG par
    ca_U*cb_U*k^(cb_U-1) - δ_hat * Y^α_Y * l^α_hat * k^(δ_hat - 1)
end

"""
    First order condition for capital employment young firms.

    Buying one unit more costs p_I / (1-τ).
"""
function FOC_K_Y(k, l, par)
    @unpack_Pars_OLG par
    p_I / (1 - τ) - δ_hat * Y^α_Y * l^α_hat * k^(δ_hat - 1)
end

"""
    Given how much capital the sophisticated firm buys, how much capital and labor does it employ.
"""
function capital_and_labor_S(k_T, par::Pars_OLG)
    k = total_capital(k_T, k_prod_S(par))
    l = labor_demand_S(k, par)
    return k, l
end

function capital_and_labor_U(k_T, par::Pars_OLG)
    if k_T > 0
        capital_and_labor_U_buy(k_T, par)
    elseif k_T < 0
        capital_and_labor_U_sell(k_T, par)
    elseif k_T == 0
        capital_and_labor_U_prod(par::Pars_OLG)
    end
end

function capital_and_labor_U_buy(k_T, par::Pars_OLG)
    k = total_capital(k_T, k_prod_U_buy(par))
    l = labor_demand_U(k, par)
    return k, l
end

function capital_and_labor_U_sell(k_T, par::Pars_OLG)
    k = total_capital(k_T, k_prod_U_sell(par))
    l = labor_demand_U(k, par)
    return k, l
end

function capital_and_labor_U_prod(par::Pars_OLG)
    f(k) = begin
        l = labor_demand_U(k, par)
        FOC_K_U_sell(k,l,par)
    end

    k = find_zero(f, [0.00001, 20])

    return k, labor_demand_U(k, par)
end

function capital_and_labor_Y(k_T, par::Pars_OLG)
    k = total_capital(k_T, 0)
    l = labor_demand_Y(k, par)
    return k, l
end

# capital_and_labor_U(k_T, par) = total_capital(k_T, k_prod_U(par)), labor_demand_U(k, par)
# capital_and_labor_Y(k_T, par) = total_capital(k_T, 0), labor_demand_Y(k, par)

function find_price_info(par)
    f(p_I) = begin
        @set! par.p_I = p_I^2
        k_T_S, k_T_U, k_T_Y = solve_firm_problems(par)
        return market_clearing(k_T_S, k_T_U, k_T_Y, par) # What one firm sells needs to be equal to what the other firm buys
    end
    find_zero(f, 0.25, tol = 1e-16)^2
end

function solve_agg_Y_iterative_unc(par; tol = 1e-10)
    diff = 10
    Y_guess = par.Y
    w_guess = par.w
    while diff > tol
        @set! par.Y = Y_guess
        #@show par.C
        Y_sol, w_sol = find_agg_Y(par)
        @set! par.w = w_sol
        diff = copy((Y_sol - Y_guess)^2 + (w_sol - w_guess)^2)
        if diff > tol
            Y_guess = copy(Y_sol)
            w_guess = copy(w_sol)
        else
        @set! par.w = w_sol
        @set! par.Y = Y_guess
        @set! par.p_I = find_price_info(par)
        @set! par.k_P_S_ss = k_prod_S(par)

    #@set! par.k_P_U_ss = k_prod_U(par)

        f_S(x) = loss_S(x, par)
        @set! par.k_T_S_ss = find_zero(f_S, [-par.k_P_S_ss + 0.00001, 5])

    #f_U(x) = loss_U(x, par)
    #@set! par.k_T_U_ss = find_zero(f_U, [-par.k_P_U_ss + 0.00001, 5])

        @set! par.k_T_U_ss = solve_U_problem(par)
        k_P_U_ss, l_U_ss = capital_and_labor_U(par.k_T_U_ss, par)
        @set! par.k_P_U_ss = k_P_U_ss

        f_Y(x) = loss_Y(x, par)
        @set! par.k_T_Y_ss = find_zero(f_Y, [0.00001, 5])

        @set! par.k_S_ss = total_capital(par.k_T_S_ss, par.k_P_S_ss)
        @set! par.k_U_ss = total_capital(par.k_T_U_ss, par.k_P_U_ss)
        @set! par.k_Y_ss = total_capital(par.k_T_Y_ss, 0)

        @set! par.l_S_ss = labor_demand_S(par.k_S_ss, par)
        @set! par.l_U_ss = labor_demand_U(par.k_U_ss, par)
        @set! par.l_Y_ss = labor_demand_Y(par.k_Y_ss, par)
        end     
    end
    return par
end

firm_prod(l, k; par) = l^par.α_hat * k^par.δ_hat

function find_agg_Y(par)
    @set! par.p_I = find_price_info(par)
    k_T_S, k_T_U, k_T_Y = solve_firm_problems(par)

    k_S, l_S = capital_and_labor_S(k_T_S, par)
    k_U, l_U = capital_and_labor_U(k_T_U, par)
    k_Y, l_Y = capital_and_labor_Y(k_T_Y, par)

    Y_S = firm_prod(l_S, k_S; par)
    Y_U = firm_prod(l_U, k_U; par)
    Y_Y = firm_prod(l_S, k_S; par)

    return agg_prod(Y_S, Y_U, Y_Y, par), wage_rate(k_S, k_U, k_Y; par)
end

function agg_prod(Y_S, Y_U, Y_Y, par)
    @unpack_Pars_OLG par
    a = ((σ - 1) / σ)
    if α_Y > 0 # strategic complementarities
        (μ * Y_Y^a + μ * γ * Y_S^a + μ * (1 - γ) * Y_U^a)^(a^-1)
    elseif α_Y == 0 # no strategic complementarities
        μ * Y_Y^a + μ * γ * Y_S^a + μ * (1 - γ) * Y_U^a
    else
        error("α_Y < 0 was not planned for right now")
    end
end

"""
    Find the mass of entrants
"""
function find_mass(par; tol = 1e-16)

    f(x) = begin
        @set! par.μ = x^2
        par = solve_agg_Y_iterative_unc(par)
        profit_S, profit_U, profit_Y = profit(par)
        exp_profit = profit_Y + par.β * (par.γ * profit_S + (1 - par.γ) * profit_U)
        return (exp_profit - par.f)^2
    end

    find_zero(f, 1.5)^2
end