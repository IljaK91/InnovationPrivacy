"""
    Find the price of capital if the firm is not selling everything it produces.
"""
function find_price_info_unc(par)
    f(p_I) = begin
        @set! par.p_I = p_I^2
        k_I_b_A, k_I_b_B = solve_firm_problems_unc(par)
        diff = (1 - par.γ) * k_I_b_B + (1 - par.τ) * par.γ * k_I_b_A
        return diff # What one firm sells needs to be equal to what the other firm buys
    end
    find_zero(f, 0.25, tol = 1e-16)^2
end


function solve_firm_problems_unc(par)
    k_I_P_A = k_I_prod_A(par)
    if par.ϕ > 0 && loss_A(-k_I_P_A, par) > 0
        k_I_b_A = -k_I_P_A
    else
        f_A(x) = loss_A(x, par)
        k_I_b_A = find_zero(f_A, [-k_I_P_A + 0.00001, 5])
    end

    # B is always unconstrained, because they are buying capital
    k_I_P_B = k_I_prod_B(par)
    f_B(x) = loss_B(x, par)
    k_I_b_B = find_zero(f_B, [-k_I_P_B + 0.0001, 5])

    return k_I_b_A, k_I_b_B
end

"""
    This is much more robust!
"""
function solve_agg_C_iterative_unc(par; tol = 1e-10)
    diff = 10
    C_guess = par.C
    w_guess = par.w
    while diff > tol
        @set! par.C = C_guess
        #@show par.C
        C_sol, w_sol = find_agg_C_unc(par)
        @set! par.w = w_sol
        diff = copy((C_sol - C_guess)^2 + (w_sol - w_guess)^2)
        if diff > tol
            C_guess = copy(C_sol)
            w_guess = copy(w_sol)
        else
            @set! par.w = w_sol
            @set! par.C = C_guess
            @set! par.p_I = find_price_info_unc(par)
            @set! par.k_I_P_A_ss = k_I_prod_A(par)
            @set! par.k_I_P_B_ss = k_I_prod_B(par)

            f_A(x) = loss_A(x, par)
            @set! par.k_I_b_A_ss = find_zero(f_A, [-par.k_I_P_A_ss + 0.00001, 5])

            f_B(x) = loss_B(x, par)
            @set! par.k_I_b_B_ss = find_zero(f_B, [-par.k_I_P_B_ss + 0.00001, 5])

            @set! par.k_I_A_ss = total_capital_A(par.k_I_b_A_ss, par.k_I_P_A_ss, par)
            @set! par.k_I_B_ss = total_capital_B(par.k_I_b_B_ss, par.k_I_P_B_ss)

            @set! par.l_A_ss = labor_demand_A(par.k_I_A_ss, par)
            @set! par.l_B_ss = labor_demand_B(par.k_I_B_ss, par)
        end
    end
    return par
end

function find_agg_C_unc(par)
    @set! par.p_I = find_price_info_unc(par)
    k_I_b_A, k_I_b_B = solve_firm_problems_unc(par)

    k_I_A, l_A = capital_and_labor_unc_A(k_I_b_A, par)
    k_I_B, l_B = capital_and_labor_B(k_I_b_B, par)

    C_A = firm_prod(l_A, k_I_A; par)
    C_B = firm_prod(l_B, k_I_B; par)

    return agg_prod(C_A, C_B, par), wage_rate(k_I_A, k_I_B; par)
end

"""
    Profit for both type of firms.
"""
function profit_unc(par)
    @unpack_Pars_Simple par

    profit_A = γ * v * C^α_C * l_A_ss^α_hat * k_I_A_ss^δ_hat - w * l_A_ss - p_I * k_I_b_A_ss - c_A(k_I_P_A_ss, par)
    profit_B = (1 - γ) * v * C^α_C * l_B_ss^α_hat * k_I_B_ss^δ_hat - w * l_B_ss - p_I / (1 - τ) * k_I_b_B_ss - c_B(k_I_P_B_ss, par)

    return profit_A, profit_B
end


#! Now, go to the constrained problem, but both firms can sell and buy!
"""
    Solve the firms problem if 
"""
function solve_firm_problems_con(par)
    k_I_P_A = k_I_prod_A(par)
    if par.ϕ > 0 && loss_A(-k_I_P_A, par) > 0
        k_I_b_A = -k_I_P_A
    else
        f_A(x) = loss_A(x, par)
        k_I_b_A = find_zero(f_A, [-k_I_P_A + 0.00001, 5])
    end

    # B is always unconstrained, because they are buying capital
    k_I_P_B = k_I_prod_B(par)
    f_B(x) = loss_B(x, par)
    k_I_b_B = find_zero(f_B, [-k_I_P_B + 0.0001, 5])

    return k_I_b_A, k_I_b_B
end