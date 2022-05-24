


function loss_A(x, par)
    @unpack_Pars_Simple par
    k_I_b = x
    k_I, l = capital_and_labor_unc_A(k_I_b, par)
    # k_I = total_capital_A(k_I_b, k_I_P, par)
    # l = labor_demand_A(k_I, par)
    # Residual
    return loss_A_only(k_I, l, par)
end

"""
    A is selling, such that selling an additional unit yields p_I, but costs only
    (1-ϕ) units of the marginal product of capital, as ϕ units are kept.
"""
function loss_A_only(k_I, l, par)
    @unpack_Pars_Simple par
    p_I - (1 - ϕ) * γ * δ_hat * v * C^α_C * l^α_hat * k_I^(δ_hat - 1)
end

"""
    B is buying, so buying costs p_I/(1-τ) and yields the full marginal product.
"""
function loss_B_only(k_I, l, par)
    @unpack_Pars_Simple par
    p_I / (1 - τ) - (1 - γ) * δ_hat * v * C^α_C * l^α_hat * k_I^(δ_hat - 1)
end

function loss_A_con(k_I, par)
    @unpack_Pars_Simple par
    l = labor_demand_A(k_I, par)

    return mc_A(k_I, par) - p_I - ϕ^δ_hat * γ * δ_hat * v * C^α_C * l^α_hat * k_I^(δ_hat - 1)
end

function loss_B(x, par)
    @unpack_Pars_Simple par
    k_I_b = x
    # k_I = total_capital_B(k_I_b, k_I_P, par)
    # l = labor_demand_B(k_I, par)
    k_I, l = capital_and_labor_B(k_I_b, par)
    # Residual
    # very important: B firms face a higher marginal cost of investment!
    loss_B_only(k_I, l, par)
end



"""
    Find the equilibrium for a given C using Optim.jl
    Advantage: For some p_I, there might be no solution to the firm's problem.
    However, it might still be that there is a solution for a different p_I
    without resorting to the constrained solution. Only if this solution is
    unavailable, should I try the constrained solution.
"""
function find_price_info_optim(par)

    
    lower = [1.25, -2.1]
    upper = [Inf, Inf]
    initial_x = [2.0, 2.0]
    f(x) = residuals(x, par)
    inner_optimizer = GradientDescent()
    results = optimize(f, g!, lower, upper, initial_x, Fminbox(inner_optimizer))

    #! Problem, bounds depend on one parameter
end

function residuals( x, par)
        @unpack_Pars_Simple par
        p_I = x[1]
        @set! par.p_I = p_I
        k_I_b_A = x[2]
        k_I_b_B = x[3]
        # How much capital is produced
        k_I_P_A = k_I_prod_A(par)
        #k_I_P_B = k_I_prod_B(par)
    
        k_I_A, l_A = capital_and_labor_unc_A(k_I_b_A, par)
        k_I_B, l_B = capital_and_labor_unc_B(k_I_b_B, par)
    
        res_1 = loss_A_only(k_I_A, l_A, par)
        res_2 = loss_B_only(k_I_B, l_B, par)
        res_3 = market_clearing(k_I_b_B, k_I_b_A, par) # market clearing!
        res_4 = capicity_constraint(k_I_b_A, k_I_P_A) # This is the constraint that A cannot sell more than it produces!
    
        return res_1^2 + res_2^2 + res_3^2 + res_4^2
end

#! Try the following: Write your own solution algorithm.
# 1. Guess p_I
# 2. Check whether a sign change occurs
# 2.1. If yes, good, find solution.
# 2.2., If no, try lower p_I.
# 3. If there is no sign change even with some minimal p_I, switch to the constraint problem.

#! Can the constraint problem ever be a solution?





function capicity_constraint(k_I_b_A, k_I_P_A)
    # note that k_I_b_A is negative
    if k_I_P_A + k_I_b_A < 0
        1000*(k_I_P_A + k_I_b_A)^2
    else
        0
    end
end




function solve_agg_C(par)
    f(C) = begin
        @set! par.C = C^2
        return C - find_agg_C(par)
    end
    find_zero(f, 0.05)^2
end

function find_agg_C(par)
    @set! par.p_I = find_price_info(par)
    k_I_b_A, k_I_b_B = solve_firm_problems(par)

    if par.ϕ > 0 && loss_A(k_I_b_A, par) > 1e-8
        k_I_A, l_A = capital_and_labor_con_A(k_I_b_A, par)
    else
        k_I_A, l_A = capital_and_labor_unc_A(k_I_b_A, par)
    end
    k_I_B, l_B = capital_and_labor_B(k_I_b_B, par)

    C_A = firm_prod(l_A, k_I_A; par)
    C_B = firm_prod(l_B, k_I_B; par)

    return agg_prod(C_A, C_B, par), wage_rate(k_I_A, k_I_B; par)
end


function find_wage_rate(par)
    k_I_A, k_I_B = capital_prod(par)
    wage_rate(k_I_A, k_I_B; par)
end
function agg_prod(C_A, C_B, par)
    @unpack_Pars_Simple par
    (γ * C_A^((σ - 1) / σ) + (1 - γ) * C_B^((σ - 1) / σ))^(v * σ / (σ - 1))
end
function firm_prod(l, k_I; par)
    @unpack_Pars_Simple par
    l^α_hat * k_I^δ_hat
end

function solve_firm_problems(par)
    k_I_P_A = k_I_prod_A(par)
    if loss_A(-k_I_P_A, par) ≈ 0
        k_I_b_A = -k_I_P_A
    elseif par.ϕ > 0 && loss_A(-k_I_P_A, par) > 0
        #println("Constrained Optimization")
        # constrained, sell everything that is produced
        f_A_con(x) = loss_A_con(x^2, par)
        k_I_A = find_zero(f_A_con, 1)^2 # this is the production decision - always positive!
        k_I_P_A = k_I_A / (1 - par.ϕ) # more is produced than employed
        k_I_b_A = -k_I_P_A # production equals sharing
    else # unconstrained, sell less than what is produced
        k_I_P_A = k_I_prod_A(par)
        f_A(x) = loss_A(x, par)
        @show f_A(-k_I_P_A)
        @show f_A(100)
        @show par.p_I
        k_I_b_A = find_zero(f_A, [-k_I_P_A, 1000])
    end

    # B is always unconstrained, because they are buying capital
    k_I_P_B = k_I_prod_B(par)
    f_B(x) = loss_B(x, par)
    k_I_b_B = find_zero(f_B, [-k_I_P_B + 0.0001, 20])

    return k_I_b_A, k_I_b_B
end

function capital_prod(par)
    k_I_b_A, k_I_b_B = solve_firm_problems(par)
    if par.ϕ > 0 && loss_A(k_I_b_A, par) > 1e-8
        println("Constrained Optimization")
        k_I_A, l_A = capital_and_labor_con_A(k_I_b_A, par)
    else
        println("Unconstrained Optimization")
        k_I_A, l_A = capital_and_labor_unc_A(k_I_b_A, par)
    end
    k_I_B, l_B = capital_and_labor_B(k_I_b_B, par)

    return k_I_A, k_I_B
end


"""
    This is much more robust!
"""
function solve_agg_C_iterative(par; tol = 1e-10)
    diff = 10
    C_guess = par.C
    w_guess = par.w
    while diff > tol
        @set! par.C = C_guess
        #@show par.C
        C_sol, w_sol = find_agg_C(par)
        @set! par.w = w_sol
        diff = copy((C_sol - C_guess)^2 + (w_sol - w_guess)^2)
        if diff > tol
            C_guess = copy(C_sol)
            w_guess = copy(w_sol)
        else
            @set! par.w = w_sol
            @set! par.C = C_guess
            @set! par.p_I = find_price_info(par)
            return par
        end
    end
end

capital_of_tau(τ) = capital_prod(solve_agg_C_iterative(Pars_Simple(τ = τ)))



function capital_and_labor_unc_A(k_I_b_A, par)
    #! Cannot put this here, as I otherwise call loss_A recursively
    # if par.ϕ > 0 && loss_A(k_I_b_A, par) > 0
    #     # constrained, sell everything that is produced
    #     f_A_con(x) = loss_A_con(x^2, par)
    #     k_I_A = find_zero(f_A_con, 1)^2 # this is the production decision - always positive!
    # else
    # This is for an unconstrained firm
    k_I_P_A = k_I_prod_A(par)
    k_I_A = total_capital_A(k_I_b_A, k_I_P_A, par)
    # end
    l_A = labor_demand_A(k_I_A, par)
    return k_I_A, l_A
end

function capital_and_labor_con_A(k_I_b_A, par)
    #! Cannot put this here, as I otherwise call loss_A recursively
    @assert par.ϕ > 0
    # if par.ϕ > 0 && loss_A(k_I_b_A, par) > 0
    #     # constrained, sell everything that is produced
    f_A_con(x) = loss_A_con(x^2, par)
    k_I_A = find_zero(f_A_con, 1)^2 # this is the production decision - always positive!
    # else
    # This is for an unconstrained firm
    # k_I_P_A = k_I_prod_A(par)
    # k_I_A = total_capital_A(k_I_b_A, k_I_P_A, par)
    # # end
    l_A = labor_demand_A(k_I_A, par)
    return k_I_A, l_A
end

function capital_and_labor_B(k_I_b_B, par)
    k_I_P_B = k_I_prod_B(par)
    k_I_B = total_capital_B(k_I_b_B, k_I_P_B)
    l_B = labor_demand_B(k_I_B, par)
    return k_I_B, l_B
end

# function profit_no_effort(par)
#     @unpack_Pars_Simple par
#     sol_A, sol_B = solve_firm_problems(par)

#     l_B = sol_B[1]
#     k_I_b_B = sol_B[2]
#     k_I_P_B = k_I_prod_B(par)
#     k_I_B = k_I_b_B + k_I_P_B

#     l_A = sol_A[1]
#     k_I_b_A = sol_A[2]
#     k_I_P_A = k_I_prod_A(par)
#     k_I_A = k_I_b_A + k_I_P_A

#     profit_A = γ * v * C^α_C * l_A^α_hat * k_I_A^δ_hat - w * l_A - p_I * k_I_b_A - c_A(k_I_P_A, par)
#     profit_B = (1 - γ) * v * C^α_C * l_B^α_hat * k_I_B^δ_hat - w * l_B - p_I * k_I_b_B/(1-τ) - c_B(k_I_P_B, par)

#     return profit_A, profit_B
# end