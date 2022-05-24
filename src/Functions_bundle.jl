"""
    Find the prices of capital if the firm is not selling everything it produces.
"""
function find_price_info_unc(par)

    lower = [0.0 + eps(), 0.0 + eps()]
    upper = [5, 5]

    function loss(xy)
        sum((f(xy...) .- c) .^ 2)
    end

    opt = optimize(loss, lower, upper, [0.5, 0.5], Fminbox(LBFGS()))
    f(p_I) = begin
        @set! par.p_I = p_I^2
        k_I_b_A, k_I_b_B = solve_firm_problems_unc(par)
        diff = (1 - par.γ) * k_I_b_B + (1 - par.τ) * par.γ * k_I_b_A
        return diff # What one firm sells needs to be equal to what the other firm buys
    end
    find_zero(f, 0.25, tol = 1e-16)^2
end

"""
    Find capital allocation, starting with guesses for the prices of different kinds of capital.
"""
function solve_firm_problems_bundle(par)
    k_I_P_A = k_I_prod_A_bundle(par)
    k_I_P_B = k_I_prod_B_bundle(par)

    K_A_guess = 2 * k_I_P_A
    K_B_guess = 2 * k_I_P_B

    #! Now I can infer the other variables
    
    f_A(x) = loss_A(x, par)
    k_I_b_A = find_zero(f_A, [-k_I_P_A + 0.00001, 5])

    # B is always unconstrained, because they are buying capital
    k_I_P_B = k_I_prod_B(par)
    f_B(x) = loss_B(x, par)
    k_I_b_B = find_zero(f_B, [-k_I_P_B + 0.0001, 5])

    return k_I_b_A, k_I_b_B
end

"""
   How much capital A firms produce at an interior solution. 
"""
k_I_prod_A_bundle(par) = ((par.p_I_A / (1 - par.ϕ)) / (par.ca_A * par.cb_A))^(1 / (par.cb_A - 1))

"""
   How much capital A firms produce at an interior solution. 
"""
k_I_prod_A_bundle(par) = ((par.p_I_B / (1 - par.ϕ)) / (par.ca_B * par.cb_B))^(1 / (par.cb_B - 1))

"""
    How much capital is produced by B firms if they are indifferent between
    producing one more unit or selling one unit less.

        AND

    The firm does not sell all of its capital! Otherwise producing one more unit has a different payoff (sell and keep the rest instead of keeping it fully).
"""
k_I_prod_B(par) = ((par.p_I / (1 - par.τ)) / (par.ca_B * par.cb_B))^(1 / (par.cb_B - 1))