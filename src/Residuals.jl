"""
    This might not be the right way to go about it to solve the problem for small firms.
"""
function residuals_unc(x; par::Pars_v3, type::Symbol)
    l_P = x[1]^2
    l_G = x[2]^2
    D_E = x[3]^2
    D_I = x[4]^2

    D = bundle(D_I, D_E; par, type)
    l = labor_demand(l_P, D; par, type)

    #! FOCs
    res = zeros(4)
    res[1] = FOC_data_buying(l_P, D, l, D_E; par, type)
    res[2] = FOC_data_gen(l_G, l_P, D, l, D_I; par, type)
    res[3] = FOC_data_proc(l_P, D, l; par, type)
    res[4] = FOC_data_sharing(l_P, D, D_I, l; par, type)

    return sum(abs.(res)) # Try to force the algorithm to come up with something more precise
end

function residuals_unc_nl(x, res; par::Pars_v3, type::Symbol)
    l_G = x[1]^2
    D_E = x[2]^2
    D_I = x[3]^2

    D = bundle(D_I, D_E; par, type)
    l = labor_demand(D; par, type)

    #! FOCs
    res[1] = FOC_data_buying(D, l, D_E; par, type)
    res[2] = FOC_data_gen(l_G, D, l, D_I; par, type)
    #res[3] = FOC_data_proc(l_P, D, l; par, type, sec_period)
    res[3] = FOC_data_sharing(D, D_I, l; par, type)

    return res # Try to force the algorithm to come up with something more precise
end

function residuals_unc_nl2(x, res; par::Pars_v3, type::Symbol, sec_period::Symbol = :no)
    #! FOCs
    if sec_period == :yes
        l_P = x[1]^2
        l_G = x[2]^2
        D_E = x[3]^2
        D_I = x[4]^2
    
        D = bundle(D_I, D_E; par, type)
        l = labor_demand(l_P, D; par, type, sec_period)
        res[1] = FOC_data_buying(l_P, D, l, D_E; par, type, sec_period)
        res[2] = FOC_data_gen(l_G, l_P, D, l, D_I; par, type, sec_period)
        res[3] = FOC_data_proc(l_P, D, l; par, type, sec_period)
        res[4] = FOC_data_sharing(l_P, D, D_I, l; par, type, sec_period)
    elseif sec_period == :no
        l_P = x[1]^2
        l_G = x[2]^2
        D_E = x[3]^2
        D_I = x[4]^2
    
        l_P2 = x[5]^2
        l_G2 = x[6]^2
        D_E2 = x[7]^2
        D_I2 = x[8]^2
        D = bundle(D_I, D_E; par, type)
        D2 = bundle(D_I2, D_E2; par, type)
        l = labor_demand(l_P, D; par, type)
        l2 = labor_demand(l_P2, D2; par, type)
        #! This is the solution to the decision in the second period
        res[1] = FOC_data_buying(l_P2, D2, l2, D_E2; par, type, sec_period=:yes)
        res[2] = FOC_data_gen(l_G2, l_P2, D2, l2, D_I2; par, type, sec_period=:yes)
        res[3] = FOC_data_proc(l_P2, D2, l2; par, type, sec_period=:yes)
        res[4] = FOC_data_sharing(l_P2, D2, D_I2, l2; par, type, sec_period=:yes)
    
        #! This is the solution to the decision in the first period
        res[5] = FOC_data_buying2(l_P, D, l, D_E, l_P2, D2, l2; par, type)
        res[6] = FOC_data_gen2(l_G, l_P, D, l, D_I, l_P2, D2, l2; par, type)
        res[7] = FOC_data_proc2(l_P, D, l, l_P2, D2, l2; par, type)
        res[8] = FOC_data_sharing2(l_P, D, D_I, l, l_P2, D2, l2; par, type)
    end
    return res # Try to force the algorithm to come up with something more precise
end


"""
    Constrained residuals, D_S = 0, no data selling.
"""
function residuals_con(x; par::Pars_v3, type::Symbol)
    l_P = x[1]^2
    l_G = x[2]^2
    D_E = x[3]^2

    D_I = data_gen(l_G; par, type)
    D = bundle(D_I, D_E; par, type)
    l = labor_demand(l_P, D; par, type)

    #! FOCs
    res = zeros(3)
    res[1] = FOC_data_buying(l_P, D, l, D_E; par, type)
    res[2] = FOC_data_gen(l_G, l_P, D, l, D_I; par, type)
    res[3] = FOC_data_proc(l_P, D, l; par, type)

    return sum(abs.(res)) # Try to force the algorithm to come up with something more precise
end

"""
    Constrained residuals, D_S = 0, no data selling.
"""
function residuals_con_nl(x, res; par::Pars_v3, type::Symbol)
    l_G = x[1]^2
    D_E = x[2]^2

    D_I = data_gen(l_G; par, type)
    D = bundle(D_I, D_E; par, type)
    l = labor_demand(D; par, type)

    #! FOCs
    res[1] = FOC_data_buying(D, l, D_E; par, type)
    res[2] = FOC_data_gen(l_G, D, l, D_I; par, type)
    #res[3] = FOC_data_proc(D, l; par, type)

    return res
end

"""
    Main residuals for finding the general equilibrium steady state
"""
function residuals_SS(x; par::Pars_v3)
    w_guess = x[1]^2
    w_P_guess = x[2]^2
    w_G_guess = x[3]^2

    @set! par.w = w_guess
    @set! par.w_P = w_P_guess
    @set! par.w_G = w_G_guess
    @set! par.p_D = find_p_D(par)

    l_P_BS, l_G_BS, l_BS, D_E_BS, D_I_BS, D_BS = sol_f_problem(par; type = :BS)
    l_P_LS, l_G_LS, l_LS, D_E_LS, D_I_LS, D_LS = sol_f_problem(par; type = :LS)
    l_P_LU, l_G_LU, l_LU, D_E_LU, D_I_LU, D_LU = sol_f_problem(par; type = :LU)
    l_P_BU, l_G_BU, l_BU, D_E_BU, D_I_BU, D_BU = sol_f_problem(par; type = :BU)

    Π_BS = firm_revenue(l_P_BS, D_BS, l_BS; par, type = :BS)
    Π_LS = firm_revenue(l_P_LS, D_LS, l_LS; par, type = :LS)
    Π_LU = firm_revenue(l_P_LU, D_LU, l_LU; par, type = :LU)
    Π_BU = firm_revenue(l_P_BU, D_BU, l_BU; par, type = :BU)

    D = [D_BS, D_LS, D_LU, D_BU]
    D_I = [D_I_BS, D_I_LS, D_I_LU, D_I_BU]
    Π = [Π_BS, Π_LS, Π_LU, Π_BU]
    l_G = [l_G_BS, l_G_LS, l_G_LU, l_G_BU]

    w = wage_int_goods(Π; par)
    w_G = wage_gen_data(Π, D, l_G, D_I; par)
    w_P = wage_gen_proc(Π, D; par)

    res = zeros(3)

    res[1] = w_guess - w
    res[2] = w_P_guess - w_P
    res[3] = w_G_guess - w_G

    return sum(res .^ 2)
end


"""
    Including a guess on Y.
"""
function residuals_SS2(x; par::Pars_v3)
    Y_guess = x[1]^2
    w_guess = x[2]^2
    w_P_guess = x[3]^2
    w_G_guess = x[4]^2
    #p_D_guess = x[5]^2

    @set! par.Y = Y_guess
    @set! par.w = w_guess
    @set! par.w_P = w_P_guess
    @set! par.w_G = w_G_guess
    @set! par.p_D = find_p_D(par)

    sol_BS = solve_firm_problem(par; type = :BS)
    sol_LS = solve_firm_problem(par; type = :LS)
    sol_LU = solve_firm_problem(par; type = :LU)
    sol_BU = solve_firm_problem(par; type = :BU)

    l_P_BS, l_G_BS, l_BS, D_E_BS, D_I_BS, D_BS = get_solution(sol_BS; par, type = :BS)
    l_P_LS, l_G_LS, l_LS, D_E_LS, D_I_LS, D_LS = get_solution(sol_LS; par, type = :LS)
    l_P_LU, l_G_LU, l_LU, D_E_LU, D_I_LU, D_LU = get_solution(sol_LU; par, type = :LU)
    l_P_BU, l_G_BU, l_BU, D_E_BU, D_I_BU, D_BU = get_solution(sol_BU; par, type = :BU)

    #D_BS = bundle(D_I_BS, D_E_BS; par)
    #D_LS = bundle(D_I_LS, D_E_LS; par)
    #D_LU = bundle(D_I_LU, D_E_LU; par)
    #D_BU = bundle(D_I_BU, D_E_BU; par)

    Y_BS = firm_output(l_P_BS, D_BS, l_BS; par, type = :BS)
    Y_LS = firm_output(l_P_LS, D_LS, l_LS; par, type = :LS)
    Y_LU = firm_output(l_P_LU, D_LU, l_LU; par, type = :LU)
    Y_BU = firm_output(l_P_BU, D_BU, l_BU; par, type = :BU)

    Y = aggregate_output(Y_BS, Y_LS, Y_LU, Y_BU; par)

    res = zeros(4)

    res[1] = (Y - Y_guess)
    res[2] = labor_market_clearing(l_BS, l_LS, l_LU, l_BU; par)
    res[3] = labor_market_clearing(l_P_BS, l_P_LS, l_P_LU, l_P_BU; par)
    res[4] = labor_market_clearing(l_G_BS, l_G_LS, l_G_LU, l_G_BU; par)

    # l_G = [l_G_BS, l_G_LS, l_G_LU, l_G_BU]
    # D_E = [D_E_BS, D_E_LS, D_E_LU, D_E_BU]
    # D_I = [D_I_BS, D_I_LS, D_I_LU, D_I_BU]

    # res[5] = data_market_clearing(l_G, D_E, D_I; par)

    return sum(res .^ 2)
end


"""
    Data market clearing condition
"""
function data_market_clearing(l_G, D_E, D_I; par::Pars_v3)
    l_G_BS = l_G[1]
    l_G_LS = l_G[2]
    l_G_LU = l_G[3]
    l_G_BU = l_G[4]

    D_E_BS = D_E[1]
    D_E_LS = D_E[2]
    D_E_LU = D_E[3]
    D_E_BU = D_E[4]

    D_I_BS = D_I[1]
    D_I_LS = D_I[2]
    D_I_LU = D_I[3]
    D_I_BU = D_I[4]

    D_S_BS = sold_data(D_I_BS, l_G_BS; par, type = :BS)
    D_S_LS = sold_data(D_I_LS, l_G_LS; par, type = :LS)
    D_S_LU = sold_data(D_I_LU, l_G_LU; par, type = :LU)
    D_S_BU = sold_data(D_I_BU, l_G_BU; par, type = :BU)

    w_BS = weight_of_type(:BS; par)
    w_LS = weight_of_type(:LS; par)
    w_LU = weight_of_type(:LU; par)
    w_BU = weight_of_type(:BU; par)

    D_S = w_BS * D_S_BS + w_LS * D_S_LS + w_LU * D_S_LU + w_BU * D_S_BU
    D_B = w_BS * D_E_BS + w_LS * D_E_LS + w_LU * D_E_LU + w_BU * D_E_BU
    #! 1-τ times all data sold is equal to what is bought
    return D_S * (1 - par.τ) - D_B
end

function residuals_p_D(p_D; par::Pars_v3)

    @set! par.p_D = p_D

    l_P_BS, l_G_BS, l_BS, D_E_BS, D_I_BS, D_BS = sol_f_problem_nl(par, type = :BS)
    l_P_LS, l_G_LS, l_LS, D_E_LS, D_I_LS, D_LS = sol_f_problem_nl(par, type = :LS)
    l_P_LU, l_G_LU, l_LU, D_E_LU, D_I_LU, D_LU = sol_f_problem_nl(par, type = :LU)
    l_P_BU, l_G_BU, l_BU, D_E_BU, D_I_BU, D_BU = sol_f_problem_nl(par, type = :BU)

    l_G = [l_G_BS, l_G_LS, l_G_LU, l_G_BU]
    D_E = [D_E_BS, D_E_LS, D_E_LU, D_E_BU]
    D_I = [D_I_BS, D_I_LS, D_I_LU, D_I_BU]

    data_market_clearing(l_G, D_E, D_I; par)
end