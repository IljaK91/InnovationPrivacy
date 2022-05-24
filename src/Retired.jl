function residuals_SS_bounded(x; par::Pars_v3)
    Y_guess = x[1]
    w_guess = x[2]
    w_P_guess = x[3]
    w_G_guess = x[4]
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
    @show diff = sum(res .^ 2)
    return sum(res .^ 2)
end

# function residuals_SS_single(x; par::Pars_v3)
#     Y_guess = x[1]^2
#     w_guess = x[2]^2
#     w_P_guess = x[3]^2
#     w_G_guess = x[4]^2
#     p_D_guess = x[5]^2

#     @set! par.Y = Y_guess
#     @set! par.w = w_guess
#     @set! par.w_P = w_P_guess
#     @set! par.w_G = w_G_guess
#     @set! par.p_D = p_D_guess

#     sol_BS = solve_firm_problem(par; type = :BS)
#     sol_LS = solve_firm_problem(par; type = :LS)
#     sol_LU = solve_firm_problem(par; type = :LU)
#     sol_BU = solve_firm_problem(par; type = :BU)

#     l_P_BS, l_G_BS, l_BS, D_E_BS, D_I_BS, D_BS = get_solution(sol_BS; par, type =:BS)
#     l_P_LS, l_G_LS, l_LS, D_E_LS, D_I_LS, D_LS = get_solution(sol_LS; par, type =:LS)
#     l_P_LU, l_G_LU, l_LU, D_E_LU, D_I_LU, D_LU = get_solution(sol_LU; par, type =:LU)
#     l_P_BU, l_G_BU, l_BU, D_E_BU, D_I_BU, D_BU = get_solution(sol_BU; par, type =:BU)

#     #D_BS = bundle(D_I_BS, D_E_BS; par)
#     #D_LS = bundle(D_I_LS, D_E_LS; par)
#     #D_LU = bundle(D_I_LU, D_E_LU; par)
#     #D_BU = bundle(D_I_BU, D_E_BU; par)

#     Y_BS = firm_output(l_P_BS, D_BS, l_BS; par, type = :BS)
#     Y_LS = firm_output(l_P_LS, D_LS, l_LS; par, type = :LS)
#     Y_LU = firm_output(l_P_LU, D_LU, l_LU; par, type = :LU)
#     Y_BU = firm_output(l_P_BU, D_BU, l_BU; par, type = :BU)

#     Y = aggregate_output(Y_BS, Y_LS, Y_LU, Y_BU; par)

#     res = zeros(5)

#     res[1] = (Y - Y_guess)
#     res[2] = labor_market_clearing(l_BS, l_LS, l_LU, l_BU; par)
#     res[3] = labor_market_clearing(l_P_BS, l_P_LS, l_P_LU, l_P_BU; par)
#     res[4] = labor_market_clearing(l_G_BS, l_G_LS, l_G_LU, l_G_BU; par)

#     l_G = [l_G_BS, l_G_LS, l_G_LU, l_G_BU]
#     D_E = [D_E_BS, D_E_LS, D_E_LU, D_E_BU]
#     D_I = [D_I_BS, D_I_LS, D_I_LU, D_I_BU]

#     res[5] = data_market_clearing(l_G, D_E, D_I; par)

#     return res
# end

"""
    Solving the model with a general solution algorithm does not work performantly.
"""
function find_steady_state(par::Pars_v3)
    f(x) = residuals_SS(x; par)
    initial_x = [par.Y, par.w, par.w_P, par.w_G]
    optimize(f, initial_x)
end


function constraint_data_selling(D_S)
    if D_S > 0
        return 0
    elseif D_S <= 0
        100000 * D_S^2
    else
        error("I should never be here")
    end
end


function wage_int_goods2(K; par)
    K_BS = K[1]
    K_LS = K[2]
    K_LU = K[3]
    K_BU = K[4]

    w_BS = weight_of_type(:BS; par)
    w_LS = weight_of_type(:LS; par)
    w_LU = weight_of_type(:LU; par)
    w_BU = weight_of_type(:BU; par)

    return par.Y^par.α_Y * par.α_L_hat * (w_BS * K_BS + w_LS * K_LS + w_LU * K_LU + w_BU * K_BU)
end



"""
    Wage for data generation
"""
function wage_gen_data2(Π, D; par)
    Π_BS = Π[1] # firm revenue
    Π_LS = Π[2] # firm revenue
    Π_LU = Π[3] # firm revenue
    Π_BU = Π[4] # firm revenue

    D_BS = D[1] # data bundle
    D_LS = D[2] # data bundle
    D_LU = D[3] # data bundle
    D_BU = D[4] # data bundle

    w_BS = weight_of_type(:BS; par) # weights
    w_LS = weight_of_type(:LS; par) # weights
    w_LU = weight_of_type(:LU; par) # weights
    w_BU = weight_of_type(:BU; par) # weights

    A_G_BS = A_G_of_type(:BS; par) # prod of data generation
    A_G_LS = A_G_of_type(:LS; par) # prod of data generation
    A_G_LU = A_G_of_type(:LU; par) # prod of data generation
    A_G_BU = A_G_of_type(:BU; par) # prod of data generation

    c_BS = Π_BS / D_BS * A_G_BS # Put some things together
    c_LS = Π_LS / D_LS * A_G_LS # Put some things together
    c_LU = Π_LU / D_LU * A_G_LU # Put some things together
    c_BU = Π_BU / D_BU * A_G_BU # Put some things together

    return par.α_K_hat * ((1 - par.γ_S) * w_BS * c_BS + (1 - par.γ_S) * w_LS * c_LS + (1 - par.γ_U) * w_LU * c_LU + (1 - par.γ_U) * w_BU * c_BU)
end

"""
    Return a vector containing all types
"""
function weights_of_types(par::Pars_v3)
    w_BS = weight_of_type(:BS; par)
    w_LS = weight_of_type(:LS; par)
    w_LU = weight_of_type(:LU; par)
    w_BU = weight_of_type(:BU; par)

    return [w_BS, w_LS, w_LU, w_BU]
end