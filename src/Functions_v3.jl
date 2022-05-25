#! Here come all the functions for the model with trade in structured data and non-rival data

"""
    Solve the problem of an unconstrained firm that can choose D_S freely (should stay positive).
"""
function solve_firm_problem_unc(par::Pars_v3; type::Symbol)
    f(x) = residuals_unc(x; par, type)
    if type == :BS
        initial_x = [1.090961407034043,
            1.7305342659571796,
            0.8017602275570367,
            0.846302424724036]
    elseif type == :BU
        initial_x = [1.2797075968276403,
            1.7305374802420503,
            0.7399594318155204,
            0.7810681256207107]
    elseif type == :LS
        initial_x = [1.0909610262794553,
            -0.11836701171722647,
            0.8017601473223311,
            0.8463025640344691]
    elseif type == :LU
        initial_x = [1.2797075968276403,
            1.7305374802420503,
            0.7399594318155204,
            0.7810681256207107]
    end
    #res = optimize(f, initial_x, Optim.Options(g_tol = 1e-14))
    res = optimize(f, initial_x)
    # try
    #     @assert Optim.converged(res) == true
    # catch
    #     τ = par.τ
    #     println("Solution may be inaccurate for $τ and $type")
    # end
    return res
end

function solve_firm_problem_unc_nl(par::Pars_v3; type::Symbol)
    res = zeros(3)
    f!(x) = residuals_unc_nl(x, res; par, type)

    initial_x = [1.090961407034043,
        1.7305342659571796,
        0.8017602275570367]

    #res = optimize(f, initial_x, Optim.Options(g_tol = 1e-14))
    res = nlsolve(f!, initial_x)
    # try
    #     @assert Optim.converged(res) == true
    # catch
    #     τ = par.τ
    #     println("Solution may be inaccurate for $τ and $type")
    # end
    return res
end

function solve_firm_problem_unc_nl2(par::Pars_v3; type::Symbol, sec_period::Symbol=:no)
    if sec_period == :yes
        res = zeros(4)
        initial_x = [1.090961407034043,
            1.7305342659571796,
            0.8017602275570367,
            0.846302424724036]
    elseif sec_period == :no
        res = zeros(8)
        initial_x = [1.090961407034043,
            1.7305342659571796,
            0.8017602275570367,
            0.846302424724036,
            1.090961407034043,
            1.7305342659571796,
            0.8017602275570367,
            0.846302424724036]
    end
    f!(x) = residuals_unc_nl2(x, res; par, type, sec_period)
    res = nlsolve(f!, initial_x)
    return res
end

"""
    Solve the problem of a constrained firm with D_S = 0 (would like to buy internal data).
"""
function solve_firm_problem_con_nl(par::Pars_v3; type::Symbol)
    #f(x) = residuals_con(x; par, type)
    res = zeros(2)
    f!(x) = residuals_con_nl(x, res; par, type)
    initial_x = [0.5, 0.5]
    #res = optimize(f, initial_x, Optim.Options(g_tol = 1e-14))
    res = nlsolve(f!, initial_x)
    #@assert Optim.converged(res) == true
    return res
end

function solve_firm_problem_con(par::Pars_v3; type::Symbol)
    f(x) = residuals_con(x; par, type)
    initial_x = [0.5, 0.5, 0.5]
    #res = optimize(f, initial_x, Optim.Options(g_tol = 1e-14))
    res = optimize(f, initial_x)
    #@assert Optim.converged(res) == true
    return res
end

"""
    Solve the problem of the firm.
"""
function solve_firm_problem(par::Pars_v3; type::Symbol)
    sol = solve_firm_problem_unc(par; type)
    l_P, l_G, l, D_E, D_I, D = get_solution(sol; par, type)

    D_S = sold_data(D_I, l_G; par, type)

    if D_S < 0
        sol_con = solve_firm_problem_con(par; type)
        return sol_con
    elseif D_S >= 0
        return sol
    else
        error("I never should be here")
    end
end

function solve_firm_problem_nl(par::Pars_v3; type::Symbol)
    sol = solve_firm_problem_unc_nl(par; type)
    l_G, l, D_E, D_I, D = get_solution_nl_unc(sol; par, type)

    D_S = sold_data(D_I, l_G; par, type)
    #! Idea: Solve first unconstrained problem. If firm tries to buy internal data, constrain to not sharing data.
    if D_S < 0
        sol_con = solve_firm_problem_con_nl(par; type)
        return get_solution_nl_con(sol_con; par, type)
    elseif D_S >= 0
        return get_solution_nl_unc(sol; par, type)
    else
        error("I never should be here")
    end
end

# function firm_profits(x; par::Pars_v3, type::Symbol)
#     l_P = x[1]^2
#     l_G = x[2]^2
#     l   = x[3]^2
#     D_E = x[4]^2
#     D_S = x[5]

#     # All other variables we can directly constrain to be 

# end

"""
    Returns the labor market clearing condition residuals
"""
function labor_market_clearing(l_BS, l_LS, l_LU, l_BU; par)
    w_BS = weight_of_type(:BS; par)
    w_LS = weight_of_type(:LS; par)
    w_LU = weight_of_type(:LU; par)
    w_BU = weight_of_type(:BU; par)

    #! Labor of each type is fixed to one for the moment.
    return l_BS * w_BS + l_LS * w_LS + l_LU * w_LU + l_BU * w_BU - 1
end

"""
    Aggregate Production Function.
"""
function aggregate_output(Y_BS, Y_LS, Y_LU, Y_BU; par::Pars_v3)
    w_BS = weight_of_type(:BS; par)
    w_LS = weight_of_type(:LS; par)
    w_LU = weight_of_type(:LU; par)
    w_BU = weight_of_type(:BU; par)

    c = (par.σ - 1) / par.σ

    if par.compl_prod == :yes
        # Aggregat Production function with strategic complementarities in production
        # May lead to instability, could not solve the model so far with this specification.
        Y_agg = (w_BS * Y_BS^c + w_LS * Y_LS^c + w_LU * Y_LU^c + w_BU * Y_BU^c)^(1 / c)
    elseif par.compl_prod == :no
        # Aggregat Production function without strategic complementarities in production.
        # Demand for one good is independent of demand for other goods.
        Y_agg = (w_BS * Y_BS^c + w_LS * Y_LS^c + w_LU * Y_LU^c + w_BU * Y_BU^c) / c
    else
        error("I should never be here!")
    end
end

function aggregate_labor(l_BS, l_LS, l_LU, l_BU; par::Pars_v3)
    w_BS = weight_of_type(:BS; par)
    w_LS = weight_of_type(:LS; par)
    w_LU = weight_of_type(:LU; par)
    w_BU = weight_of_type(:BU; par)

    return w_BS * l_BS + w_LS * l_LS + w_LU * l_LU + w_BU * l_BU
end

sol_f_problem(par::Pars_v3; type::Symbol) = get_solution(solve_firm_problem(par; type); par, type)
sol_f_problem_nl(par::Pars_v3; type::Symbol) = solve_firm_problem_nl(par; type)

function get_ss_solution(sol)
    Y = sol.minimizer[1]^2
    w = sol.minimizer[2]^2
    w_P = sol.minimizer[3]^2
    w_G = sol.minimizer[4]^2
    p_D = sol.minimizer[5]^2

    return Y, w, w_P, w_G, p_D
end

"""
    Find the steady state given a set of parameters
"""
function find_steady_state(par::Pars_v3; tol=1e-10, show_steps::Symbol=:yes)
    diff = 10
    counter = 0
    while diff > tol
        counter += 1
        p_D_sol = find_p_D(par)
        @set! par.p_D = p_D_sol

        l_G_BS, l_BS, D_E_BS, D_I_BS, D_BS = sol_f_problem_nl(par; type=:BS)
        l_G_LS, l_LS, D_E_LS, D_I_LS, D_LS = sol_f_problem_nl(par; type=:LS)
        l_G_LU, l_LU, D_E_LU, D_I_LU, D_LU = sol_f_problem_nl(par; type=:LU)
        l_G_BU, l_BU, D_E_BU, D_I_BU, D_BU = sol_f_problem_nl(par; type=:BU)

        l_G = [l_G_BS, l_G_LS, l_G_LU, l_G_BU]
        l = [l_BS, l_LS, l_LU, l_BU]
        D_E = [D_E_BS, D_E_LS, D_E_LU, D_E_BU]
        D_I = [D_I_BS, D_I_LS, D_I_LU, D_I_BU]
        D = [D_BS, D_LS, D_LU, D_BU]

        #D_BS = bundle(D_I_BS, D_E_BS; par)
        #D_LS = bundle(D_I_LS, D_E_LS; par)
        #D_LU = bundle(D_I_LU, D_E_LU; par)
        #D_BU = bundle(D_I_BU, D_E_BU; par)

        Y_BS = firm_output(D_BS, l_BS; par, type=:BS)
        Y_LS = firm_output(D_LS, l_LS; par, type=:LS)
        Y_LU = firm_output(D_LU, l_LU; par, type=:LU)
        Y_BU = firm_output(D_BU, l_BU; par, type=:BU)

        Π_BS = firm_revenue(D_BS, l_BS; par, type=:BS)
        Π_LS = firm_revenue(D_LS, l_LS; par, type=:LS)
        Π_LU = firm_revenue(D_LU, l_LU; par, type=:LU)
        Π_BU = firm_revenue(D_BU, l_BU; par, type=:BU)
        
        Π = [Π_BS, Π_LS, Π_LU, Π_BU]

        Y = aggregate_output(Y_BS, Y_LS, Y_LU, Y_BU; par)
        w = wage_int_goods(Π; par)
        w_G = wage_gen_data(Π, D, l_G, D_I; par)

        if par.α_Y == 0
            res = zeros(3)

            res[1] = w - par.w
            res[2] = 10 * (w_G - par.w_G) # Make this constraint more important?


            diff = sum(res .^ 2)
        else
            par.α_Y != 0.0
            res = zeros(4)

            res[1] = Y - par.Y
            res[2] = w - par.w
            res[3] = 10 * (w_G - par.w_G) # Make this constraint more important?

            diff = sum(res .^ 2)
        end

        if show_steps == :yes
            #@show Y
            #@show w
            #@show w_G
            #@show w_P
            println("Iteration: $counter, Residual: $diff")
        end

        #! Try to make the convergence a bit more smooth, choose a small lambda
        @set! par.Y = par.λ * Y + (1 - par.λ) * par.Y
        @set! par.w = par.λ * w + (1 - par.λ) * par.w
        @set! par.w_G = par.λ * w_G + (1 - par.λ) * w_G

        if diff < tol || counter > 100
        
            #! Check the solution with other equations! Did I write the equations for wages correctly?
            # try
            #     @assert abs(aggregate_labor(l_P_BS, l_P_LS, l_P_LU, l_P_BU; par) - par.P) < 1e-2
            #     @assert abs(aggregate_labor(l_G_BS, l_G_LS, l_G_LU, l_G_BU; par) - par.G) < 1e-2
            #     @assert abs(aggregate_labor(l_BS, l_LS, l_LU, l_BU; par) - par.L) < 1e-2
            # catch
            #     @show aggregate_labor(l_P_BS, l_P_LS, l_P_LU, l_P_BU; par)
            #     @show aggregate_labor(l_G_BS, l_G_LS, l_G_LU, l_G_BU; par)
            #     @show aggregate_labor(l_BS, l_LS, l_LU, l_BU; par)
        
            #     error("Labor Market Clearing did not work!")
            # end
        
            #! If convergence has been achieved, save the solution
            D_S_BS = sold_data(D_I_BS, l_G_BS; par, type=:BS)
            D_S_LS = sold_data(D_I_LS, l_G_LS; par, type=:LS)
            D_S_LU = sold_data(D_I_LU, l_G_LU; par, type=:LU)
            D_S_BU = sold_data(D_I_BU, l_G_BU; par, type=:BU)
        
            D_G_BS = data_gen(l_G_BS; par, type=:BS)
            D_G_LS = data_gen(l_G_LS; par, type=:LS)
            D_G_LU = data_gen(l_G_LU; par, type=:LU)
            D_G_BU = data_gen(l_G_BU; par, type=:BU)
        
            prof_BS = firm_profits(l_BS, l_G_BS, D_S_BS, D_E_BS; par, type=:BS)
            prof_LS = firm_profits(l_LS, l_G_LS, D_S_LS, D_E_LS; par, type=:LS)
            prof_LU = firm_profits(l_LU, l_G_LU, D_S_LU, D_E_LU; par, type=:LU)
            prof_BU = firm_profits(l_BU, l_G_BU, D_S_BU, D_E_BU; par, type=:BU)
        
            agg_D_S = D_S_BS + D_S_LS + D_S_LU + D_S_BU
            agg_D_G = D_G_BS + D_G_LS + D_G_LU + D_G_BU
            agg_D = D_BS + D_LS + D_LU + D_BU
        
            @set! par.D_G_SS = agg_D_G
            @set! par.D_S_SS = agg_D_S
            @set! par.Ω = data_multiplier(agg_D_G, agg_D_S; par)
            @set! par.Ω_bundle = data_multiplier_bundle(agg_D, agg_D_G)
        
            @set! par.prof_LU_ss = prof_LU # Labor good production small customer base, unsophisticated
            @set! par.prof_BU_ss = prof_BU # Labor good production big customer base, unsophisticated
            @set! par.prof_LS_ss = prof_LS # Labor good production small customer base, sophisticated
            @set! par.prof_BS_ss = prof_BS # Labor good production big customer base, sophisticated
        
            @set! par.l_LU_ss = l_LU # Labor good production small customer base, unsophisticated
            @set! par.l_BU_ss = l_BU # Labor good production big customer base, unsophisticated
            @set! par.l_LS_ss = l_LS # Labor good production small customer base, sophisticated
            @set! par.l_BS_ss = l_BS # Labor good production big customer base, sophisticated
        
            @set! par.l_G_LU_ss = l_G_LU # Labor good production small customer base, unsophisticated
            @set! par.l_G_BU_ss = l_G_BU # Labor good production big customer base, unsophisticated
            @set! par.l_G_LS_ss = l_G_LS # Labor good production small customer base, sophisticated
            @set! par.l_G_BS_ss = l_G_BS # Labor good production big customer base, sophisticated
        
            @set! par.D_S_LU_ss = D_S_LU # data sharing small customer base, unsophisticated
            @set! par.D_S_BU_ss = D_S_BU # data sharing big customer base, unsophisticated
            @set! par.D_S_LS_ss = D_S_LS # data sharing small customer base, sophisticated
            @set! par.D_S_BS_ss = D_S_BS # data sharing big customer base, sophisticated
        
            @set! par.D_G_LU_ss = D_G_LU # data generated small customer base, unsophisticated
            @set! par.D_G_BU_ss = D_G_BU # data generated big customer base, unsophisticated
            @set! par.D_G_LS_ss = D_G_LS # data generated small customer base, sophisticated
            @set! par.D_G_BS_ss = D_G_BS # data generated big customer base, sophisticated
        
            @set! par.D_LU_ss = D_LU # data generated small customer base, unsophisticated
            @set! par.D_BU_ss = D_BU # data generated big customer base, unsophisticated
            @set! par.D_LS_ss = D_LS # data generated small customer base, sophisticated
            @set! par.D_BS_ss = D_BS # data generated big customer base, sophisticated
        
            @set! par.D_E_LU_ss = D_E_LU # data generated small customer base, unsophisticated
            @set! par.D_E_BU_ss = D_E_BU # data generated big customer base, unsophisticated
            @set! par.D_E_LS_ss = D_E_LS # data generated small customer base, sophisticated
            @set! par.D_E_BS_ss = D_E_BS # data generated big customer base, sophisticated
        
            @set! par.D_I_LU_ss = D_I_LU # data generated small customer base, unsophisticated
            @set! par.D_I_BU_ss = D_I_BU # data generated big customer base, unsophisticated
            @set! par.D_I_LS_ss = D_I_LS # data generated small customer base, sophisticated
            @set! par.D_I_BS_ss = D_I_BS # data generated big customer base, sophisticated
        
            @set! par.Firm_sol = Save_Firm_Solution(l, l_G, D_I, D_E, D; par)
        end
    end
    return par
end

"""
    Solve all firm problems 
"""
function find_steady_state_PE(par::Pars_v3; solve_price::Symbol = :yes)
    if solve_price == :yes
        p_D_sol = find_p_D(par)
        @set! par.p_D = p_D_sol
    end

    l_P_BS, l_G_BS, l_BS, D_E_BS, D_I_BS, D_BS = sol_f_problem_nl(par; type = :BS)
    l_P_LS, l_G_LS, l_LS, D_E_LS, D_I_LS, D_LS = sol_f_problem_nl(par; type = :LS)
    l_P_LU, l_G_LU, l_LU, D_E_LU, D_I_LU, D_LU = sol_f_problem_nl(par; type = :LU)
    l_P_BU, l_G_BU, l_BU, D_E_BU, D_I_BU, D_BU = sol_f_problem_nl(par; type = :BU)

    D_S_BS = sold_data(D_I_BS, l_G_BS; par, type = :BS)
    D_S_LS = sold_data(D_I_LS, l_G_LS; par, type = :LS)
    D_S_LU = sold_data(D_I_LU, l_G_LU; par, type = :LU)
    D_S_BU = sold_data(D_I_BU, l_G_BU; par, type = :BU)

    D_G_BS = data_gen(l_G_BS; par, type = :BS)
    D_G_LS = data_gen(l_G_LS; par, type = :LS)
    D_G_LU = data_gen(l_G_LU; par, type = :LU)
    D_G_BU = data_gen(l_G_BU; par, type = :BU)

    K_BS = knowledge(l_P_BS, D_BS; par, type = :BS) 
    K_LS = knowledge(l_P_LS, D_LS; par, type = :LS) 
    K_LU = knowledge(l_P_LU, D_LU; par, type = :LU) 
    K_BU = knowledge(l_P_BU, D_BU; par, type = :BU)

    prof_BS = firm_profits(l_BS, l_P_BS, l_G_BS, D_S_BS, D_E_BS; par, type = :BS) 
    prof_LS = firm_profits(l_LS, l_P_LS, l_G_LS, D_S_LS, D_E_LS; par, type = :LS) 
    prof_LU = firm_profits(l_LU, l_P_LU, l_G_LU, D_S_LU, D_E_LU; par, type = :LU) 
    prof_BU = firm_profits(l_BU, l_P_BU, l_G_BU, D_S_BU, D_E_BU; par, type = :BU) 

    @set! par.l_LU_ss = l_LU # Labor good production small customer base, unsophisticated
    @set! par.l_BU_ss = l_BU # Labor good production big customer base, unsophisticated
    @set! par.l_LS_ss = l_LS # Labor good production small customer base, sophisticated
    @set! par.l_BS_ss = l_BS # Labor good production big customer base, sophisticated
    
    @set! par.prof_LU_ss = prof_LU # Labor good production small customer base, unsophisticated
    @set! par.prof_BU_ss = prof_BU # Labor good production big customer base, unsophisticated
    @set! par.prof_LS_ss = prof_LS # Labor good production small customer base, sophisticated
    @set! par.prof_BS_ss = prof_BS # Labor good production big customer base, sophisticated
    
    @set! par.l_P_LU_ss = l_P_LU # Labor good production small customer base, unsophisticated
    @set! par.l_P_BU_ss = l_P_BU # Labor good production big customer base, unsophisticated
    @set! par.l_P_LS_ss = l_P_LS # Labor good production small customer base, sophisticated
    @set! par.l_P_BS_ss = l_P_BS # Labor good production big customer base, sophisticated

    @set! par.l_G_LU_ss = l_G_LU # Labor good production small customer base, unsophisticated
    @set! par.l_G_BU_ss = l_G_BU # Labor good production big customer base, unsophisticated
    @set! par.l_G_LS_ss = l_G_LS # Labor good production small customer base, sophisticated
    @set! par.l_G_BS_ss = l_G_BS # Labor good production big customer base, sophisticated

    @set! par.D_S_LU_ss = D_S_LU # data sharing small customer base, unsophisticated
    @set! par.D_S_BU_ss = D_S_BU # data sharing big customer base, unsophisticated
    @set! par.D_S_LS_ss = D_S_LS # data sharing small customer base, sophisticated
    @set! par.D_S_BS_ss = D_S_BS # data sharing big customer base, sophisticated

    @set! par.D_G_LU_ss = D_G_LU # data generated small customer base, unsophisticated
    @set! par.D_G_BU_ss = D_G_BU # data generated big customer base, unsophisticated
    @set! par.D_G_LS_ss = D_G_LS # data generated small customer base, sophisticated
    @set! par.D_G_BS_ss = D_G_BS # data generated big customer base, sophisticated

    @set! par.D_LU_ss = D_LU # data generated small customer base, unsophisticated
    @set! par.D_BU_ss = D_BU # data generated big customer base, unsophisticated
    @set! par.D_LS_ss = D_LS # data generated small customer base, sophisticated
    @set! par.D_BS_ss = D_BS # data generated big customer base, sophisticated

    @set! par.D_E_LU_ss = D_E_LU # data generated small customer base, unsophisticated
    @set! par.D_E_BU_ss = D_E_BU # data generated big customer base, unsophisticated
    @set! par.D_E_LS_ss = D_E_LS # data generated small customer base, sophisticated
    @set! par.D_E_BS_ss = D_E_BS # data generated big customer base, sophisticated

    @set! par.D_I_LU_ss = D_I_LU # data generated small customer base, unsophisticated
    @set! par.D_I_BU_ss = D_I_BU # data generated big customer base, unsophisticated
    @set! par.D_I_LS_ss = D_I_LS # data generated small customer base, sophisticated
    @set! par.D_I_BS_ss = D_I_BS # data generated big customer base, sophisticated

    @set! par.K_LU_ss = K_LU # knowledge employed small customer base, unsophisticated
    @set! par.K_BU_ss = K_BU # knowledge employed big customer base, unsophisticated
    @set! par.K_LS_ss = K_LS # knowledge employed small customer base, sophisticated
    @set! par.K_BS_ss = K_BS # knowledge employed big customer base, sophisticated

    return par
end

"""
    Find first the price for data that clears the data market.
"""
function find_p_D(par::Pars_v3; verbose = false)
    f(p_D) = residuals_p_D(p_D^2; par)
    #f(p_D) = residuals_p_D(p_D; par)
    find_zero(f, 0.3, verbose = verbose)^2
    #find_zero(f, [0.001, 2], verbose = verbose, xatol = 1e-10, xrtol = 1e-10, atol = 1e-10, rtol = 1e-10)
end

"""
    Compute comparatics statics using the simplified version of the model
"""
function comp_statics(τ_set, ν_set, A_G_common_set; tol = 1e-10, ζ = 1.0, α_S = 0.8, α_U = 0.6)
    grid = gridmake(τ_set, ν_set, A_G_common_set)
    sol = Array{Any}(undef, 0)
    for i in 1:size(grid)[1]
        @show i
        @show grid[i, :]
        if i == 1
            push!(sol, find_steady_state(Pars_v3(τ = grid[i, 1], ν = grid[i, 2], A_G_common = grid[i, 3], ζ = ζ, α_S = α_S, α_U = α_U); tol))
        else
            #try
            push!(sol, find_steady_state(Pars_v3(τ = grid[i, 1], ν = grid[i, 2], A_G_common = grid[i, 3], Y = sol[i-1].Y, w = sol[i-1].w, w_G = sol[i-1].w_G, ζ = ζ, α_S = α_S, α_U = α_U); tol))
            #catch
            #    push!(sol, find_steady_state(Pars_v3(τ = grid[i, 1], ν = grid[i, 2], A_G_common = grid[i, 3], A_P_common = grid[i, 4]); tol))
            #end
        end
    end
    return sol
end

function find_value_func(par::Pars_v3; τ)
    @set! par.τ = τ
    @set! par.second_period = :yes
    K_last_set = 0.0:0.1:2
    D_last_set = 0.0:0.1:2
    grid = gridmake(K_last_set, D_last_set)
    V = zeros(size(grid)[1], 4)
    par = find_steady_state(par, tol = 1e-12)
    i = 0
    f(K_last, D_last) = begin
        @set! par.D_last_LU = D_last # data generated small customer base, unsophisticated
        @set! par.D_last_BU = D_last # data generated big customer base, unsophisticated
        @set! par.D_last_LS = D_last # data generated small customer base, sophisticated
        @set! par.D_last_BS = D_last # data generated big customer base, sophisticated

        @set! par.K_last_LU = K_last # knowledge employed small customer base, unsophisticated
        @set! par.K_last_BU = K_last # knowledge employed big customer base, unsophisticated
        @set! par.K_last_LS = K_last # knowledge employed small customer base, sophisticated
        @set! par.K_last_BS = K_last # knowledge employed big customer base, sophisticated

        par = find_steady_state_PE(par, solve_price = :no)

        par.prof_BS_ss
    end

    A = [f(x, y) for x in K_last_set, y in D_last_set]

    # linear interpolation
    return LinearInterpolation((K_last_set, D_last_set), A)

    # cubic spline interpolation
    #interp_cubic = CubicSplineInterpolation((xs, ys), A)
    #DataFrame(K_last = grid[:, 1], D_last = grid[:, 2], V_BS = V[:, 1], V_LS = V[:, 2], V_LU = V[:, 3], V_BU = V[:, 4])
end