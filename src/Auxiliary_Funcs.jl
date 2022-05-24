"""
    Calculate firm variables and populate struct with steady state values.
"""
function Save_Firm_Solution(l, l_P, l_G, D_I, D_E, D; par::Pars_v3)
    type = [:BS, :LS, :LU, :BU]

    D_S = sold_data.(D_I, l_G, type; par)
    D_G = data_gen.(l_G, type; par)
    prof = firm_profits.(l, l_P, l_G, D_S, D_E, type; par)
    K = knowledge.(l_P, D, type; par)
    res = Firm_sol()
    @set! res.prof_BS_ss = prof[1] # Labor good production big customer base, sophisticated
    @set! res.prof_LS_ss = prof[2] # Labor good production small customer base, sophisticated
    @set! res.prof_LU_ss = prof[3] # Labor good production small customer base, unsophisticated
    @set! res.prof_BU_ss = prof[4] # Labor good production big customer base, unsophisticated

    @set! res.l_BS_ss = l[1] # Labor good production big customer base, sophisticated
    @set! res.l_LS_ss = l[2] # Labor good production small customer base, sophisticated
    @set! res.l_LU_ss = l[3] # Labor good production small customer base, unsophisticated
    @set! res.l_BU_ss = l[4] # Labor good production big customer base, unsophisticated

    @set! res.l_P_BS_ss = l_P[1] # Labor good production big customer base, sophisticated
    @set! res.l_P_LS_ss = l_P[2] # Labor good production small customer base, sophisticated
    @set! res.l_P_LU_ss = l_P[3] # Labor good production small customer base, unsophisticated
    @set! res.l_P_BU_ss = l_P[4] # Labor good production big customer base, unsophisticated

    @set! res.l_G_BS_ss = l_G[1] # Labor good production big customer base, sophisticated
    @set! res.l_G_LS_ss = l_G[2] # Labor good production small customer base, sophisticated
    @set! res.l_G_LU_ss = l_G[3] # Labor good production small customer base, unsophisticated
    @set! res.l_G_BU_ss = l_G[4] # Labor good production big customer base, unsophisticated

    @set! res.D_S_BS_ss = D_S[1] # data sharing big customer base, sophisticated
    @set! res.D_S_LS_ss = D_S[2] # data sharing small customer base, sophisticated
    @set! res.D_S_LU_ss = D_S[3] # data sharing small customer base, unsophisticated
    @set! res.D_S_BU_ss = D_S[4] # data sharing big customer base, unsophisticated

    @set! res.D_G_BS_ss = D_G[1] # data generated big customer base, sophisticated
    @set! res.D_G_LS_ss = D_G[2] # data generated small customer base, sophisticated
    @set! res.D_G_LU_ss = D_G[3] # data generated small customer base, unsophisticated
    @set! res.D_G_BU_ss = D_G[4] # data generated big customer base, unsophisticated

    @set! res.D_BS_ss = D[1] # data generated big customer base, sophisticated
    @set! res.D_LS_ss = D[2] # data generated small customer base, sophisticated
    @set! res.D_LU_ss = D[3] # data generated small customer base, unsophisticated
    @set! res.D_BU_ss = D[4] # data generated big customer base, unsophisticated

    @set! res.D_E_BS_ss = D_E[1] # data generated big customer base, sophisticated
    @set! res.D_E_LS_ss = D_E[2] # data generated small customer base, sophisticated
    @set! res.D_E_LU_ss = D_E[3] # data generated small customer base, unsophisticated
    @set! res.D_E_BU_ss = D_E[4] # data generated big customer base, unsophisticated

    @set! res.D_I_BS_ss = D_I[1] # data generated big customer base, sophisticated
    @set! res.D_I_LS_ss = D_I[2] # data generated small customer base, sophisticated
    @set! res.D_I_LU_ss = D_I[3] # data generated small customer base, unsophisticated
    @set! res.D_I_BU_ss = D_I[4] # data generated big customer base, unsophisticated

    @set! res.K_BS_ss = K[1] # knowledge employed big customer base, sophisticated
    @set! res.K_LS_ss = K[2] # knowledge employed small customer base, sophisticated
    @set! res.K_LU_ss = K[3] # knowledge employed small customer base, unsophisticated
    @set! res.K_BU_ss = K[4] # knowledge employed big customer base, unsophisticated

    return res
end

function get_alpha_K_type(type::Symbol;par::Pars_v3)
    if type == :LU # little, unsophisticated
        par.α_U_K_hat
    elseif type == :LS # little, sophisticated
        par.α_S_K_hat
    elseif type == :BU # big, unsophisticated
        par.α_U_K_hat
    elseif type == :BS # big, sophisticated
        par.α_S_K_hat
    else
        error("Firm type needs to be either LU, LS, BU or BS")
    end
end

function get_alpha_L_type(type::Symbol;par::Pars_v3)
    if type == :LU # little, unsophisticated
        par.α_U_L_hat
    elseif type == :LS # little, sophisticated
        par.α_S_L_hat
    elseif type == :BU # big, unsophisticated
        par.α_U_L_hat
    elseif type == :BS # big, sophisticated
        par.α_S_L_hat
    else
        error("Firm type needs to be either LU, LS, BU or BS")
    end
end

function weight_of_type(type::Symbol; par::Pars_v3)
    if type == :LU # little, unsophisticated
        (1 - par.w_B) * (1 - par.w_S)
    elseif type == :LS # little, sophisticated
        (1 - par.w_B) * par.w_S
    elseif type == :BU # big, unsophisticated
        par.w_B * (1 - par.w_S)
    elseif type == :BS # big, sophisticated
        par.w_B * par.w_S
    else
        error("Firm type needs to be either LU, LS, BU or BS")
    end
end

"""
    Data Generation Productivity / Customer Base depending on type.
"""
function A_G_of_type(type::Symbol; par::Pars_v3)
    if type == :LU # little, unsophisticated
        par.A_G_S
    elseif type == :LS # little, sophisticated
        par.A_G_S
    elseif type == :BU # big, unsophisticated
        par.A_G_B
    elseif type == :BS # big, sophisticated
        par.A_G_B
    else
        error("Firm type needs to be either LU, LS, BU or BS")
    end
end

"""
    Data Processing Productivity depending on type.
"""
function A_P_of_type(type::Symbol; par::Pars_v3)
    if type == :LU # little, unsophisticated
        par.A_P_U
    elseif type == :LS # little, sophisticated
        par.A_P_S
    elseif type == :BU # big, unsophisticated
        par.A_P_U
    elseif type == :BS # big, sophisticated
        par.A_P_S
    else
        error("Firm type needs to be either LU, LS, BU or BS")
    end
end


"""
    Squares the minimizer and returns the solution values
"""
function get_solution(sol; par::Pars_v3, type::Symbol)
    if length(sol.minimizer) == 4
        l_P = sol.minimizer[1]^2
        l_G = sol.minimizer[2]^2
        D_E = sol.minimizer[3]^2
        D_I = sol.minimizer[4]^2
        D = bundle(D_I, D_E; par, type)
        l = labor_demand(l_P, D; par, type)
    elseif length(sol.minimizer) == 3
        l_P = sol.minimizer[1]^2
        l_G = sol.minimizer[2]^2
        D_E = sol.minimizer[3]^2
        D_I = data_gen(l_G; par, type)
        D = bundle(D_I, D_E; par, type)
        l = labor_demand(l_P, D; par, type)
    else
        error("I should never be here")
    end

    return l_P, l_G, l, D_E, D_I, D
end

function get_solution_nl(sol; par::Pars_v3, type::Symbol)
    if length(sol.zero) == 4
        l_P = sol.zero[1]^2
        l_G = sol.zero[2]^2
        D_E = sol.zero[3]^2
        D_I = sol.zero[4]^2
        D = bundle(D_I, D_E; par, type)
        l = labor_demand(l_P, D; par, type)
    elseif length(sol.zero) == 3
        l_P = sol.zero[1]^2
        l_G = sol.zero[2]^2
        D_E = sol.zero[3]^2
        D_I = data_gen(l_G; par, type)
        D = bundle(D_I, D_E; par, type)
        l = labor_demand(l_P, D; par, type)
    else
        error("I should never be here")
    end

    return l_P, l_G, l, D_E, D_I, D
end

function get_solution_nl2(sol; par::Pars_v3, type::Symbol)
    if length(sol.zero) == 4
        l_P = sol.zero[1]^2
        l_G = sol.zero[2]^2
        D_E = sol.zero[3]^2
        D_I = sol.zero[4]^2
        D = bundle(D_I, D_E; par, type)
        l = labor_demand(l_P, D; par, type)
        return l_P, l_G, l, D_E, D_I, D
    elseif length(sol.zero) == 3
        l_P = sol.zero[1]^2
        l_G = sol.zero[2]^2
        D_E = sol.zero[3]^2
        D_I = data_gen(l_G; par, type)
        D = bundle(D_I, D_E; par, type)
        l = labor_demand(l_P, D; par, type)
        return l_P, l_G, l, D_E, D_I, D
    elseif length(sol.zero) == 8
        l_P = sol.zero[1]^2
        l_G = sol.zero[2]^2
        D_E = sol.zero[3]^2
        D_I = sol.zero[4]^2
        l_P2 = sol.zero[5]^2
        l_G2 = sol.zero[6]^2
        D_E2 = sol.zero[7]^2
        D_I2 = sol.zero[8]^2
        D = bundle(D_I, D_E; par, type)
        D2 = bundle(D_I2, D_E2; par, type)
        l = labor_demand(l_P, D; par, type, sec_period=:no)
        l2 = labor_demand(l_P2, D2; par, type, sec_period=:yes)
        return l_P, l_G, l, D_E, D_I, D, l_P2, l_G2, l2, D_E2, D_I2, D2
    elseif length(sol.zero) == 6 #! Need to take also into account an intermediate case in which firms switch from data sharing to not sharing... but unlikely to arise in our setting
        l_P = sol.zero[1]^2
        l_G = sol.zero[2]^2
        D_E = sol.zero[3]^2
        l_P2 = sol.zero[4]^2
        l_G2 = sol.zero[5]^2
        D_E2 = sol.zero[6]^2
        D_I = data_gen(l_G; par, type)
        D_I2 = data_gen(l_G2; par, type)
        D = bundle(D_I, D_E; par, type)
        D2 = bundle(D_I2, D_E2; par, type)
        l = labor_demand(l_P, D; par, type, sec_period=:no)
        l2 = labor_demand(l_P2, D2; par, type, sec_period=:yes)
    else
        error("I should never be here")
    end
end

function set_a_Y(compl_prod, σ, ζ)
    if compl_prod == :yes
        α_Y = (σ * ζ - σ + 1) / (σ * ζ) # auxiliary paramter
    elseif compl_prod == :no
        α_Y = 0
    else
        error("Field compl_prod needs to be either :no or :yes")
    end
    return α_Y
end

"""
    Return the correct gamma
"""
function gamma_of_type(par::Pars_v3; type::Symbol)
    if type in [:BS, :LS]
        par.γ_S
    elseif type in [:BU, :LU]
        par.γ_U
    else
        error("type needs to be either :BS, :LS, :BU or :LU")
    end
end

function K_last_of_type(par::Pars_v3; type::Symbol)
    if type == :LU # little, unsophisticated
        par.K_last_LU
    elseif type == :LS # little, sophisticated
        par.K_last_LS
    elseif type == :BU # big, unsophisticated
        par.K_last_BU
    elseif type == :BS # big, sophisticated
        par.K_last_BS
    else
        error("Firm type needs to be either LU, LS, BU or BS")
    end
end

function D_last_of_type(par::Pars_v3; type::Symbol)
    if type == :LU # little, unsophisticated
        par.D_last_LU
    elseif type == :LS # little, sophisticated
        par.D_last_LS
    elseif type == :BU # big, unsophisticated
        par.D_last_BU
    elseif type == :BS # big, sophisticated
        par.D_last_BS
    else
        error("Firm type needs to be either LU, LS, BU or BS")
    end
end