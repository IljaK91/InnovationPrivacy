"""
    First order condition for data generation
"""
function FOC_data_gen(l_G, D, l, D_I; par::Pars_v3, type::Symbol)
    @unpack_Pars_v3 par
    @assert type in [:LU, :LS, :BU, :BS]

    D_G = data_gen(l_G; par, type)
    Π = firm_revenue(D, l; par, type)
    α_K_hat = get_alpha_K_type(type; par)

    α_K_hat * (1 - ϕ) * D_G / l_G * Π / D * dD_dDI(D, D_I; par) - w_G
end

function FOC_data_gen2(l_G, l_P, D, l, D_I, l_P_next, D_next, l_next; par::Pars_v3, type::Symbol)
    @unpack_Pars_v3 par
    @assert type in [:LU, :LS, :BU, :BS]

    D_G = data_gen(l_G; par, type)
    Π = firm_revenue(l_P, D, l; par, type)
    Π_next = firm_revenue(l_P_next, D_next, l_next; par, type)
    γ = gamma_of_type(par; type)

    α_K_hat * (1 - γ) * (1 - ϕ) * D_G / l_G * Π / D * dD_dDI(D, D_I; par) +
    β * (1 - δ_K) * α_K_hat * (1 - γ) * (1 - ϕ) * D_G / l_G * Π_next / D_next * dD_dDI(D, D_I; par) -
    w_G
end

"""
    First order condition for intermediate good labor
"""
function FOC_data_proc(l_P, D, l; par::Pars_v3, type::Symbol, sec_period::Symbol = :no)
    @unpack_Pars_v3 par
    @assert type in [:LU, :LS, :BU, :BS]

    γ = gamma_of_type(par; type)
    Π = firm_revenue(l_P, D, l; par, type)
    if sec_period == :no
        α_K_hat * γ * Π / l_P - w_P
    elseif sec_period == :yes
        α_K_hat * γ * Π / l_P - w_P2
    end
end

function FOC_data_proc2(l_P, D, l, l_P_next, D_next, l_next; par::Pars_v3, type::Symbol)
    @unpack_Pars_v3 par
    @assert type in [:LU, :LS, :BU, :BS]

    γ = gamma_of_type(par; type)
    Π = firm_revenue(l_P, D, l; par, type)
    Π_next = firm_revenue(l_P_next, D_next, l_next; par, type)
    K_next = knowledge(l_P_next, D_next; par, type, sec_period = :yes)
    K = knowledge(l_P, D; par, type, sec_period = :no)
    α_K_hat * γ * Π / l_P + β * (1 - δ_K) * γ * α_K_hat * Π_next / K_next * K / l_P - w_P
end

"""
    First order condition for data sharing
"""
function FOC_data_sharing(D, DI, l; par::Pars_v3, type::Symbol)
    @unpack_Pars_v3 par
    @assert type in [:LU, :LS, :BU, :BS]

    Π = firm_revenue(D,l; par, type)
    α_K_hat = get_alpha_K_type(type;par)
    p_D - α_K_hat * Π / D * dD_dDS(D, DI; par)
end

function FOC_data_sharing2(l_P, D, DI, l, l_P_next, D_next, l_next; par::Pars_v3, type::Symbol)
    @unpack_Pars_v3 par
    @assert type in [:LU, :LS, :BU, :BS]

    Π = firm_revenue(l_P, D, l; par, type)
    Π_next = firm_revenue(l_P_next, D_next, l_next; par, type)
    γ = gamma_of_type(par; type)

    p_D - α_K_hat * (1 - γ) * Π / D * dD_dDS(D, DI; par) - β * (1 - δ_K) * α_K_hat * (1 - γ) * Π_next / D_next * dD_dDS(D, DI; par)
end

"""
    First order condition for data buying
"""
function FOC_data_buying(D, l, DE; par::Pars_v3, type::Symbol)
    @unpack_Pars_v3 par
    @assert type in [:LU, :LS, :BU, :BS]

    Π = firm_revenue(D, l; par, type)
    α_K_hat = get_alpha_K_type(type; par)
    α_K_hat * Π / D * dD_dDE(D, DE; par) - p_D / (1 - τ)
end

"""
    First order condition for data buying, 2 period case.
"""
function FOC_data_buying2(l_P, D, l, DE, l_P_next, D_next, l_next; par::Pars_v3, type::Symbol)
    @unpack_Pars_v3 par
    @assert type in [:LU, :LS, :BU, :BS]

    γ = gamma_of_type(par; type)
    Π = firm_revenue(l_P, D, l; par, type)
    Π_next = firm_revenue(l_P_next, D_next, l_next; par, type)
    α_K_hat * (1 - γ) * Π / D * dD_dDE(D, DE; par) + β * (1 - δ_D) * α_K_hat * (1 - γ) * Π_next / D_next * dD_dDE(D, DE; par) - p_D / (1 - τ)
end

dD_dDI(D, DI; par::Pars_v3) = (D / DI)^(1 / par.ε)
dD_dDS(D, DI; par::Pars_v3) = (1 - par.ν) * (D / DI)^(1 / par.ε)
dD_dDE(D, DE; par::Pars_v3) = par.ξ * (D / DE)^(1 / par.ε)


"""
    First order condition for intermediate good labor
"""
function FOC_int_labor(l_P, D, l; par::Pars_v3, type::Symbol, sec_period::Symbol=:no)
    @unpack_Pars_v3 par
    @assert type in [:LU, :LS, :BU, :BS]

    Π = firm_revenue(l_P, D, l; par, type)
    if sec_period == :no
        α_L_hat * Π / l - w
    elseif sec_period == :yes
        α_L_hat * Π / l - w2
    end
    #K = knowledge(l_P, D; par, type)
    #Y^α_Y * α_L_hat * K^α_K_hat * l^(α_L_hat - 1) - w
end

"""
    The value of the data bundle using both internal and external data.
"""
function bundle(D_I, D_E; par::Pars_v3, type::Symbol)
    @unpack_Pars_v3 par
    (D_I^((ε - 1) / ε) + ξ * D_E^((ε - 1) / ε))^(ε / (ε - 1))
end

data_multiplier(D_G, D_S; par::Pars_v3) = 1 - (par.τ - par.ν) * D_S / D_G
data_multiplier2(D_G, D_S; par::Pars_v3) = 1 - (par.τ2 - par.ν2) * D_S / D_G
data_multiplier_bundle(D, D_G) = D / D_G


function knowledge(l_P, D; par::Pars_v3, type::Symbol, sec_period::Symbol = :no)
    A_P = A_P_of_type(type; par)
    γ = gamma_of_type(par; type)
    K_last = K_last_of_type(par; type)
    if sec_period == :yes
        A_P * l_P^γ * D^(1 - γ) + (1 - par.δ_K) * K_last
    elseif sec_period == :no
        A_P * l_P^γ * D^(1 - γ)
    else
        error("second_period needs to be either :yes or :no")
    end
end

knowledge(l_P, D, type::Symbol; par::Pars_v3, sec_period::Symbol = :no) = knowledge(l_P, D; par, type, sec_period)

function data_gen(l_G; par::Pars_v3, type::Symbol)
    A_G = A_G_of_type(type; par)
    A_G * l_G^(1 - par.ϕ)
end
data_gen(l_G, type::Symbol; par::Pars_v3) = data_gen(l_G; par, type)
data_intern(D_G, D_S; par::Pars_v3) = D_G - (1 - par.ν) * D_S


function sold_data(D_I, l_G; par::Pars_v3, type::Symbol)
    D_G = data_gen(l_G; par, type)

    return (D_G - D_I) / (1 - par.ν)
end

sold_data(D_I, l_G, type; par::Pars_v3) = sold_data(D_I, l_G; par, type)

"""
    output = Y_i = K^α*l^(1-α)
"""
function firm_output(l_P, D, l; par::Pars_v3, type::Symbol)
    K = knowledge(l_P, D; par, type)
    K^par.α * l^(1 - par.α)
end
firm_output(l_P, D, l, type; par::Pars_v3) = firm_output(l_P, D, l; par, type)

"""
    revenue = p_i * Y_i
"""
function firm_revenue(D, l; par::Pars_v3, type::Symbol, sec_period::Symbol = :no)
    α_K_hat = get_alpha_K_type(type; par)
    α_L_hat = get_alpha_L_type(type; par)
    par.ζ * par.Y^par.α_Y * D^α_K_hat * l^α_L_hat
end

firm_revenue(D, l, type; par::Pars_v3, sec_period::Symbol = :no) = firm_revenue(D, l; par, type, sec_period)

#! this does not work anymore, need to see where it throws an error...
# function firm_revenue(K, l; par::Pars_v3, sec_period::Symbol = :no)
#     if sec_period == :no
#         par.ζ * par.Y^par.α_Y * K^par.α_K_hat * l^par.α_L_hat
#     elseif sec_period == :yes
#         par.ζ * par.Y2^par.α_Y * K^par.α_K_hat * l^par.α_L_hat
#     end
# end 

"""
    Analytical solution for labor demand
"""
function labor_demand(D; par::Pars_v3, type::Symbol)
    @unpack_Pars_v3 par
    α_K_hat = get_alpha_K_type(type; par)
    α_L_hat = get_alpha_L_type(type; par)
    (ζ * α_L_hat * Y^(α_Y) * D^(α_K_hat) / w)^(1 / (1 - α_L_hat))
end
labor_demand(D, type; par::Pars_v3, sec_period::Symbol = :no) = labor_demand(D; par, type, sec_period)

"""
    Analytical Solution when D_S > 0.
"""
function labor_gen(par::Pars_v3; type::Symbol, sec_period::Symbol=:no)
    @unpack_Pars_v3 par
    A_G = A_G_of_type(type; par)
    if sec_period == :no
        ((1 - ϕ) * A_G / (1 - ν) * p_D / w_G)^(1 / ϕ)
    elseif sec_period == :yes
        ((1 - ϕ) * A_G / (1 - ν) * p_D2 / w_G2)^(1 / ϕ)
    end
end


"""
    Firm profits
"""
function firm_profits(l, l_P, l_G, D_S, D_E; par::Pars_v3, type::Symbol)
    D_G = data_gen(l_G; par, type)
    D_I = data_intern(D_G, D_S; par)
    D = bundle(D_I, D_E; par, type)
    K = knowledge(l_P, D; par, type)
    Π = firm_revenue(K, l; par)
    Π - par.w * l - par.w_G * l_G - par.w_P * l_P + par.p_D * D_S - par.p_D / (1 - par.τ) * D_E
end
firm_profits(l, l_P, l_G, D_S, D_E, type::Symbol; par::Pars_v3) = firm_profits(l, l_P, l_G, D_S, D_E; par, type)