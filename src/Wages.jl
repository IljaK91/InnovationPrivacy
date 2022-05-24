"""
    Closed form expression for the wage for labor that produces intermediate goods
"""
function wage_int_goods(Π; par)
    Π_BS = Π[1] # firm revenue
    Π_LS = Π[2] # firm revenue
    Π_LU = Π[3] # firm revenue
    Π_BU = Π[4] # firm revenue

    w_BS = weight_of_type(:BS; par)
    w_LS = weight_of_type(:LS; par)
    w_LU = weight_of_type(:LU; par)
    w_BU = weight_of_type(:BU; par)

    return par.ζ * par.α_L_hat * (w_BS * Π_BS + w_LS * Π_LS + w_LU * Π_LU + w_BU * Π_BU)
end

function wage_gen_data(Π, D, l_G, D_I; par)
    Π_BS = Π[1] # firm revenue
    Π_LS = Π[2] # firm revenue
    Π_LU = Π[3] # firm revenue
    Π_BU = Π[4] # firm revenue

    D_BS = D[1] # data bundle
    D_LS = D[2] # data bundle
    D_LU = D[3] # data bundle
    D_BU = D[4] # data bundle

    D_I_BS = D_I[1] # data bundle
    D_I_LS = D_I[2] # data bundle
    D_I_LU = D_I[3] # data bundle
    D_I_BU = D_I[4] # data bundle

    l_G_BS = l_G[1] # data bundle
    l_G_LS = l_G[2] # data bundle
    l_G_LU = l_G[3] # data bundle
    l_G_BU = l_G[4] # data bundle

    w_BS = weight_of_type(:BS; par) # weights
    w_LS = weight_of_type(:LS; par) # weights
    w_LU = weight_of_type(:LU; par) # weights
    w_BU = weight_of_type(:BU; par) # weights

    D_G_BS = data_gen(l_G_BS; par, type = :BS)
    D_G_LS = data_gen(l_G_LS; par, type = :LS)
    D_G_LU = data_gen(l_G_LU; par, type = :LU)
    D_G_BU = data_gen(l_G_BU; par, type = :BU)

    c_BS = D_G_BS * Π_BS / D_BS * (D_BS / D_I_BS)^(1 / par.ε) # Put some things together
    c_LS = D_G_LS * Π_LS / D_LS * (D_LS / D_I_LS)^(1 / par.ε) # Put some things together
    c_LU = D_G_LU * Π_LU / D_LU * (D_LU / D_I_LU)^(1 / par.ε) # Put some things together
    c_BU = D_G_BU * Π_BU / D_BU * (D_BU / D_I_BU)^(1 / par.ε) # Put some things together

    return par.ζ * par.α_K_hat * (1 - par.ϕ) * ((1 - par.γ_S) * w_BS * c_BS + (1 - par.γ_S) * w_LS * c_LS + (1 - par.γ_U) * w_LU * c_LU + (1 - par.γ_U) * w_BU * c_BU)
end

"""
    Wage for data processing.
"""
function wage_gen_proc(Π, D; par)
    Π_BS = Π[1] # firm revenue
    Π_LS = Π[2] # firm revenue
    Π_LU = Π[3] # firm revenue
    Π_BU = Π[4] # firm revenue

    w_BS = weight_of_type(:BS; par) # weights
    w_LS = weight_of_type(:LS; par) # weights
    w_LU = weight_of_type(:LU; par) # weights
    w_BU = weight_of_type(:BU; par) # weights

    return par.ζ * par.α_K_hat * (par.γ_S * w_BS * Π_BS + par.γ_S * w_LS * Π_LS + par.γ_U * w_LU * Π_LU + par.γ_U * w_BU * Π_BU)
end