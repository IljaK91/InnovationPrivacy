"""
    Constructor for the Bayesian Parameters (weights on private signal and prior)
"""
function BPars(n, par::ParsGrowth)
    @unpack_ParsGrowth par
    V_n  = 1/(σ²_θ^-1 + n*σ²_ϵ^-1)
    w_pn = σ²_θ^-1/(σ²_θ^-1 + n*σ²_ϵ^-1)
    w_an = n*σ²_ϵ^-1/(σ²_θ^-1 + n*σ²_ϵ^-1)
    BParsGrowth(V_n  = V_n,
                w_pn = w_pn,
                w_an = w_an)
end
BPars(n; par::ParsGrowth) = BParsGrowth(n, par)