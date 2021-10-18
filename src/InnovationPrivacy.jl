#! This is the package in which we put all our dependencies, functions, types etc.

module InnovationPrivacy

using Parameters

using QuantEcon: qnwnorm, gridmake

include("Types.jl")

include("Functions.jl")
include("Functions_Replication.jl")

export  Pars, ParsGrowth, BPars, BParsGrowth,
        @unpack_Pars, @unpack_ParsGrowth, @unpack_BPars, @unpack_BParsGrowth,
        firm_profit,
        capital_demand,
        test_func,
        Eπ,
        extractK,
        signal,
        Tau_of_K_old,
        solve_model,
        K_solve
end # module
