#! This is the package in which we put all our dependencies, functions, types etc.

module InnovationPrivacy

using Parameters, Setfield

using QuantEcon: qnwnorm, gridmake
using StatsFuns: normcdf
using Roots: find_zero
using Distributions: cdf, Pareto
using StatsBase: ecdf


include("Types.jl")

include("Functions.jl")
include("Functions_Replication.jl")
include("TestFuncs.jl")

export  Pars, ParsGrowth, BPars, BParsGrowth,
        @unpack_Pars, @unpack_ParsGrowth, @unpack_BPars, @unpack_BParsGrowth,
        firm_profit,
        capital_demand,
        test_func,
        Eπ, Er, EA,
        extractK,
        signal,
        Tau_of_K_old,
        solve_model,
        K_solve,
        Transition_Probability,
        find_limit,
        find_V_n_max,
        find_V_n,
        find_all_n,
        find_V_zero,
        find_u_star,
        u_to_z,
        z_to_u,
        Firm_Distribution,
        iterate_population,
        Eπ_Distribution,
        Er_Distribution,
        labor_costs,
        μ_n,
        Transition_Matrix2,
        Equilibrium_Distribution
end # module
