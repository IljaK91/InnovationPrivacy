#! This is the package in which we put all our dependencies, functions, types etc.

module InnovationPrivacy

using Parameters, Setfield

using QuantEcon: qnwnorm, gridmake
using StatsFuns: normcdf
using Roots: find_zero
using Distributions: cdf, Pareto


include("Types.jl")

include("Functions.jl")
include("Functions_Replication.jl")

export  Pars, ParsGrowth, BPars, BParsGrowth,
        @unpack_Pars, @unpack_ParsGrowth, @unpack_BPars, @unpack_BParsGrowth,
        firm_profit,
        capital_demand,
        test_func,
        Eπ, Er,
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
        find_V_0,
        find_u_star,
        u_to_z,
        z_to_u,
        Firm_Distribution,
        iterate_population,
        Eπ_Distribution,
        Er_Distribution
end # module
