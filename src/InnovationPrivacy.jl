#! This is the package in which we put all our dependencies, functions, types etc.
#! In desparate need of cleaning up!
module InnovationPrivacy
using Suppressor
@suppress using Parameters, Setfield, Optim, NLsolve, DataFrames, Interpolations

@suppress using QuantEcon: qnwnorm, gridmake
#using StatsFuns: normcdf
@suppress using Roots: find_zero
@suppress using Distributions: cdf, Pareto
@suppress using StatsBase: ecdf
@suppress using Statistics: quantile

@suppress using Plots
#pgfplotsx()

include("Types.jl")
include("Functions_olg.jl")
include("Functions_v3.jl")
include("FOCS_v3.jl")
#include("Retired.jl")
#include("PlotFuncs.jl")
include("Auxiliary_Funcs.jl")
include("Wages.jl")
include("Residuals.jl")
include("Functions_plot_comp_statics.jl")

#include("Functions.jl")
#include("FunctionsDeprInfo.jl")
#include("Functions_Replication.jl")
#include("TestFuncs.jl")
#include("Functions_simple.jl")
#include("Functions_simple_solving.jl")
#include("Functions_simple_model.jl")
#include("Functions_bundle.jl")

export Pars, ParsGrowth, Pars_Bundle, BPars, BParsGrowth, Pars_Simple, Pars_OLG, Pars_v3, residuals, solve_firm_problem, FOC_data_sharing, bundle, dD_dDI, dD_dDE, FOC_data_gen, FOC_data_buying, FOC_int_labor, data_gen, get_solution, weight_of_type, sold_data, firm_output, aggregate_output, residuals_SS, find_steady_state, get_ss_solution, residuals_SS_single, find_p_D, find_steady_state_iterative, comp_statics, residuals_p_D, firm_revenue, sol_f_problem, find_steady_state_PE, firm_profits, solve_firm_problem_unc, solve_firm_problem_con, solve_firm_problem_unc2, solve_firm_problem_nl, sol_f_problem_nl, find_value_func, Firm_sol, Save_Firm_Solution, solve_firm_problem_unc_nl2, solve_firm_problem_unc_nl, get_solution_nl2, A_G_of_type, solve_firm_problem_con_nl,
    @unpack_Pars, @unpack_ParsGrowth, @unpack_BPars, @unpack_BParsGrowth, @unpack_Pars_OLG, @unpack_Pars_Simple, @unpack_Pars_Bundle, @unpack_Pars_v3, @unpack_Firm_sol,
    firm_profit, MPD, get_alpha_K_type, plot_comp_static_firm, plot_comp_static_firm_export,
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
find_all_V_n,
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
Equilibrium_Distribution,
find_V_0,
mass_and_entry,
mean_preserving_noise,
CalcAllUncertainty,
BPars_dep,
Eπ_dep,
Transition_Probability_dep,
ParsV2,
find_price_info,
solve_agg_C,
find_agg_C,
wage_rate,
find_wage_rate,
capital_prod,
solve_agg_C_iterative,
capital_of_tau,
profit,
k_I_prod_A,
k_I_prod_B,
solve_firm_problems,
k_I_prod_A,
k_I_prod_B,
loss_A,
loss_B,
loss_A_con,
solve_agg_Y_iterative_unc,
profit_unc,
find_mass,
solve_U_problem
end # module
