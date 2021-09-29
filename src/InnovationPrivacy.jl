#! This is the package in which we put all our dependencies, functions, types etc.

module InnovationPrivacy

using Parameters

using QuantEcon: qnwnorm

include("Types.jl")

include("Functions.jl")

export  Pars,
        BPars,
        firm_profit,
        capital_demand,
        @unpack_Pars,
        @unpack_BPars
end # module
