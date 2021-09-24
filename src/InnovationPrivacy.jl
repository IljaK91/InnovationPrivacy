#! This is the package in which we put all our dependencies, functions, types etc.

module InnovationPrivacy

using Parameters

include("Functions.jl")
include("Types.jl")

export  Pars,
        BPars,
        firm_profit
end # module
