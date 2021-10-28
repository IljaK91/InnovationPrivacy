using Test, InnovationPrivacy, StatsBase

@testset "Check Transition Probabilities sum up to one" begin
    for N in 5:1:15, n in 0:1:10
       for i in 1:N
           @test sum(Transition_Probability(n, ParsGrowth(N = N), normalize = :no)[i, :]) ≈ 1.0
       end
    end
end

@testset "Check Bayesian Identities" begin
    for σ²_θ in 0.1:0.2:1.5, σ²_ϵ in 0.1:0.2:1.5, n in 0:1:10
        bpar = BPars(n, ParsGrowth(σ²_θ = σ²_θ, σ²_ϵ = σ²_ϵ))
        @test bpar.w_an + bpar.w_pn ≈ 1.0
    end
end

@testset "Check Expected Profit Identity" begin
    for σ²_θ in 0.1:0.5:1.6, σ²_ϵ in 0.1:0.5:1.6, n in 0:2:10, a_bar in -2:0.2:2
        par  = ParsGrowth(σ²_θ = σ²_θ, σ²_ϵ = σ²_ϵ)
        bpar = BPars(n, par)

        e_profits = Eπ(a_bar; par, bpar)
        e_revenue = Er(a_bar; par, bpar)
        wage_bill = labor_costs(a_bar; par, bpar)


        @test e_profits - (e_revenue - wage_bill - par.w*par.f) ≈ 0.
    end
end

@testset "Check Transition Matrix" begin
    for σ²_θ in 0.1:0.5:1.1, σ²_ϵ in 0.1:0.5:1.1, n in 0:5:10
    
        par      = ParsGrowth(σ²_θ = σ²_θ, σ²_ϵ = σ²_ϵ)
        bpar     = BPars(n, par)
        
        Π_test   = Transition_Matrix2(par, bpar)
        Π_direct = Transition_Probability(par, bpar)

        @test sum((Π_direct - Π_test).^2) < 1e-4
    end
end