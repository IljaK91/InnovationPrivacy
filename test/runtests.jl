using Test, InnovationPrivacy

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