"""
    Simpler and more direct way to compute the transition matrix. Less precise.
"""
function Transition_Matrix2(par::ParsGrowth, bpar::BParsGrowth; T = 1_000_000)
    @unpack_ParsGrowth par
    @unpack_BParsGrowth bpar
    Π     = zeros(N, N)
    for i in 1:N
        a_bar = a_grid[i]
        μ     = μ_n(a_bar, par, bpar)

        a_next     = μ .+ √V_n*randn(T) .+ √σ²_ϵ.*randn(T)
        a_bar_next = n./(n+1).*a_bar .+ 1 ./(n+1).*a_next

        f          = ecdf(a_bar_next)

        for j in 1:N
            if j == 1
                Π[i,j] = f(a_thresholds[1])
            elseif j == N
                Π[i,j] = 1 - f(a_thresholds[j-1])
            else
                Π[i,j] = f(a_thresholds[j]) - f(a_thresholds[j-1])
            end
        end 
    end
    return Π
end
