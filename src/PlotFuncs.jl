
function plot_comp_static_firm(x, y; xlabel = "", ylabel = "", Filename::String, Folder::String)
    y_LU = y[1]
    y_BU = y[2]
    y_LS = y[3]
    y_BS = y[4]

    plot(x, y_LU ./ y_LU[1], xlabel = xlabel, ylabel = ylabel, label = "Small, unsophisticated", color = :black, ls = :dashdot)
    plot!(x, y_BU ./ y_BU[1], color = :blue, ls = :dot, label = "Big, unsophisticated")
    plot!(x, y_LS ./ y_LS[1], color = :green, ls = :dash, label = "Small, sophisticated")
    plot!(x, y_BS ./ y_BS[1], color = :red, label = "Big, sophisticated")

    savefig("../Graphs/v3/" *Folder* "/" * Filename * "_normalized.pdf")
    savefig("../Graphs/v3/" *Folder* "/" * Filename * "_normalized.png")

    plot(x, y_LU, color = :black, ls = :dashdot, xlabel = xlabel, ylabel = ylabel, label = "Small, unsophisticated")
    plot!(x, y_BU, color = :blue, ls = :dot, label = "Big, unsophisticated")
    plot!(x, y_LS, color = :green, ls = :dash, label = "Small, sophisticated")
    plot!(x, y_BS, color = :red, label = "Big, sophisticated")

    savefig("../Graphs/v3/" *Folder* "/" * Filename * ".pdf")
    savefig("../Graphs/v3/" *Folder* "/" * Filename * ".png")
end