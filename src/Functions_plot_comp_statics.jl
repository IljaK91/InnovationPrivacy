
function plot_comp_static_firm(x, y; xlabel = "", ylabel = "", Filename::String, Folder::String)
    y_LU = y[1]
    y_BU = y[2]
    y_LS = y[3]
    y_BS = y[4]

    plot(x, y_LU ./ y_LU[1], xlabel = xlabel, ylabel = ylabel, label = "Small, unsophisticated", color = :black, ls = :dash)
    plot!(x, y_BU ./ y_BU[1], color = :black, label = "Big, unsophisticated")
    plot!(x, y_LS ./ y_LS[1], color = :red, ls = :dash, label = "Small, sophisticated")
    plot!(x, y_BS ./ y_BS[1], color = :red, label = "Big, sophisticated")

    savefig("../Graphs/v3/" * Folder * "/" * Filename * "_normalized.pdf")
    savefig("../Graphs/v3/" * Folder * "/" * Filename * "_normalized.png")

    plot(x, y_LU, color = :black, ls = :dash, xlabel = xlabel, ylabel = ylabel, label = "Small, unsophisticated")
    plot!(x, y_BU, color = :black, label = "Big, unsophisticated")
    plot!(x, y_LS, color = :red, ls = :dash, label = "Small, sophisticated")
    plot!(x, y_BS, color = :red, label = "Big, sophisticated")

    savefig("../Graphs/v3/" * Folder * "/" * Filename * ".pdf")
    savefig("../Graphs/v3/" * Folder * "/" * Filename * ".png")
end
function plot_comp_static_firm_export(x, y; xlabel = "", ylabel = "", Filename::String, Folder::String, legend = :no, norm = :no)
    if norm == :yes
        y_LU = y[1]./y[1][1]
        y_BU = y[2]./y[2][1]
        y_LS = y[3]./y[3][1]
        y_BS = y[4]./y[4][1]
    elseif norm == :no
        y_LU = y[1]
        y_BU = y[2]
        y_LS = y[3]
        y_BS = y[4]
    else
        error("Keyword norm needs to be either :yes or :no")
    end

    if legend == :no
        p = plot(x, y_LU, color = :black, ls = :dash, xlabel = xlabel, ylabel = ylabel, label = "")
        plot!(x, y_BU, color = :black, label = "")
        plot!(x, y_LS, color = :red, ls = :dash, label = "")
        plot!(x, y_BS, color = :red, label = "")
    elseif legend == :yes
        p = plot(x, y_LU, color = :black, ls = :dash, xlabel = xlabel, ylabel = ylabel, label = "Small, unsophisticated")
        plot!(x, y_BU, color = :black, label = "Big, unsophisticated")
        plot!(x, y_LS, color = :red, ls = :dash, label = "Small, sophisticated")
        plot!(x, y_BS, color = :red, label = "Big, sophisticated")
    end
    return p
end
function plot_comp_static_firm1d(x, y; xlabel="", ylabel="", Filename::String, Folder::String)
    y_L = y[1]
    y_B = y[2]

    plot(x, y_L ./ y_L[1], xlabel=xlabel, ylabel=ylabel, label="low data collection ability", color=:black)
    plot!(x, y_B ./ y_B[1], color=:red, label="high data collection ability")

    savefig("../Graphs/v3/" * Folder * "/" * Filename * "_normalized_1d.pdf")
    savefig("../Graphs/v3/" * Folder * "/" * Filename * "_normalized_1d.png")

    plot(x, y_L, xlabel=xlabel, ylabel=ylabel, label="low data collection ability", color=:black)
    plot!(x, y_B, color=:red, label="high data collection ability")

    savefig("../Graphs/v3/" * Folder * "/" * Filename * "_1d.pdf")
    savefig("../Graphs/v3/" * Folder * "/" * Filename * "_1d.png")
end
function plot_comp_static_firm_export1d(x, y; xlabel = "", ylabel = "", Filename::String, Folder::String, legend = :no, norm = :no)
    if norm == :yes
        y_L = y[1]./y[1][1]
        y_B = y[2]./y[2][1]
    elseif norm == :no
        y_L = y[1]
        y_B = y[2]
    else
        error("Keyword norm needs to be either :yes or :no")
    end

    if legend == :no
        p = plot(x, y_L, color=:black, xlabel=xlabel, ylabel=ylabel, label="")
        plot!(x, y_B, color=:red, label="")
    elseif legend == :yes
        p = plot(x, y_L, color=:black, xlabel=xlabel, ylabel=ylabel, label="low data collection ability")
        plot!(x, y_B, color=:red, label="high data collection ability")
    end
    return p
end