basis_order(basis) = Dict('d' => 2, 't' => 3, 'q' => 4, '5' => 5, '6' => 6)[lowercase(basis[end-1])]

function maglegend(mags, bases, num_electrons; fraction = (m, n) -> m / n)
    markers = [
        [
            [
                MarkerElement(
                    marker = magmarker(
                        n;
                        fraction = fraction(mag, num_electrons),
                        radius = 0.5,
                    ),
                    color = mag == 0 ? :transparent : :black,
                ),
                MarkerElement(
                    marker = magmarker(n; fraction = 1, radius = 0.5),
                    strokewidth = 1,
                    color = :transparent,
                    strokecolor = :black,
                ),
            ] for mag in mags
        ] for n in basis_order.(bases)
    ]
    maglabel(m) = m == 0 ? L"m_\mathrm{S} = 0" : L"m_\mathrm{S} =\!±\!%$(round(Int,m/2))"
    labels = [
        [i == length(bases) ? maglabel(mag) : " " for mag in mags] for i in eachindex(bases)
    ]

    titles = bases
    return markers, labels, titles
end


function plot_hring(ax, df; colorrange, colormap)
    num_electrons = df.num_electrons[1]

    df = @by(df, [:dimerization_angle, :λ, :basis], begin
        :mag = :mag[argmin(:energy_hf)]
        :energy_hf = minimum(:energy_hf)
    end)

    bases = map(x -> uppercase(x[end-1:end]), unique(df.basis))
    mags = sort(unique(df.mag))

    sort!(df, :λ, rev = true)
    for dfg in groupby(df, [:λ, :basis])
        corners = basis_corners(dfg.basis[1])

        markerfull = magmarker(corners; fraction = 1, radius = 0.49)

        scatterlines!(
            ax,
            dfg.dimerization_angle,
            dfg.energy_hf;
            markercolor = :white,
            color = dfg.λ[1],
            marker = markerfull,
            strokecolor = colormap[-(dfg.λ[1] - colorrange[1])/-(colorrange...)],
            strokewidth = 1,
            colorrange,
            colormap,
        )
        gf = @rsubset(dfg, :mag != 0)
        if size(gf, 1) != 0
            marker = [
                magmarker(corners; fraction = mag / num_electrons, radius = 0.5) for
                mag in gf.mag
            ]
            scatter!(
                ax,
                gf.dimerization_angle,
                gf.energy_hf;
                color = gf.λ[1],
                colorrange,
                colormap,
                marker = marker,
            )
        end
    end
    @show bases

    markers, labels, titles = maglegend(mags, bases, num_electrons)
    axislegend(
        ax,
        markers,
        labels,
        titles,
        titlegap = 3,
        padding = 0,
        orientation = :horizontal,
        nbanks = length(mags),
        position = (:right, :bottom),
        groupgap = 0,
        margin = (3, 3, 3, 3),
    )
end

function fig_hrings()
    df = @rsubset(load_data("hrings_hf"), :λ < 0.02)
    @rtransform!(df, :dimerization_angle = 360 * :dimerization / :num_electrons)
    fig = Figure(size = (240 * 3, 280), figure_padding = 0)
    axs =
        [Axis(fig[1, i], xlabel = L"Dimerization $δ$", ylabel = "Energy (Ha)") for i = 1:3]
    dfs = [@rsubset(df, :num_electrons == num_electrons && :basis == (num_electrons == 8 ? "aug-cc-pvtz" : "aug-cc-pvqz")) for num_electrons in [4, 6, 8]]

    for (i, ax) in enumerate(axs)
        if i > 1
            hideydecorations!(ax, label = true, ticks = false, ticklabels = false)
        end
        n = dfs[i].num_electrons[1]
        xlims!(ax, -4 / n, 180 / n^1.3 * 4^0.3)
    end

    couplings = sort(unique(df.λ))
    δλ = couplings[2] - couplings[1]
    coupling_boundaries = vcat(couplings .- δλ / 2, [couplings[end] + δλ / 2])
    colorrange = extrema(coupling_boundaries)
    colormap = cgrad(
        [Makie.wong_colors()[1], Makie.wong_colors()[2]],
        length(couplings),
        categorical = true,
    )

    plot_hring.(axs, dfs; colormap, colorrange)
    Colorbar(fig[1, 4]; label = L"$λ$ (a. u.)", colorrange, colormap, ticks = couplings)

    ylims!(axs[1], -1.91, -1.601)
    ylims!(axs[2], -3.07, -2.69)
    ylims!(axs[3], -4.0, -3.1)

    axs[3].yticks = [-3.4, -3.7, -4.0, -4.3]

    xlim = [(-1 / 2, 45), (-1 / 2, 24), (-1 / 2, 22.5)]

    s = 70
    x, y = 10, 140
    scatter!(
        axs[1],
        [x + s / 2],
        [y + s / 2],
        marker = :circle,
        color = :transparent,
        markersize = 1.14s,
        strokecolor = :grey,
        strokewidth = 1,
        space = :pixel,
    )
    scatter!(
        axs[1],
        [x + s / 2],
        [y + s / 2],
        marker = '⊗',
        color = colormap[3],
        markersize = 12,
        space = :pixel,
    )
    scatter!(
        axs[1],
        [x + s / 2 - 5],
        [y + s / 2 + 10],
        marker = '𝝺',
        color = colormap[3],
        markersize = 14,
        space = :pixel,
    )

    for (i, (df, ax, lim)) in enumerate(zip(dfs, axs, xlim))
        xlims!(ax, lim)
        ax.xticks = ([0, 10, 20, 30, 40], ["0°", "10°", "20°", "30°", "40°"])
        text!(
            ax,
            0,
            1,
            text = rich(
                rich("$(Char(codepoint('a')-1+i))", font = :bold),
                " H",
                subscript(string(df.num_electrons[1])),
                ", R = $(df.R[1]) Å",
            ),
            align = (:left, :top),
            offset = (4, -3),
            space = :relative,
        )

        img = load("$(path())/assets/h$(df.num_electrons[1]).png")
        image!(ax, (x, x + s), (y, y + s), img, space = :pixel)
    end

    f = 0.805
    t = (-30 .- range(0, 25, 40))
    lines!(
        axs[1],
        x + s / 2 .+ f * s / 2 .* cosd.(t),
        y + s / 2 .+ f * s / 2 .* sind.(t),
        space = :pixel,
        color = :black,
        linewidth = 2,
    )
    text!(
        axs[1],
        x + s / 2 .+ 1.4f * s / 2 * cosd(-55),
        y + s / 2 .+ 1.4 * f * s / 2 * sind(-55),
        space = :pixel,
        text = "δ",
        fontsize = 12,
    )
    angle = -90
    d = f * s / 2 * [cosd(angle), sind(angle)]
    lines!(
        axs[1],
        [x + s / 2 + 0.17d[1], x + s / 2 + d[1]],
        [y + s / 2 + 0.17d[2], y + s / 2 + d[2]],
        space = :pixel,
        color = :black,
        linewidth = 2,
    )
    text!(
        axs[1],
        x + s / 2 - 2,
        y + s / 2 - f * s / 4,
        space = :pixel,
        text = "R",
        align = (:right, :center),
        fontsize = 12,
    )

    colgap!(fig.layout, 6)


    fig
end

function plot_hrings_qmc(ax, df)

    R = df.R[1]
    df = @rsubset(dropmissing(df, :energy_qmc), !isnothing(:energy_qmc))

    df = @by(
        df,
        [:R, :dimerization, :λ, :num_electrons, :basis],
        begin
            :mag = :mag[argmin(:energy_qmc)]
            :energy_qmc = minimum(:energy_qmc) ± :energy_qmc_error[argmin(:energy_qmc)]
            :energy_hf = :energy_hf[argmin(:energy_qmc)]
        end
    )

    num_electrons = Int(df.num_electrons[1])

    function fraction(mag, num_electrons)
        f = mag / num_electrons
        if num_electrons == 8
            f = f^0.6
        end
        return f
    end

    sort!(df, :dimerization)
    for gbas in groupby(df, :basis)
        λs = nothing
        energies = []

        for g in groupby(gbas, :dimerization)
            color = Makie.wong_colors()[1+(g.dimerization[1]!=0)]
            corners = basis_corners(g.basis[1])
            markerfull = magmarker(corners; fraction = 1, radius = 0.5)
            lines!(ax, g.λ, g.energy_qmc, color = color)
            errorbars!(ax, g.λ, g.energy_qmc, color = color)
            scatter!(ax, g.λ, g.energy_qmc, marker = markerfull, color = :white)

            gf = @rsubset(g, :mag != 0)
            scatter!(
                ax,
                gf.λ,
                gf.energy_qmc,
                marker = [
                    magmarker(
                        corners;
                        fraction = fraction(mag, num_electrons),
                        radius = 0.5,
                    ) for mag in gf.mag
                ],
                color = color,
            )
            scatter!(
                ax,
                g.λ,
                g.energy_qmc;
                marker = markerfull,
                color = :transparent,
                strokecolor = color,
                strokewidth = 1,
            )
            lines!(ax, g.λ, g.energy_hf, color=color, linestyle=:dash)
            λs = g.λ
            push!(energies, g.energy_qmc)
        end

        ΔE = -(energies...)
        for i in eachindex(λs)[1:end-1]
            if ΔE[i] * ΔE[i+1] < 0
                λcross = λs[i] - (λs[i+1] - λs[i]) / (ΔE[i+1] - ΔE[i]) * ΔE[i]
                println("H_$(num_electrons)/$(gbas.basis[1]): λc = $λcross")
                break
            end
        end
    end

    dimerizations = sort(unique(df.dimerization))

    originmarkers = [PolyElement(color = Makie.wong_colors()[n]) for n = 1:2]
    originlabels = [L"δ = 0°", L"δ = %$(360/df.num_electrons[1]*dimerizations[2])°"]
    axislegend(ax, originmarkers, originlabels, titlegap = 3, position = (:center, :bottom))


    mags = sort(unique(df.mag))
    bases = [uppercase(b[end-1:end]) for b in unique(df.basis)]
    markers, labels, titles = maglegend(mags, bases, num_electrons; fraction)

    axislegend(
        ax,
        markers,
        labels,
        [""],
        titlegap = 3,
        orientation = :horizontal,
        nbanks = length(mags),
        position = (:left, :bottom),
        groupgap = 0,
        margin = (3, 3, 3, 3),
    )
end

function fig_hrings_qmc()
    df = @rsubset(load_data("hrings_qmc"), !isnothing(:energy_qmc), :λ < 0.02 && :basis == (:num_electrons == 8 ? "aug-cc-pvtz" : "aug-cc-pvqz"))
    df[!, :energy_qmc] = float.(df.energy_qmc)
    df[!, :energy_qmc_error] = float.(df.energy_qmc_error)

    fig = Figure(size = (390, 370), figure_padding = 3)
    axs = [Axis(fig[i, 1], xlabel = L"$λ$ (a.u.)", ylabel = "Energy (Ha)") for i = 1:2]

    dfs = [@rsubset(df, :num_electrons == n) for n in [4, 8]]
    plot_hrings_qmc(axs[1], dfs[1])
    plot_hrings_qmc(axs[2], dfs[2])

    for (i, (ax, df)) in enumerate(zip(axs, dfs))
        text!(
            ax,
            0,
            1,
            text = rich(
                rich("$(Char(codepoint('a')-1+i))", font = :bold),
                " H",
                subscript(string(Int(df.num_electrons[1]))),
                ", ",
                rich("R", font = :italic),
                " = $(df.R[1]) Å, aug-cc-p$(uppercase(df.basis[1][end-2:end]))",
            ),
            space = :relative,
            align = (:left, :top),
            offset = (4, -2),
        )
    end
    linkxaxes!(axs...)
    hidexdecorations!(axs[1], ticks = false)

    methodmarkers = [LineElement(linestyle=:dash), LineElement()]
    methodlabels = ["UHF", "AFQMC"]
    axislegend(axs[1], methodmarkers, methodlabels, titlegap = 3, position = (:right, :bottom))
    axislegend(axs[2], methodmarkers, methodlabels, titlegap = 3, position = (0.9, 0))

    # ylims!(axs[1], nothing, -1.86)
    # ylims!(axs[2], nothing, -3.9)
    ylims!(axs[1], -2.09, -1.77)
    ylims!(axs[2], -4.3, -3.7)
    rowgap!(fig.layout, 6)
    fig
end

function fig_h2()
    df = load_data("h2")
    fig = Figure(size = (390, 260), figure_padding = 5)
    ax = Axis(fig[1, 1], xlabel = L"$d_\mathrm{HH}$ (\AA)", ylabel = "Energy (Ha)")


    couplings = unique(df.λ)
    δ = couplings[2] - couplings[1]
    coupling_boundaries = vcat(couplings .- δ / 2, [couplings[end] + δ / 2])

    colorrange = extrema(coupling_boundaries)
    colormap = cgrad(
        [Makie.wong_colors()[1], Makie.wong_colors()[2]],
        length(couplings),
        categorical = true,
    )

    dff = @rsubset(df, :basis == "aug-cc-pv5z")
    sort!(dff, :mag; rev = true)

    marker_full = magmarker(5; fraction = 1, radius = 0.47)
    for g in groupby(@rsubset(dff, :mag == 2), [:λ])
        marker = magmarker(5; fraction = g.mag[1] / 2, radius = 0.5)
        lines!(ax, g.R, g.energy_fci; color = g.λ[1], colormap, colorrange)
        if g.mag[1] == 0
            scatter!(
                ax,
                g.R,
                g.energy_fci;
                color = :white,
                colormap,
                colorrange,
                marker = marker_full,
            )
        end
        scatterlines!(
            ax,
            g.R,
            g.energy_fci;
            color = g.λ[1],
            colormap,
            colorrange,
            marker = :pentagon,
            markersize = 10,
        )
    end

    dfs0 = @rsubset(dff, :mag == 0 && :λ == 0)
    dfs1 = @rsubset(dff, :mag == 0 && :λ == maximum(couplings))
    scatterlines!(
        ax,
        dfs0.R,
        dfs0.energy_fci;
        color = colormap.colors[1],
        markercolor = :white,
        strokewidth = 1,
        strokecolor = colormap.colors[1],
        colormap,
        colorrange,
        marker = :pentagon,
        markersize = 10,
    )
    scatterlines!(
        ax,
        dfs1.R,
        dfs1.energy_fci;
        color = colormap.colors[5],
        linestyle = :dash,
        markercolor = :white,
        strokewidth = 1,
        strokecolor = colormap.colors[5],
        colormap,
        colorrange,
        marker = :pentagon,
        markersize = 9,
    )

    markers = [
        MarkerElement(
            markersize = 9,
            marker = :pentagon,
            markercolor = :white,
            strokewidth = 1,
        ),
        MarkerElement(markersize = 10, marker = :pentagon, markercolor = :black),
    ]
    labels = [L"m_\mathrm{S} = 0", L"m_\mathrm{S} = ±1"]

    # axins = Axis(
    #     fig[1, 1],
    #     width = Relative(0.3),
    #     height = Relative(0.25),
    #     halign = 0.9,
    #     valign = 0.56,
    #     ylabel = L"$E-E_\mathrm{5Z}$ (mHa)",
    #     xlabel = "basis",
    #     xticks = ([3, 4, 5], ["TZ", "QZ", "5Z"]),
    #     ylabelpadding = 1,
    #     ylabelsize = 12,
    #     xlabelsize = 12,
    #     xlabelpadding = 1,
    #     xticklabelsize = 12,
    #     yticklabelsize = 12,
    #     xticksize = 3,
    #     yticksize = 3,
    # )
    # dft = select(
    #     @rsubset(df, :basis == "aug-cc-pv5z"),
    #     :mag,
    #     :λ,
    #     :R,
    #     :energy_fci => :energy_fci_t,
    # )
    # df = leftjoin(df, dft, on = [:mag, :λ, :R])
    # df0 = @rsubset(df, :mag == 0 && :basis != "aug-cc-pv5z")
    # df1 = @rsubset(df, :mag == 2 && (:basis != "aug-cc-pv5z" || :R == 2))
    # scatter!(
    #     axins,
    #     basis_order.(df0.basis) + 10 * (df0.λ .- mean(couplings)),
    #     1000 * (df0.energy_fci - df0.energy_fci_t);
    #     marker = magmarker.(basis_order.(df0.basis); fraction = 1.0, radius = 0.5),
    #     color = :white,
    #     strokecolor = [colormap[b/maximum(couplings)] for b in df0.λ],
    #     strokewidth = 1,
    #     colormap,
    #     colorrange,
    # )
    # scatter!(
    #     axins,
    #     basis_order.(df1.basis) + 10 * (df1.λ .- mean(couplings)),
    #     1000 * (df1.energy_fci - df1.energy_fci_t);
    #     marker = magmarker.(basis_order.(df1.basis); fraction = 1.0, radius = 0.5),
    #     color = df1.λ,
    #     colormap,
    #     colorrange,
    # )

    s = 80
    x, y = 40, 100
    img = load("$(path())/assets/h2.png")
    image!(ax, (x, x + s), (y, y + s), img', space = :pixel)
    arrows2d!(
        ax,
        Point2f(x + s / 2 + 10, y + s / 2),
        Point2f(0, 0.25 * s);
        space = :pixel,
        align = :center,
        color = colormap[3],
        tipwidth = 6,
        tiplength = 6,
        shaftwidth = 2,
    )
    text!(
        ax,
        (x + s / 2 + 17, y + s / 2 - 1),
        text = "𝝺",
        space = :pixel,
        align = (:center, :center),
        color = colormap[3],
    )

    text!(
        ax,
        0.12,
        1,
        text = rich("H", subscript("2")),
        align = (:left, :top),
        offset = (3, -6),
        space = :relative,
    )

    axislegend(
        ax,
        markers,
        labels,
        "FCI/aug-cc-pV5Z",
        align = :right,
        gridshalign = :right,
        padding = 0,
    )
    ylims!(ax, -1.25, 0.0)
    Colorbar(fig[1, 2]; colorrange, colormap, label = L"$λ$ (a.u.)")
    fig
end

function fig_origin()
    df = load_data("ring4_origin")
    df0 = @rsubset(df, :origin == 0)
    df08 = @rsubset(df, :origin != 0)

    fig = Figure(size = (390, 280), figure_padding = 6)
    ax = Axis(fig[1, 1], xlabel = L"$λ$ (a. u.)", ylabel = L"$E$ (Ha)")


    num_electrons = 4
    for g in groupby(df, [:mag, :origin, :basis])
        color = Makie.wong_colors()[1+(g.origin[1]!=0)]
        scatter!(
            ax,
            g.λ,
            g.energy_hf,
            marker = magmarker(basis_corners(g.basis[1]); fraction = 1),
            color = :white,
            strokecolor = color,
            strokewidth = 1,
        )
        scatter!(
            ax,
            g.λ,
            g.energy_hf,
            marker = magmarker(
                basis_corners(g.basis[1]);
                fraction = g.mag[1] / num_electrons,
                radius = 0.54,
            ),
            color = color,
        )
    end

    mags = unique(df.mag)
    bases = unique(df.basis)
    markers, labels, titles = maglegend(mags, [uppercase(b[end-1:end]) for b in bases], num_electrons)

    originmarkers = [PolyElement(color = Makie.wong_colors()[n]) for n = 1:2]
    originlabels = [L"(0,0,0)", L"(0.8,0,0)"]

    axislegend(
        ax,
        markers,
        labels,
        titles,
        titlegap = 3,
        orientation = :horizontal,
        nbanks = length(mags),
        position = (:left, :bottom),
        groupgap = 0,
        margin = (3, 3, 3, 3),
    )

    axislegend(ax, originmarkers, originlabels, "Origin (Å)", position = (:center, :bottom))
    text!(
        ax,
        0,
        1,
        text = L"\mathrm{H}_4,~R=%$(df.R[1])\,Å,~δ = 0",
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
    )
    text!(
        ax,
        1,
        1,
        text = "UHF",
        font = :bold,
        align = (:right, :top),
        offset = (-4, -2),
        space = :relative,
    )
    ylims!(ax, nothing, -1.5)

    return fig
end


function fig_horigin()
    df = load_data("h_origin")
    fig = Figure(size = (390, 280), figure_padding = 4)
    ax = Axis(
        fig[1, 1],
        ylabel = L"$E-E_{\mathbf{r}_0 = 0}$ (kcal/mol)",
        xlabel = L"$λ$ (a.u.)",
    )

    df0 = @rsubset(df, :R == 0)
    df = @rsubset(leftjoin(df, df0, on = [:λ, :basis], renamecols = "" => "0"), :R > 0)
    @rtransform!(df, :energy_diff = :energy_fci - :energy_fci0)

    @. model(x, p) = p[2] / x^3 + p[1]

    extrapolated = @by df [:λ, :R] begin
        :extrapolation = Ref(
            curve_fit(
                model,
                float.(basis_order.(:basis))[2:end],
                :energy_diff[2:end],
                [0, :energy_diff[end]],
            ).param,
        )
    end
    λRmodel(x, p) = @views p[1] .* (x[:, 1] .* x[:, 2]) .^ 2

    fit = curve_fit(
        λRmodel,
        hcat(extrapolated.λ, extrapolated.R),
        first.(extrapolated.extrapolation),
        [1.0],
    )


    for g in groupby(extrapolated, :R)
        x = range(extrema(g.λ)..., 50)
        lines!(
            x,
            λRmodel(hcat(x, fill(g.R[1], length(x))), fit.param) * kcalpermolHa,
            color = :black,
            label = L"c |λ×\mathbf{r}_0|^2",
        )
        scatter!(
            ax,
            g.λ,
            first.(g.extrapolation) * kcalpermolHa,
            label = L"\mathbf{r}_0 = %$(@sprintf(\"%.2f\", g.R[1]))~Å~\mathbf{e}_x",
        )
    end

    axins = Axis(
        fig[1, 1],
        width = Relative(0.3),
        height = Relative(0.25),
        halign = 0.15,
        valign = 0.78,
        ylabel = L"$E-E_{\mathbf{r}_0=0}$ (kcal/mol)",
        xlabel = "aug-cc-pV•",
        xticks = ([3, 4, 5, 6], ["TZ", "QZ", "5Z", "6Z"]),
        ylabelpadding = 1,
        ylabelsize = 12,
        xlabelsize = 12,
        xlabelpadding = 1,
        xticklabelsize = 12,
        yticklabelsize = 12,
        xticksize = 3,
        yticksize = 3,
    )
    text!(
        axins,
        1,
        1,
        align = (:right, :top),
        space = :relative,
        text = L"$λ=%$(df.λ[end])$ a.u.",
        offset = (-4, -2),
        fontsize = 12,
    )
    dff = leftjoin(@rsubset(df, :λ == unique(df.λ)[end]), extrapolated, on = [:λ, :R])
    for g in groupby(dff, :R)
        plot!(axins, basis_order.(g.basis), g.energy_diff * kcalpermolHa)
        x = range(3, 6.5, 50)
        lines!(axins, x, model(x, g.extrapolation[1]) * kcalpermolHa)
    end

    @show fit.param

    text!(
        ax,
        1,
        1,
        text = "FCI/extrapolated",
        font = :bold,
        align = (:right, :top),
        offset = (-4, -2),
        space = :relative,
    )
    text!(
        ax,
        0.25,
        1,
        text = "H atom",
        align = (:right, :top),
        offset = (-4, -2),
        space = :relative,
    )
    axislegend(ax, position = (0.6, 0.885), unique = true)
    ylims!(ax, -0.05, 1.8)
    fig

end

function fig_basis_set()
    df = load_data("hrings_hf")
    df = @rsubset(df, :dimerization == 0 && :num_electrons == 6 && :mag in [0,6])
    fig = Figure(size = (390, 280), figure_padding = 4)
    ax = Axis(
        fig[1, 1],
        ylabel = L"$E$ (Ha)",
        xlabel = L"$λ$ (a.u.)",
    )

    sort!(df, :mag)
    for g in groupby(df, [:basis, :mag])
        corners = basis_corners(g.basis[1])
        color = Makie.wong_colors()[g.mag[1]÷g.num_electrons[1]+1]
        lines!(ax, g.λ, g.energy_hf;  color=color)
    end
    for g in groupby(df, [:basis, :mag])
        corners = basis_corners(g.basis[1])
        color = Makie.wong_colors()[g.mag[1]÷g.num_electrons[1]+1]
        @show corners, g.mag[1]
        markerfull = magmarker(corners; fraction=1)
        mag = g.mag[1] != 0
        white = corners == 3 ? :white : :transparent
        
        scatter!(ax, g.λ, g.energy_hf; marker=markerfull, color=mag ? color : white , strokecolor=mag ? :transparent : color, strokewidth = !mag)
    end

    text!(
        ax,
        1,
        1,
        text = "UHF",
        font = :bold,
        align = (:right, :top),
        offset = (-4, -2),
        space = :relative,
    )
    text!(
        ax,
        0.25,
        1,
        text = "H₆",
        align = (:right, :top),
        offset = (-4, -2),
        space = :relative,
    )

    bases = map(x -> uppercase(x[end-1:end]), unique(df.basis))
    mags = sort(unique(df.mag))
    
    markers, labels, titles = maglegend(mags, bases, df.num_electrons[1])



    for m in markers
        m[2][1].markercolor[] = Makie.wong_colors()[2]
        m[2][2].markerstrokecolor[] = Makie.wong_colors()[2]
        m[1][2].markerstrokecolor[] = Makie.wong_colors()[1]
    end
    
    axislegend(
        ax,
        markers,
        labels,
        titles,
        titlegap = 3,
        padding = 0,
        orientation = :horizontal,
        nbanks = length(mags),
        position = (:left, :center),
        groupgap = 0,
        titlehalign = :left,
        margin = (12, 3, 3, 3),
    )
    fig
end

function fig_phases()
    df = load_data("hring_phases")
    @rtransform!(df, :dimerization_angle = 360 * :dimerization / :num_electrons)
    fig = Figure(size = (390, 250))

    ax = Axis(fig[1,1], xlabel=L"$λ_\mathrm{eff}$ (a.u.)", ylabel=L"$M/N$ (a.u.)")

    df4 = @rsubset(df, :num_electrons == 4 && :λ > 0.0015 && :λ < 0.00415)

    df8 = @rsubset(df, :num_electrons == 8 && :λ > 0.0015 && :λ < 0.0042)
    sort!(df8, :λ)
    df8 = @by(df8, [:dimerization_angle, :λ, :basis], begin #remove duplicates
        :mag = :mag[1]
        :magnetic_moment = :magnetic_moment[1]
    end)
    
    lines!(df8.λ, abs.(df8.magnetic_moment), label=rich("H₈, R=0.9 Å, δ=$(df8.dimerization_angle[1])°, ", rich("aug-cc-pVTZ")), color=Makie.wong_colors()[2])
    
    df4_0 = @rsubset(df4, :mag == 0)
    df4_1 = @rsubset(df4, :mag == 2)
    mask0 = df4_0.energy_hf .< df4_1.energy_hf
    mask1 = df4_0.energy_hf .> df4_1.energy_hf

    lines!(ax, df4_0.λ[mask0], abs.(df4_0.magnetic_moment[mask0]), label=rich("H₄, R=0.5 Å, δ=$(df4.dimerization_angle[1])°, ", rich("aug-cc-pVQZ")), color=Makie.wong_colors()[1])
    lines!(ax, df4_1.λ[mask1], abs.(df4_1.magnetic_moment[mask1]), color=Makie.wong_colors()[1])
    scatter!(ax, [df4_0.λ[1], df4_1.λ[end]],[0.0, df8.magnetic_moment[end]], marker=:utriangle, strokecolor=Makie.wong_colors()[2], strokewidth=1, color=:white)
    
    scatter!(ax, collect(extrema(df4_0.λ[mask0])), [0.,0.], marker=:diamond, strokecolor=Makie.wong_colors()[1], strokewidth=1, color=:white)


    halfmarker = magmarker(4; fraction = 0.5, radius = 0.5)
    scatter!(ax, collect(extrema(df4_1.λ[mask1])), df4_1.magnetic_moment[mask1][[1,sum(mask1)]], marker=:diamond, strokecolor=Makie.wong_colors()[1], strokewidth=1, color=:white)
    scatter!(ax, collect(extrema(df4_1.λ[mask1])), df4_1.magnetic_moment[mask1][[1,sum(mask1)]], marker=halfmarker, color=Makie.wong_colors()[1])

    # xlims!(ax, 0.0014, 0.0042)

    text!(0.0018, 0.1, text=L"m_\mathrm{S} = 0")
    text!(0.0028, 2.1, text=L"m_\mathrm{S} = 1")
    text!(0.00365, 3.54, text=L"m_\mathrm{S} = 0")

    axislegend(ax, position=(:left,:top))

    fig
end
