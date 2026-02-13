
function fig_cyclobutadiene()
    df = load_data("cyclobutadiene")

    labels = Dict("d2h" => rich("D", subscript("2h")), "d4h" => rich("D", subscript("4h")))
    colors = Dict("d2h" => Makie.wong_colors()[1], "d4h" => Makie.wong_colors()[2])
    corners = 3
    fig = Figure(size = (390, 280), figure_padding = 5)
    ax = Axis(
        fig[1, 1],
        xlabel = L"$λ$ (a. u.)",
        ylabel = L"$E - E^{D_\mathrm{2h}}_{λ=0,m_\mathrm{S} = 0}$ (kcal/mol)",
    )

    E0 = @rsubset(df, :geometry == "d2h" && :mag == 0 && :λ == 0).energy_qmc |> only

    plots = Dict{Tuple{String,Int},Any}()

    df[!, :energy_diff] = (df.energy_qmc .± df.energy_qmc_error .- E0) * kcalpermolHa
    sort!(df, :mag)
    for g in groupby(df, [:mag, :geometry])
        color = colors[g.geometry[1]]
        sl = scatterlines!(
            ax,
            g.λ,
            g.energy_diff;
            color,
            marker = magmarker(3; fraction = 1, radius = 0.5),
            markercolor = :white,
            strokecolor = color,
            strokewidth = 1,
            label,
        )

        if g.mag[1] == 0 && g.geometry[1] == "d4h"
            Ω = 0.05
        end

        if g.mag[1] != 0
            marker = magmarker(3; fraction = g.mag[1] / 2, radius = 0.48)
            s = scatter!(
                ax,
                g.λ,
                g.energy_diff;
                color,
                marker,
                strokecolor = color,
                strokewidth = 1,
                label,
            )
            sl = [sl, s]
        end
        plots[(g.geometry[1], g.mag[1])] = sl
    end

    text!(
        ax,
        1,
        1,
        text = "AFQMC/aug-cc-pVTZ",
        font = :bold,
        space = :relative,
        align = (:right, :top),
        offset = (-4, -2),
    )
    text!(
        ax,
        0,
        1,
        text = rich("C", subscript("4"), "H", subscript("4")),
        space = :relative,
        align = (:left, :top),
        offset = (4, -2),
    )

    gap0 = -(@rsubset(df, :mag == 0 && :λ == 0).energy_diff...)
    gapλ = -(
        @rsubset(
            df,
            (:geometry == "d2h" && :mag == 0 || :geometry == "d4h" && :mag == 1) &&
            :λ == 0.02
        ).energy_diff...
    )
    @show gap0
    @show gapλ

    df1 = @rsubset(df, :geometry == "d2h" && :mag == 0)
    df2 = @rsubset(df, :geometry == "d4h" && :mag == 1)
    @assert df1.λ == df2.λ
    ΔE = df1.energy_diff - df2.energy_diff

    for i in eachindex(ΔE)[1:end-1]
        if ΔE[i] * ΔE[i+1] < 0
            λcross = df1.λ[i] - (df1.λ[i+1] - df1.λ[i]) / (ΔE[i+1] - ΔE[i]) * ΔE[i]
            @show λcross
            break
        end
    end

    mags = unique(df.mag)
    geoms = ["d2h", "d4h"]
    pm(m) = m == 0 ? "" : "\\!±\\!"
    axislegend(
        ax,
        [[plots[(geom, m)] for m in mags] for geom in geoms],
        [[L"m_\mathrm{S}\!=\!%$(pm(m))%$(round(Int,m))" for m in mags] for geom in geoms],
        getindex.(Ref(labels), geoms),
        position = (1, 0.9),
        orientation = :horizontal,
        nbanks = 3,
        titleposition = :left,
        padding = 0,
    )
    xlims!(ax, -0.0005, 0.0208)
    ylims!(ax, -5, 231)

    img = load("$(path())/assets/c4h4.png")
    s = 70
    x = 20
    y = 60
    image!(ax, (x, x + s), (y, y + s), img, space = :pixel)
    colormap =
        cgrad([Makie.wong_colors()[1], Makie.wong_colors()[2]], 5, categorical = true)
    scatter!(
        ax,
        [x + s / 2],
        [y + s / 2],
        marker = '⊗',
        color = colormap[3],
        markersize = 12,
        space = :pixel,
    )
    scatter!(
        ax,
        [x + s / 2 + 8],
        [y + s / 2 - 2],
        marker = '𝝺',
        color = colormap[3],
        markersize = 14,
        space = :pixel,
    )
    fig
end


# Note on the symmetric case:
#
# We assume a paraboloid that is symmetric under dcc1<->dcc2,
# which is reasonable due to the symmetry of the molecule.
# However, this means that our fit is insensitive to the possibility that,
# due to deviations from a paraboloid, the minimum is away from the symmetric line.
# We will take the fact that the paraboloid fits the data as an indirect proof that this is not the case.
function paraboloid_model(xy, p; symmetric)
    e0 = p[1]
    M = [p[2]; p[3];; p[3]; symmetric ? p[2] : p[6]]
    dxy0 = symmetric ? [p[4],p[4]] : p[4:5] 

    return [e0 + @views dot(xy[i, :] - dxy0, M, xy[i, :] - dxy0) for i in axes(xy, 1)]
end

function fit_paraboloid(df; symmetric)
    function model(xy, p)
        paraboloid_model(xy, p; symmetric)
    end
    idx = argmin(df.energy_diff)
    dcc1 = df.dcc1[idx]
    dcc2 = df.dcc2[idx]

    if symmetric
        p0 = [minimum(Measurements.value, df.energy_qmc), 1.45, 0.1, dcc1]
    else
        p0 = [minimum(Measurements.value, df.energy_qmc), 1.45, 0.1, dcc1, dcc2, 1.45]
    end
    fit = curve_fit(
        model,
        hcat(df.dcc1, df.dcc2),
        Measurements.value.(df.energy_qmc),
        1 ./ Measurements.uncertainty.(df.energy_qmc) .^ 2,
        p0,
    )
    param = fit.param .± stderror(fit)
    
    # see model definition
    if symmetric
        push!(param, param[4])
    end

    χsq = rss(fit)

    return param, χsq / dof(fit)
end

function plot_surface(ax, df; color, marker)
    colormap = cgrad([color, :white])


    emin, emax = Measurements.value.(extrema(df.energy_diff))
    colors =
        [colormap[((e-emin)/(emax-emin))^0.6] for e in Measurements.value.(df.energy_diff)]

    points = Point3d.(df.dcc1, df.dcc2, Measurements.value.(df.energy_diff))

    xs = round.(Int, sum.(df.perturbation))
    ys = round.(Int, map(x -> x[1] - x[2], df.perturbation))

    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)

    xs .+= -xmin + 1
    ys .+= -ymin + 1

    idxs = fill(-1, xmax - xmin + 1, ymax - ymin + 1)
    for (i, (x, y)) in enumerate(zip(xs, ys))
        idxs[x, y] = i
    end

    faces = Vector{Int}[]

    for y = 1:size(idxs, 2)-1, x = 1:size(idxs, 1)-1
        i00 = idxs[x, y]
        i10 = idxs[x+1, y]
        i01 = idxs[x, y+1]
        i11 = idxs[x+1, y+1]

        covered = filter(!=(-1), [i00, i10, i11, i01])
        if length(covered) < 3
            continue
        elseif length(covered) == 3
            push!(faces, covered)
        else
            # choose convex triangularization
            if (df.energy_diff[i00] + df.energy_diff[i11]) / 2 <
               (df.energy_diff[i01] + df.energy_diff[i10]) / 2
                tri1a = [i00, i10, i11]
                tri1b = [i00, i11, i01]
                push!(faces, tri1a)
                push!(faces, tri1b)
            else
                tri2a = [i00, i10, i01]
                tri2b = [i11, i01, i10]
                push!(faces, tri2a)
                push!(faces, tri2b)
            end
        end

    end

    faces = mapreduce(transpose, vcat, faces)


    mesh!(points, faces, color = colors, shading = false)

    symmetric = select(@rsubset(df, :dcc1 == :dcc2), :dcc1, :dcc2, :energy_diff)
    if !isempty(symmetric)
        sort!(symmetric, :dcc1; rev = true)
        #interpolate
        f = 0.62
        insert!(
            symmetric,
            2,
            Dict(
                name => f * symmetric[1, name] + (1 - f) * symmetric[2, name] for
                name in [:dcc1, :dcc2, :energy_diff]
            ),
        )

        lines!(
            ax,
            symmetric.dcc1[1:2],
            symmetric.dcc2[1:2],
            Measurements.value.(symmetric.energy_diff[1:2]);
            color,
        )
        lines!(
            ax,
            symmetric.dcc1[2:6],
            symmetric.dcc2[2:6],
            Measurements.value.(symmetric.energy_diff[2:6]);
            color,
            linestyle = :dash,
        )
        lines!(
            ax,
            symmetric.dcc1[6:end],
            symmetric.dcc2[6:end],
            Measurements.value.(symmetric.energy_diff[6:end]);
            color,
        )
    end
    center = @rsubset(df, iszero(:perturbation))

    # dcc1 > dcc2 is a copy of dcc1 < dcc2 that was added for the plot and thus correlated
    # therefore we delete it again before doing the fit, which assumes symmetry in the model.
    dffit = df
    if df.geometry[1] == "d4h"
        dffit = @rsubset(df, :dcc1 <= :dcc2 + 1e-8)
    end
    
    param, χ²dof = fit_paraboloid(dffit; symmetric = df.geometry[1] == "d4h")
    Emin, dcc1fit, dcc2fit = param[1:3] 
    println("$(df.geometry[1]): Emin = $Emin, dcc1 = $dcc1fit, dcc2 = $dcc2fit")
    @show χ²dof

    xs = LinRange(extrema(df.dcc1)..., 100)
    ys = LinRange(extrema(df.dcc2)..., 100)

    scatter!(ax, points; overdraw = true, color, marker)
    scatter!(
        ax,
        center.dcc1,
        center.dcc2,
        Measurements.value.(center.energy_diff);
        color,
        marker,
        markersize = 10,
        overdraw = true,
        strokewidth = 1,
    )
    return Emin
end

function fig_cyclobutadiene_perturbation()
    df = load_data("cyclobutadiene_perturbation")

    df.energy_qmc .= df.energy_qmc .± df.energy_qmc_error
    E0 = minimum(df.energy_qmc)
    df[!, :energy_diff] = (df.energy_qmc .- E0) * kcalpermolHa

    fig = Figure(size = (390, 280), figure_padding = (20, 20, 0, 0))


    ax = Axis3(
        fig[1, 1],
        xlabel = L"d_\mathrm{CC,1}~(\AA)",
        ylabel = L"d_\mathrm{CC,2}~(\AA)",
        zlabel = L"$E-E_\mathrm{min}$ (kcal/mol)",
        xlabeloffset = 30,
        ylabeloffset = 30,
        zlabeloffset = 30,
        azimuth = 1.15π,
        elevation = 0.05π,
    )

    df = vcat(
        df,
        @rtransform(
            @rsubset(df, :geometry == "d4h" && !allequal(:perturbation)),
            :dcc1 = :dcc2,
            :dcc2 = :dcc1,
            :perturbation = reverse(:perturbation)
        )
    )

    EminD2h = plot_surface(
        ax,
        @rsubset(df, :geometry == "d2h"),
        color = Makie.wong_colors()[1],
        marker = :diamond,
    )
    EminD4h = plot_surface(
        ax,
        @rsubset(df, :geometry == "d4h"),
        color = Makie.wong_colors()[2],
        marker = :circle,
    )

    EDh40 = minimum(df.energy_qmc)
    EDh20 = minimum(@rsubset(df, :geometry == "d2h").energy_qmc)
    Ebarrier = (EminD2h - EminD4h) * kcalpermolHa
    @show Ebarrier

    number =
        @sprintf("%.2g(%1.0f)\nkcal/mol", Ebarrier, Measurements.uncertainty(10 * Ebarrier))
    bracket!(
        ax,
        0.23,
        0.46,
        0.23,
        0.26;
        space = :relative,
        style = :square,
        justification = :right,
        orientation = :down,
        width = 5,
        linewidth = 1,
        text = number,
        fontsize = 10,
        textcolor = Makie.wong_colors()[1],
        rotation = 0,
        textoffset = 22,
    )

    text!(
        ax,
        1.47,
        1.34,
        7.2,
        text = rich("AFQMC/aug-cc-pVTZ", font = :bold),
        space = :data,
        markerspace = :data,
        align = (:right, :baseline),
        fontsize = 0.12,
        rotation = (0.5, -0.5, -0.5, 0.5),
    )
    text!(
        ax,
        1.47,
        1.34,
        6.3,
        text = rich(rich("λ", font = :italic), " = $(df.λ[1]) a.u."),
        space = :data,
        markerspace = :data,
        align = (:right, :baseline),
        fontsize = 0.12,
        rotation = (0.5, -0.5, -0.5, 0.5),
    )
    zlims!(ax, 0, 8)
    xlims!(ax, 1.32, 1.48)
    ylims!(ax, 1.34, 1.61)
    resize_to_layout!(fig)
    return fig
end


function plot_energy_heatmap(ax, df; color, E0, marker)
    colormap = cgrad([color, :white])
    emin, emax = Measurements.value.(extrema(df.energy_diff))

    symmetric =  df.geometry[1] == "d4h"

    param, χ²dof = fit_paraboloid(df; symmetric)
    Emin, dcc1fit, dcc2fit = param[[1,4,5]] 
    println("$(df.geometry[1]): Emin = $Emin, dcc1 = $dcc1fit, dcc2 = $dcc2fit")
    @show χ²dof

    if symmetric
        df = vcat(df, select(df, :dcc1 => :dcc2, :dcc2 => :dcc1, Not([:dcc1, :dcc2])))
    end
    

    xs = range(1.31, 1.49, 100)
    ys = range(1.35, 1.63, 100)
    zs = reshape(paraboloid_model(mapreduce(collect,hcat,Iterators.product(xs,ys))', Measurements.value.(param);symmetric), length(xs), length(ys))
    zs = (zs .- E0) * kcalpermolHa

    emax += 0.7
    df.colors =
        [colormap[((e-emin)/(emax-emin))] for e in Measurements.value.(df.energy_diff)]
    if symmetric
        for iy = axes(zs,2)
            for ix = axes(zs,1)
                if xs[ix] - ys[iy] < -0.1
                    zs[ix,iy] = NaN
                end
            end
        end
    end


    xlims!(ax, extrema(xs) .+(0.005,-0.005)...)
    heatmap!(ax, xs, ys, zs; colormap=colormap, colorrange=(emin,emax), interpolate=true)
    
    scatter!(ax, df.dcc1, df.dcc2; strokecolor=:black, marker, color=df.colors, strokewidth=1)
    lines!(ax, xs[1:end], 0.1 .+ xs[1:end], color=:gray)

    center = @rsubset(df, iszero(:perturbation))
    scatter!(
        ax,
        center.dcc1,
        center.dcc2;
        color=center.colors,
        marker,
        markersize = 10,
        overdraw = true,
        strokewidth = 2,
    )

    scatter!(
        ax,
        [dcc1fit],
        [dcc2fit];
        color=:white,
        marker=:cross,
    )

    return colormap, (emin,emax), Emin
end

function fig_cyclobutadiene_perturbation_heatmap()
    df = load_data("cyclobutadiene_perturbation")

    df.energy_qmc .= df.energy_qmc .± df.energy_qmc_error
    E0 = minimum(df.energy_qmc)
    df[!, :energy_diff] = (df.energy_qmc .- E0) * kcalpermolHa

    fig = Figure(size = (390, 280))
    ax = Axis(fig[1,1],
        xlabel = L"d_\mathrm{CC,1}~(\AA)",
        ylabel = L"d_\mathrm{CC,2}~(\AA)",
    )


    cmap1, crange1, Emin1 = plot_energy_heatmap(ax, @rsubset(df, :geometry == "d2h"); color = Makie.wong_colors()[1], marker=:circle, E0=Measurements.value(E0))
    cmap2, crange2, Emin2 = plot_energy_heatmap(ax, @rsubset(df, :geometry == "d4h"); color = Makie.wong_colors()[2], marker=:diamond, E0=Measurements.value(E0))

    gap = (Emin1-Emin2)*kcalpermolHa
    @show gap

    fakecmap1 = cgrad([:white, :white, cmap1.colors...],[crange2[1], crange1[1]-0.0025,crange1[1]+0.005,crange1[2]]./crange1[2])
    fakecmap2 = cgrad([cmap2.colors..., :white],[crange2..., crange1[2]]./crange1[2])
    fakecrange = (crange2[1],crange1[2])

    Colorbar(fig[1,2], colormap=fakecmap1,colorrange=fakecrange, ticksvisible = false,ticklabelsvisible=false, nsteps=200)
    axbar = Axis(fig[1,2], limits=(0..1, fakecrange))
	hidedecorations!(axbar)
	# arrows2d!(axbar, [0.5], [0.0], [0.0], [crange1[1]], shaftwidth=2, tipwidth=8)
	lines!(axbar, [0,1], [gap.val, gap.val],color=:black)
	errorbars!(axbar, [0.5], [gap], color=:black, whiskerwidth=5)
    Colorbar(fig[1,3], colormap=fakecmap2,colorrange=fakecrange, label=L"$E-E_\mathrm{min}$ (kcal/mol)")

    text!(
        ax,
        1,
        1,
        text = rich(rich("AFQMC/aug-cc-pVTZ\n", font = :bold), rich("λ", font = :italic), " = $(df.λ[1]) a.u."),
        space = :relative,
        align = (:right, :top),
        offset = (-4,-2)
    )

    text!(
        ax,
        0.65,
        0.72,
        fontsize=14,
        text = L"m_\mathrm{S} = 0",
        space = :relative,
        align = (:center, :center),
    )
    text!(
        ax,
        0.82,
        0.62,
        text = L"m_\mathrm{S} = \pm 1",
        space = :relative,
        align = (:center, :center),
    )
    
    colgap!(fig.layout,1,Relative(0.02))
    colgap!(fig.layout,2,Relative(0.0))
    fig
end
