path() = "$(@__DIR__)/.."

function load_data(name)
    return DataFrame(JSON.parsefile("$(path())/data/$name.json"))
end

const kcalpermolHa = 627.509394

const marker_cycle =
    [:circle, :diamond, :rect, :utriangle, :ltriangle, :dtriangle, :rtriangle]

function ctheme()
    set_texfont_family!(FontFamily("TeXGyreHeros"))
    return merge(
        Theme(
            figure_padding = 1,
            Axis = (; xgridvisible = false, ygridvisible = false),
            Legend = (;
                framevisible = false,
                labelsize = 12,
                titlesize = 14,
                patchsize = (14, 14),
            ),
            palette = (color = Makie.wong_colors(), marker = marker_cycle),
            Scatter = (cycle = Cycle([:color, :marker], covary = true),),
        ),
    )
end

function basis_corners(basis)
    if endswith(basis, "tz")
        return 3
    elseif endswith(basis, "qz")
        return 4
    elseif endswith(basis, "5z")
        return 5
    end
    error("unknown basis $basis")
end

function magmarker(corners; fraction, radius = 0.53)
    outline = [
        Point2f(radius * cos(φ), radius * sin(φ)) for
        φ in π / 2 .+ range(0, 2π; length = corners + 1)[1:end-1]
    ]

    h = radius * (-1 + 2 * fraction)

    h += 0.001

    i = 1
    while i <= length(outline)
        p1 = outline[i]
        p2 = outline[mod1(i + 1, length(outline))]
        if (p1[2] - h) * (p2[2] - h) < 0
            xnew = p1[1] + (p2[1] - p1[1]) / (p2[2] - p1[2]) * (h - p1[2])
            insert!(outline, i + 1, Point2f(xnew, h))
            i += 1
        end
        i += 1
    end


    outline = filter(x -> x[2] <= h + 0.0001, outline)
    if length(outline) == 0
        return nothing
    end

    path = vcat([MoveTo(outline[end])], [LineTo(p) for p in outline], [ClosePath()])

    return BezierPath(path)
end
