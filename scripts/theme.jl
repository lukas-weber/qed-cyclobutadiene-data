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
            MeasPlot = (cycle = Cycle([:color, :marker], covary = true),),
            Scatter = (cycle = Cycle([:color, :marker], covary = true),),
        ),
    )
end
