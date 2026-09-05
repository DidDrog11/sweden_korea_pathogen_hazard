using Colors
using ColorBlendModes

function discretize(layer, n::Integer)
    return (x -> round(Int64, x)).(rescale(layer, 1, n))
end

function get_bivariate(dose, uncertainty; high1=colorant"#759d77", high2=colorant"#6676a1", nbreaks = 3)
    # Generate a bivariate palette
    function _palette(; low=colorant"#e8e8e8", high=colorant"#120fe3", breaks=3)
        breakpoints = LinRange(0.0, 1.0, breaks)
        return Makie.ColorSchemes.weighted_color_mean.(breakpoints, high, low)
    end

    # Get the palettes
    colormap1 = _palette(; high=high1, breaks=nbreaks)
    colormap2 = _palette(; high=high2, breaks=nbreaks)
    colormatrix = [ColorBlendModes.blend.(c1, c2; mode=BlendMultiply) for c1 in colormap1, c2 in colormap2]

    # Discrete maps
    m1 = discretize(quantize(dose), nbreaks)
    m2 = discretize(quantize(uncertainty), nbreaks)
    category = similar(m1)
    for i in eachindex(category)
        category[i] = LinearIndices(colormatrix)[m1[i], m2[i]]
    end
    return category, colormatrix
end

function add_bivariate_legend!(location, colormatrix, nbreaks; width=0.28, height=0.28, halign=0.9, valign=0.04,  xlabel="", ylabel="", xlabelsize=12, ylabelsize=12)
    ax_inset = Axis(
        location,
        width=Relative(width),
        height=Relative(height),
        aspect=DataAspect(),
        halign=halign,
        valign=valign,
        xticksvisible=false,
        yticksvisible=false,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        xlabel = xlabel,
        ylabel = ylabel,
        xlabelsize = xlabelsize,
        ylabelsize = ylabelsize
    )
    heatmap!(ax_inset, (reshape(colormatrix, (nbreaks, nbreaks))))
end

function add_colorbar!(location, colorscheme; titlesize=12, titlealign=:left, halign=0.7, valign=0.16, width=0.2, height=0.03, title="")
    cbar_axis = Axis(
        location,
        width=Relative(width),
        height=Relative(height),
        halign=halign,
        valign=valign,
        xticksvisible=false,
        yticksvisible=false,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        title=title,
        titlesize = titlesize,
        titlefont=:regular,
        titlealign=titlealign
    )
    cmap = [get(colorscheme, i) for i in 0:0.025:1]
    X = hcat(0:0.025:1)
    heatmap!(cbar_axis, X, colormap=cmap)
end 


function add_dashes(poly, width=0.05)
    bbox = SDT.boundingbox(poly)
    l,b,t,r = bbox.left - 2, bbox.bottom - 2, bbox.top + 2, bbox.right + 2
    dashes =
        [
            Feature(
                Polygon((l + 0.5i, b), (l + 2 +0.5i, t), (l + 2 + width + 0.5i, t),(l + width + 0.5i, b) ), Dict()
            ) for i in 1:15
        ]
    dashes = [intersect(d, poly) for d in dashes]
    filter(x->!AG.isempty(x.geometry), dashes)
end