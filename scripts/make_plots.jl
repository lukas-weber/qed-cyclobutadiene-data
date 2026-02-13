using CairoMakie
using DataFrames
using DataFramesMeta
using FileIO: load
using JSON
using LaTeXStrings
using LinearAlgebra
using LsqFit
using MathTeXEngine
using Measurements
using Printf
using Statistics



include("common.jl")
include("theme.jl")
include("hring.jl")
include("cyclobutadiene.jl")

function make_plots()
    function saveplot(name, fig)
        save("$(path())/plots/$name.pdf", fig, pt_per_unit = 0.625)
    end

    with_theme(ctheme()) do
        saveplot("hrings", fig_hrings())
        saveplot("hrings_qmc", fig_hrings_qmc())
        saveplot("h2", fig_h2())
        saveplot("cyclobutadiene", fig_cyclobutadiene())
        saveplot("cyclobutadiene_perturbation", fig_cyclobutadiene_perturbation_heatmap())
        saveplot("origin", fig_origin())
        saveplot("h_origin", fig_horigin())
        saveplot("basis_set", fig_basis_set())
    end
end

function (@main)(_)
    make_plots()
    return nothing
end
