using Pkg
Pkg.activate(".")
using JuMP, Gurobi, PowerModels
using Statistics, StatsBase, Random, Dates, Distributions
using Distributed, ParallelDataTransfer
using CSV, DataFrames, Printf
using JLD2, FileIO
using StatsPlots, PlotThemes

project_root = @__DIR__;
include(joinpath(project_root, "src", "UnitCommitment", "SDDiP", "def.jl"))




case = "case30pwl"; tightness = true; cutSelection = "SMC";
theme(:default)
## ================================================================ gap vs. Time ================================================================ ##
isddlpResult63 = load("src/UnitCommitment/numericalResults-$case/Periods6-Real3/isddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
isddlpResult65 = load("src/UnitCommitment/numericalResults-$case/Periods6-Real5/isddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
isddlpResult610 = load("src/UnitCommitment/numericalResults-$case/Periods6-Real10/isddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
isddlpResult83 = load("src/UnitCommitment/numericalResults-$case/Periods8-Real3/isddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
isddlpResult85 = load("src/UnitCommitment/numericalResults-$case/Periods8-Real5/isddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
isddlpResult810 = load("src/UnitCommitment/numericalResults-$case/Periods8-Real10/isddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
isddlpResult123 = load("src/UnitCommitment/numericalResults-$case/Periods12-Real3/isddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
isddlpResult125 = load("src/UnitCommitment/numericalResults-$case/Periods12-Real5/isddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
isddlpResult1210 = load("src/UnitCommitment/numericalResults-$case/Periods12-Real10/isddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]

sddpResult63 = load("src/UnitCommitment/numericalResults-$case/Periods6-Real3/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult65 = load("src/UnitCommitment/numericalResults-$case/Periods6-Real5/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult610 = load("src/UnitCommitment/numericalResults-$case/Periods6-Real10/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult83 = load("src/UnitCommitment/numericalResults-$case/Periods8-Real3/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult85 = load("src/UnitCommitment/numericalResults-$case/Periods8-Real5/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult810 = load("src/UnitCommitment/numericalResults-$case/Periods8-Real10/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult123 = load("src/UnitCommitment/numericalResults-$case/Periods12-Real3/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult125 = load("src/UnitCommitment/numericalResults-$case/Periods12-Real5/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult1210 = load("src/UnitCommitment/numericalResults-$case/Periods12-Real10/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]


isddlpResult63[!, :gap] = parse.(Float64, replace.(isddlpResult63.gap, r"%" => ""))
isddlpResult65[!, :gap] = parse.(Float64, replace.(isddlpResult65.gap, r"%" => ""))
isddlpResult610[!, :gap] = parse.(Float64, replace.(isddlpResult610.gap, r"%" => ""))
isddlpResult83[!, :gap] = parse.(Float64, replace.(isddlpResult83.gap, r"%" => ""))
isddlpResult85[!, :gap] = parse.(Float64, replace.(isddlpResult85.gap, r"%" => ""))
isddlpResult810[!, :gap] = parse.(Float64, replace.(isddlpResult810.gap, r"%" => ""))
isddlpResult123[!, :gap] = parse.(Float64, replace.(isddlpResult123.gap, r"%" => ""))
isddlpResult125[!, :gap] = parse.(Float64, replace.(isddlpResult125.gap, r"%" => ""))
isddlpResult1210[!, :gap] = parse.(Float64, replace.(isddlpResult1210.gap, r"%" => ""))

sddpResult63[!, :gap] = parse.(Float64, replace.(sddpResult63.gap, r"%" => ""))
sddpResult65[!, :gap] = parse.(Float64, replace.(sddpResult65.gap, r"%" => ""))
sddpResult610[!, :gap] = parse.(Float64, replace.(sddpResult610.gap, r"%" => ""))
sddpResult83[!, :gap] = parse.(Float64, replace.(sddpResult83.gap, r"%" => ""))
sddpResult85[!, :gap] = parse.(Float64, replace.(sddpResult85.gap, r"%" => ""))
sddpResult810[!, :gap] = parse.(Float64, replace.(sddpResult810.gap, r"%" => ""))
sddpResult123[!, :gap] = parse.(Float64, replace.(sddpResult123.gap, r"%" => ""))
sddpResult125[!, :gap] = parse.(Float64, replace.(sddpResult125.gap, r"%" => ""))
sddpResult1210[!, :gap] = parse.(Float64, replace.(sddpResult1210.gap, r"%" => ""))





# 对每个 DataFrame 进行采样，只选择 Iter 为 5 的倍数的行
isddlpResult63 = filter(row -> row.Iter % 5 == 1, isddlpResult63)
isddlpResult65 = filter(row -> row.Iter % 5 == 1, isddlpResult65)
isddlpResult610 = filter(row -> row.Iter % 5 == 1, isddlpResult610)
isddlpResult83 = filter(row -> row.Iter % 5 == 1, isddlpResult83)
isddlpResult85 = filter(row -> row.Iter % 5 == 1, isddlpResult85)
isddlpResult810 = filter(row -> row.Iter % 5 == 1, isddlpResult810)
isddlpResult123 = filter(row -> row.Iter % 5 == 1, isddlpResult123)
isddlpResult125 = filter(row -> row.Iter % 5 == 1, isddlpResult125)
isddlpResult1210 = filter(row -> row.Iter % 5 == 1, isddlpResult1210)

sddpResult63 = filter(row -> row.Iter % 5 == 1, sddpResult63)
sddpResult65 = filter(row -> row.Iter % 5 == 1, sddpResult65)
sddpResult610 = filter(row -> row.Iter % 5 == 1, sddpResult610)
sddpResult83 = filter(row -> row.Iter % 5 == 1, sddpResult83)
sddpResult85 = filter(row -> row.Iter % 5 == 1, sddpResult85)
sddpResult810 = filter(row -> row.Iter % 5 == 1, sddpResult810)
sddpResult123 = filter(row -> row.Iter % 5 == 1, sddpResult123)
sddpResult125 = filter(row -> row.Iter % 5 == 1, sddpResult125)
sddpResult1210 = filter(row -> row.Iter % 5 == 1, sddpResult1210)





timegap = @df isddlpResult63 plot(:Time, :gap, label="SDDℓP-(6,3)", 
                                                title = "Gap vs. Iteration", 
                                                xlab = "Iteration", 
                                                xlim = [0,2000],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                marker=(:xcross, 2, 1.), 
                                                color=:goldenrod,  # 使用蓝色
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"), legend=:outerright,
                                                linestyle=:solid)  # 实线

@df isddlpResult65 plot!(:Time, :gap, marker=(:star, 2, 1.), label="SDDℓP-(6,5)", linestyle=:solid, color=:orange)
# @df isddlpResult610 plot!(:Time, :gap, marker=(:hexagon, 2, 1.), label="SDDℓP-(6,10)", linestyle=:solid, color=:green)
@df isddlpResult83 plot!(:Time, :gap, marker=(:plus, 2, 1.), label="SDDℓP-(8,3)", linestyle=:solid, color=:Crimson)
@df isddlpResult85 plot!(:Time, :gap, marker=(:star, 2, 1.), label="SDDℓP-(8,5)", linestyle=:solid, color=:HotPink)
# @df isddlpResult810 plot!(:Time, :gap, marker=(:hexagon, 2, 1.), label="SDDℓP-(8,10)", linestyle=:solid, color=:brown)
@df isddlpResult123 plot!(:Time, :gap, marker=(:xcross, 2, 1.), label="SDDℓP-(12,3)", linestyle=:solid, color=:lightslategray)
@df isddlpResult125 plot!(:Time, :gap, marker=(:square, 2, 1.), label="SDDℓP-(12,5)", linestyle=:solid, color=:cyan)
# @df isddlpResult1210 plot!(:Time, :gap, marker=(:diamond, 2, 1.), label="SDDℓP-(12,10)", linestyle=:solid, color=:palevioletred)

# sddp 结果，虚线，颜色保持一致
@df sddpResult63 plot!(:Time, :gap, marker=(:xcross, 2, 1.), label="SDDP-(6,3)", linestyle=:dot, color=:goldenrod)
@df sddpResult65 plot!(:Time, :gap, marker=(:star, 2, 1.), label="SDDP-(6,5)", linestyle=:dot, color=:orange)
# @df sddpResult610 plot!(:Time, :gap, marker=(:hexagon, 2, 1.), label="SDDP-(6,10)", linestyle=:dot, color=:green)
@df sddpResult83 plot!(:Time, :gap, marker=(:plus, 2, 1.), label="SDDP-(8,3)", linestyle=:dot, color=:Crimson)
@df sddpResult85 plot!(:Time, :gap, marker=(:star, 2, 1.), label="SDDP-(8,5)", linestyle=:dot, color=:HotPink)
# @df sddpResult810 plot!(:Time, :gap, marker=(:hexagon, 2, 1.), label="SDDP-(8,10)", linestyle=:dot, color=:brown)
@df sddpResult123 plot!(:Time, :gap, marker=(:xcross, 2, 1.), label="SDDP-(12,3)", linestyle=:dot, color=:lightslategray)
@df sddpResult125 plot!(:Time, :gap, marker=(:square, 2, 1.), label="SDDP-(12,5)", linestyle=:dot, color=:cyan)
# @df sddpResult1210 plot!(:Time, :gap, marker=(:diamond, 2, 1.), label="SDDP-(12,10)", linestyle=:dot, color=:palevioletred)





itergap = @df isddlpResult65 plot(:Iter, :gap, label="SDDℓP-(6,5)", 
                                                title = "Relative gaps vs. Iteration", 
                                                xlab = "Iteration", 
                                                # xlim = [0,150],
                                                ylim = [0,70],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                marker=(:xcross, 2, 1.), 
                                                color=:orange,  
                                                size=(700,450),
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"), legend=:outerright,
                                                linestyle=:solid)  
@df isddlpResult63 plot!(:Iter, :gap, marker=(:star, 2, 1.), label="SDDℓP-(6,3)", linestyle=:solid, color=:black)
@df isddlpResult610 plot!(:Iter, :gap, marker=(:hexagon, 2, 1.), label="SDDℓP-(6,10)", linestyle=:solid, color=:purple)
@df isddlpResult83 plot!(:Iter, :gap, marker=(:plus, 2, 1.), label="SDDℓP-(8,3)", linestyle=:solid, color=:ForestGreen)
@df isddlpResult85 plot!(:Iter, :gap, marker=(:star, 2, 1.), label="SDDℓP-(8,5)", linestyle=:solid, color=:red)
@df isddlpResult810 plot!(:Iter, :gap, marker=(:hexagon, 2, 1.), label="SDDℓP-(8,10)", linestyle=:solid, color=:Violet)
@df isddlpResult123 plot!(:Iter, :gap, marker=(:xcross, 2, 1.), label="SDDℓP-(12,3)", linestyle=:solid, color=:lightslategray)
@df isddlpResult125 plot!(:Iter, :gap, marker=(:square, 2, 1.), label="SDDℓP-(12,5)", linestyle=:solid, color=:LightCoral)
@df isddlpResult1210 plot!(:Iter, :gap, marker=(:diamond, 2, 1.), label="SDDℓP-(12,10)", linestyle=:solid, color=:DodgerBlue)


@df sddpResult63 plot!(:Iter, :gap, marker=(:star, 2, 1.), label="SDDP-(6,3)", linestyle=:dashdot, color=:black)
@df sddpResult65 plot!(:Iter, :gap, marker=(:xcross, 2, 1.), label="SDDP-(6,5)", linestyle=:dashdotdot, color=:orange)
@df sddpResult610 plot!(:Iter, :gap, marker=(:hexagon, 2, 1.), label="SDDP-(6,10)", linestyle=:dot, color=:purple)
@df sddpResult83 plot!(:Iter, :gap, marker=(:plus, 2, 1.), label="SDDP-(8,3)", linestyle=:dashdot, color=:ForestGreen)
@df sddpResult85 plot!(:Iter, :gap, marker=(:star, 2, 1.), label="SDDP-(8,5)", linestyle=:dashdotdot, color=:red)
@df sddpResult810 plot!(:Iter, :gap, marker=(:hexagon, 2, 1.), label="SDDP-(8,10)", linestyle=:dot, color=:Violet)
@df sddpResult123 plot!(:Iter, :gap, marker=(:xcross, 2, 1.), label="SDDP-(12,3)", linestyle=:dashdot, color=:lightslategray)
@df sddpResult125 plot!(:Iter, :gap, marker=(:square, 2, 1.), label="SDDP-(12,5)", linestyle=:dashdotdot, color=:LightCoral)
@df sddpResult1210 plot!(:Iter, :gap, marker=(:diamond, 2, 1.), label="SDDP-(12,10)", linestyle=:dot, color=:DodgerBlue)
itergap |> save("/Users/aaron/Downloads/gap_iter.pdf")

