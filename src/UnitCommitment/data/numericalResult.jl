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
sddlpResult63 = load("src/UnitCommitment/numericalResults-$case/Periods6-Real3/sddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
sddlpResult65 = load("src/UnitCommitment/numericalResults-$case/Periods6-Real5/sddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
sddlpResult610 = load("src/UnitCommitment/numericalResults-$case/Periods6-Real10/sddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
sddlpResult83 = load("src/UnitCommitment/numericalResults-$case/Periods8-Real3/sddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
sddlpResult85 = load("src/UnitCommitment/numericalResults-$case/Periods8-Real5/sddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
sddlpResult810 = load("src/UnitCommitment/numericalResults-$case/Periods8-Real10/sddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
sddlpResult123 = load("src/UnitCommitment/numericalResults-$case/Periods12-Real3/sddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
sddlpResult125 = load("src/UnitCommitment/numericalResults-$case/Periods12-Real5/sddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]
sddlpResult1210 = load("src/UnitCommitment/numericalResults-$case/Periods12-Real10/sddlpResult-$cutSelection-true.jld2")["sddlpResult"][:solHistory]

sddpResult63 = load("src/UnitCommitment/numericalResults-$case/Periods6-Real3/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult65 = load("src/UnitCommitment/numericalResults-$case/Periods6-Real5/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult610 = load("src/UnitCommitment/numericalResults-$case/Periods6-Real10/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult83 = load("src/UnitCommitment/numericalResults-$case/Periods8-Real3/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult85 = load("src/UnitCommitment/numericalResults-$case/Periods8-Real5/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult810 = load("src/UnitCommitment/numericalResults-$case/Periods8-Real10/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult123 = load("src/UnitCommitment/numericalResults-$case/Periods12-Real3/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult125 = load("src/UnitCommitment/numericalResults-$case/Periods12-Real5/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]
sddpResult1210 = load("src/UnitCommitment/numericalResults-$case/Periods12-Real10/sddpResult-$cutSelection-true.jld2")["sddpResult"][:solHistory]


sddlpResult63[!, :gap] = parse.(Float64, replace.(sddlpResult63.gap, r"%" => ""))
sddlpResult65[!, :gap] = parse.(Float64, replace.(sddlpResult65.gap, r"%" => ""))
sddlpResult610[!, :gap] = parse.(Float64, replace.(sddlpResult610.gap, r"%" => ""))
sddlpResult83[!, :gap] = parse.(Float64, replace.(sddlpResult83.gap, r"%" => ""))
sddlpResult85[!, :gap] = parse.(Float64, replace.(sddlpResult85.gap, r"%" => ""))
sddlpResult810[!, :gap] = parse.(Float64, replace.(sddlpResult810.gap, r"%" => ""))
sddlpResult123[!, :gap] = parse.(Float64, replace.(sddlpResult123.gap, r"%" => ""))
sddlpResult125[!, :gap] = parse.(Float64, replace.(sddlpResult125.gap, r"%" => ""))
sddlpResult1210[!, :gap] = parse.(Float64, replace.(sddlpResult1210.gap, r"%" => ""))

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
sddlpResult63 = filter(row -> row.Iter % 5 == 1, sddlpResult63)
sddlpResult65 = filter(row -> row.Iter % 5 == 1, sddlpResult65)
sddlpResult610 = filter(row -> row.Iter % 5 == 1, sddlpResult610)
sddlpResult83 = filter(row -> row.Iter % 5 == 1, sddlpResult83)
sddlpResult85 = filter(row -> row.Iter % 5 == 1, sddlpResult85)
sddlpResult810 = filter(row -> row.Iter % 5 == 1, sddlpResult810)
sddlpResult123 = filter(row -> row.Iter % 5 == 1, sddlpResult123)
sddlpResult125 = filter(row -> row.Iter % 5 == 1, sddlpResult125)
sddlpResult1210 = filter(row -> row.Iter % 5 == 1, sddlpResult1210)

sddpResult63 = filter(row -> row.Iter % 5 == 1, sddpResult63)
sddpResult65 = filter(row -> row.Iter % 5 == 1, sddpResult65)
sddpResult610 = filter(row -> row.Iter % 5 == 1, sddpResult610)
sddpResult83 = filter(row -> row.Iter % 5 == 1, sddpResult83)
sddpResult85 = filter(row -> row.Iter % 5 == 1, sddpResult85)
sddpResult810 = filter(row -> row.Iter % 5 == 1, sddpResult810)
sddpResult123 = filter(row -> row.Iter % 5 == 1, sddpResult123)
sddpResult125 = filter(row -> row.Iter % 5 == 1, sddpResult125)
sddpResult1210 = filter(row -> row.Iter % 5 == 1, sddpResult1210)





timegap = @df sddlpResult63 plot(:Time, :gap, label="SDDℓP-(6,3)", 
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

@df sddlpResult65 plot!(:Time, :gap, marker=(:star, 2, 1.), label="SDDℓP-(6,5)", linestyle=:solid, color=:orange)
# @df sddlpResult610 plot!(:Time, :gap, marker=(:hexagon, 2, 1.), label="SDDℓP-(6,10)", linestyle=:solid, color=:green)
@df sddlpResult83 plot!(:Time, :gap, marker=(:plus, 2, 1.), label="SDDℓP-(8,3)", linestyle=:solid, color=:Crimson)
@df sddlpResult85 plot!(:Time, :gap, marker=(:star, 2, 1.), label="SDDℓP-(8,5)", linestyle=:solid, color=:HotPink)
# @df sddlpResult810 plot!(:Time, :gap, marker=(:hexagon, 2, 1.), label="SDDℓP-(8,10)", linestyle=:solid, color=:brown)
@df sddlpResult123 plot!(:Time, :gap, marker=(:xcross, 2, 1.), label="SDDℓP-(12,3)", linestyle=:solid, color=:lightslategray)
@df sddlpResult125 plot!(:Time, :gap, marker=(:square, 2, 1.), label="SDDℓP-(12,5)", linestyle=:solid, color=:cyan)
# @df sddlpResult1210 plot!(:Time, :gap, marker=(:diamond, 2, 1.), label="SDDℓP-(12,10)", linestyle=:solid, color=:palevioletred)

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





itergap = @df sddlpResult65 plot(:Iter, :gap, label="SDDℓP-(6,5)", 
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
@df sddlpResult63 plot!(:Iter, :gap, marker=(:star, 2, 1.), label="SDDℓP-(6,3)", linestyle=:solid, color=:black)
@df sddlpResult610 plot!(:Iter, :gap, marker=(:hexagon, 2, 1.), label="SDDℓP-(6,10)", linestyle=:solid, color=:purple)
@df sddlpResult83 plot!(:Iter, :gap, marker=(:plus, 2, 1.), label="SDDℓP-(8,3)", linestyle=:solid, color=:ForestGreen)
@df sddlpResult85 plot!(:Iter, :gap, marker=(:star, 2, 1.), label="SDDℓP-(8,5)", linestyle=:solid, color=:red)
@df sddlpResult810 plot!(:Iter, :gap, marker=(:hexagon, 2, 1.), label="SDDℓP-(8,10)", linestyle=:solid, color=:Violet)
@df sddlpResult123 plot!(:Iter, :gap, marker=(:xcross, 2, 1.), label="SDDℓP-(12,3)", linestyle=:solid, color=:lightslategray)
@df sddlpResult125 plot!(:Iter, :gap, marker=(:square, 2, 1.), label="SDDℓP-(12,5)", linestyle=:solid, color=:LightCoral)
@df sddlpResult1210 plot!(:Iter, :gap, marker=(:diamond, 2, 1.), label="SDDℓP-(12,10)", linestyle=:solid, color=:DodgerBlue)


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


## ================================================================ time vs. Iter ================================================================ ##
T = 8; num = 5; tightness = true; cutSelection = "SMC";
sddlpResultLC = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddlpResult-LC-$tightness.jld2")["sddlpResult"][:solHistory]
sddlpResultSMC = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddlpResult-SMC-$tightness.jld2")["sddlpResult"][:solHistory]
sddlpResultELC = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/MagnantiWong/isddlpHCResult-0.0-$tightness.jld2")["sddlpResult"][:solHistory]
sddpResultSMC = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddpResult-SMC-$tightness.jld2")["sddpResult"][:solHistory]
sddpResultLC = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddpResult-LC-$tightness.jld2")["sddpResult"][:solHistory]
sddpResultELC = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/MagnantiWong/sddpResult-0.0-$tightness.jld2")["sddpResult"][:solHistory]

# 对每个 DataFrame 进行采样，只选择 Iter 为 3 的倍数的行
sddlpResultLC = filter(row -> row.Iter % 2 == 1, sddlpResultLC)
sddlpResultSMC = filter(row -> row.Iter % 2 == 1, sddlpResultSMC)
sddlpResultELC = filter(row -> row.Iter % 2 == 1, sddlpResultELC)
sddpResultSMC = filter(row -> row.Iter % 2 == 1, sddpResultSMC)
sddpResultLC = filter(row -> row.Iter % 2 == 1, sddpResultLC)
sddpResultELC = filter(row -> row.Iter % 2 == 1, sddpResultELC)


timeiter = @df sddlpResultLC plot(:Iter, :time, label="SDDℓP-LC", 
                                                title = "Iteration time (sec.) vs. Iteration", 
                                                xlab = "Iteration", 
                                                # ylab = "Iteration time",
                                                # ylim = [20000,23000],
                                                xlim = [0,150],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                # size=(400,300),
                                                marker=(:circle, 2, 1.), 
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"), legend=:outerright);
@df sddlpResultSMC plot!(:Iter, :time, marker=(:star, 2, 1.), label="SDDℓP-SMC")
@df sddlpResultELC plot!(:Iter, :time, marker=(:hexagon, 2, 1.), label="SDDℓP-PLC")
@df sddpResultLC plot!(:Iter, :time, marker=(:circle, 2, 1.), label="SDDP-LC")
@df sddpResultSMC plot!(:Iter, :time, marker=(:circle, 2, 1.), label="SDDP-SMC")
@df sddpResultELC plot!(:Iter, :time, marker=(:circle, 2, 1.), label="SDDP-PLC")



## ================================================================ binarization ================================================================ ##
coef = 10; T = 6; num = 5; case = "case30";
sddipResult6 = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-6.jld2")["sddipResult"][:solHistory]
sddipResult7 = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-7.jld2")["sddipResult"][:solHistory]
sddipResult8 = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-8.jld2")["sddipResult"][:solHistory]
sddipResult9 = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-9.jld2")["sddipResult"][:solHistory]
sddipResult10 = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-10.jld2")["sddipResult"][:solHistory]

sddipResult5[!, :gap] = parse.(Float64, replace.(sddipResult5.gap, r"%" => ""))
sddipResult6[!, :gap] = parse.(Float64, replace.(sddipResult6.gap, r"%" => ""))
sddipResult7[!, :gap] = parse.(Float64, replace.(sddipResult7.gap, r"%" => ""))
sddipResult8[!, :gap] = parse.(Float64, replace.(sddipResult8.gap, r"%" => ""))
sddipResult9[!, :gap] = parse.(Float64, replace.(sddipResult9.gap, r"%" => ""))
sddipResult10[!, :gap] = parse.(Float64, replace.(sddipResult10.gap, r"%" => ""))

sddipResult5 = filter(row -> row.Iter % 2 == 1, sddipResult5)
sddipResult6 = filter(row -> row.Iter % 2 == 1, sddipResult6)
sddipResult7 = filter(row -> row.Iter % 2 == 1, sddipResult7)
sddipResult8 = filter(row -> row.Iter % 2 == 1, sddipResult8)
sddipResult9 = filter(row -> row.Iter % 2 == 1, sddipResult9)
sddipResult10 = filter(row -> row.Iter % 2 == 1, sddipResult10)

sddipResult5[!, :LB] = sddipResult5[!, :LB] * 6
sddipResult6[!, :LB] = sddipResult6[!, :LB] * 6
sddipResult7[!, :LB] = sddipResult7[!, :LB] * 6
sddipResult8[!, :LB] = sddipResult8[!, :LB] * 6
sddipResult9[!, :LB] = sddipResult9[!, :LB] * 6
sddipResult10[!, :LB] = sddipResult10[!, :LB] * 6

sddipResult5[!, :UB] = sddipResult5[!, :UB] * 6
sddipResult6[!, :UB] = sddipResult6[!, :UB] * 6
sddipResult7[!, :UB] = sddipResult7[!, :UB] * 6
sddipResult8[!, :UB] = sddipResult8[!, :UB] * 6
sddipResult9[!, :UB] = sddipResult9[!, :UB] * 6
sddipResult10[!, :UB] = sddipResult10[!, :UB] * 6


lbtime = @df sddipResult6 plot(:Time, :LB, label="1/64", 
                                                title = "Upper and lower bounds vs. Time", 
                                                xlab = "Time (sec.)", 
                                                ylab = "Lower bounds (× 10³)",
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=10, 
                                                # size=(400,300),
                                                yformatter = y -> y/10^3,
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"), 
                                                marker=(:circle, 2, 1.), 
                                                linewidth=1,
                                                legend=:bottomright,
                                                color=:DodgerBlue)  # Use blue color for 1/64
@df sddipResult6 plot!(:Time, :UB, label="", linewidth=1, marker=(:circle, 2, 1.), color=:DodgerBlue)  # Plot UB for 1/64, no label

# Plot for 1/128 (LB and UB with the same color)
@df sddipResult7 plot!(:Time, :LB, label="1/128", linewidth=1, marker=(:star, 2, 1.), color=:red)  # Use red color for 1/128
@df sddipResult7 plot!(:Time, :UB, label="", linewidth=1, marker=(:star, 2, 1.), color=:red)  # Plot UB for 1/128, no label

@df sddipResult8 plot!(:Time, :LB, linewidth=1, marker=(:hexagon, 2, 1.), label="1/256", color=:orange)
@df sddipResult8 plot!(:Time, :UB, linewidth=1, marker=(:hexagon, 2, 1.), label="", color=:orange)

@df sddipResult9 plot!(:Time, :LB, linewidth=1, marker=(:xcross, 2, 1.), label="1/512", color=:green)
@df sddipResult9 plot!(:Time, :UB, linewidth=1, marker=(:xcross, 2, 1.), label="", color=:green)

@df sddipResult10 plot!(:Time, :LB, linewidth=1, marker=(:plus, 2, 1.), label="1/1024", color=:purple)
@df sddipResult10 plot!(:Time, :UB, linewidth=1, marker=(:plus, 2, 1.), label="", color=:purple)
lbtime |> save("/Users/aaron/Downloads/SMC-sddip_binarization_comparsion_time.pdf")