using Pkg
Pkg.activate(".")
using JuMP, Gurobi, PowerModels
using Statistics, StatsBase, Random, Dates, Distributions
using Distributed, ParallelDataTransfer
using CSV, DataFrames, Printf
using JLD2, FileIO
using StatsPlots, PlotThemes

project_root = @__DIR__;
include(joinpath(project_root, "src", "alg", "utilities", "structs.jl"))
theme(:default)

case = "case30pwl"; 
tightness = true; 
cutSelection = :SMC; 
num = 10; T = 12; 
algorithm = :SDDPL;

for cut in [:LC, :PLC, :SMC]
    for T in [6, 8, 12] 
        for num in [3, 5, 10] 
            solHistory = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]
        end
    end
end


## ----------------------------------------------------------------------------------------------- ##
## the same cut with different instances
cutSelection = :SMC
sddlpResult63 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods6-Real3/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]
sddlpResult65 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods6-Real5/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]
sddlpResult610 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods6-Real10/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]

sddlpResult83 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods8-Real3/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]
sddlpResult85 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods8-Real5/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]
sddlpResult810 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods8-Real10/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]

sddlpResult123 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods12-Real3/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]
sddlpResult125 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods12-Real5/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]
sddlpResult1210 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods12-Real10/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]

sddlpResult63.new_gap = parse.(Float64, replace.(sddlpResult63.gap, "%" => ""))
sddlpResult65.new_gap = parse.(Float64, replace.(sddlpResult65.gap, "%" => ""))
sddlpResult610.new_gap = parse.(Float64, replace.(sddlpResult610.gap, "%" => ""))
sddlpResult83.new_gap = parse.(Float64, replace.(sddlpResult83.gap, "%" => ""))
sddlpResult85.new_gap = parse.(Float64, replace.(sddlpResult85.gap, "%" => ""))
sddlpResult810.new_gap = parse.(Float64, replace.(sddlpResult810.gap, "%" => ""))
sddlpResult123.new_gap = parse.(Float64, replace.(sddlpResult123.gap, "%" => ""))
sddlpResult125.new_gap = parse.(Float64, replace.(sddlpResult125.gap, "%" => ""))
sddlpResult1210.new_gap = parse.(Float64, replace.(sddlpResult1210.gap, "%" => ""))




timegap = @df sddlpResult63 plot(:Time, :new_gap, label="T = 6, #R = 3", 
                                                # title = "Gap vs. Time", 
                                                xlab = "Time", 
                                                ylab = "Gap (%)", 
                                                xlim = [0,500],
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

@df sddlpResult65 plot!(:Time, :new_gap, marker=(:star, 2, 1.), label="T = 6, #R = 5", linestyle=:solid, color=:orange)
@df sddlpResult610 plot!(:Time, :new_gap, marker=(:hexagon, 2, 1.), label="T = 6, #R = 10", linestyle=:solid, color=:green)
@df sddlpResult83 plot!(:Time, :new_gap, marker=(:plus, 2, 1.), label="T = 8, #R = 3", linestyle=:solid, color=:Crimson)
@df sddlpResult85 plot!(:Time, :new_gap, marker=(:star, 2, 1.), label="T = 8, #R = 5", linestyle=:solid, color=:HotPink)
@df sddlpResult810 plot!(:Time, :new_gap, marker=(:hexagon, 2, 1.), label="T = 8, #R = 10", linestyle=:solid, color=:brown)
@df sddlpResult123 plot!(:Time, :new_gap, marker=(:xcross, 2, 1.), label="T = 12, #R = 3", linestyle=:solid, color=:lightslategray)
@df sddlpResult125 plot!(:Time, :new_gap, marker=(:square, 2, 1.), label="T = 12, #R = 5", linestyle=:solid, color=:cyan)
@df sddlpResult1210 plot!(:Time, :new_gap, marker=(:diamond, 2, 1.), label="T = 12, #R = 10", linestyle=:solid, color=:palevioletred)
