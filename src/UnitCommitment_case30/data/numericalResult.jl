using Pkg
Pkg.activate(".")
using JuMP, Gurobi, PowerModels
using Statistics, StatsBase, Random, Dates, Distributions
using Distributed, ParallelDataTransfer
using CSV, DataFrames, Printf
using JLD2, FileIO
using StatsPlots, PlotThemes

project_root = @__DIR__;
include(joinpath(project_root, "src", "UnitCommitment_case30", "SDDiP", "def.jl"))




case = "case30"; tightness = true; 
T = 6; num = 5;ℓ = 0.0;
theme(:default)
## ================================================================ the same algorithm with different cuts ================================================================ ##
# LB vs Iter: 
lbiter = @df load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddipResult-SMC-$tightness-5.jld2")["sddipResult"][:solHistory] plot(:Iter, :LB, 
                                                label="SMC", 
                                                title = "Lower Bounds vs. Iteration", 
                                                xlab = "Iteration", 
                                                # ylab = "Lower Bounds",
                                                # ylim = [1000,10000],
                                                xlim = [0,60],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                size=(400,300),
                                                marker=(:circle, 3, 1.),
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"), legend=:right) # :right, :left, :top, :bottom, :inside, :best, :legend, :topright, :topleft, :bottomleft, :bottomright
@df load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddipResult-LC-$tightness-5.jld2")["sddipResult"][:solHistory] plot!(:Iter, :LB, marker=(:star, 3, 1.), label="LC")
@df load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/MagnantiWong/sddipResult-$ℓ-$tightness.jld2")["sddpResult"][:solHistory] plot!(:Iter, :LB, marker=(:hexagon, 3, 1.), label="ELC")
lbiter |> save("/Users/aaron/Downloads/figures/sddipResult-cut_comparsion.pdf")


lbiter = @df load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddpResult-SMC-$tightness.jld2")["sddpResult"][:solHistory] plot(:Iter, :LB, 
                                                label="SMC", 
                                                title = "Lower Bounds vs. Iteration", 
                                                xlab = "Iteration", 
                                                # ylab = "Lower Bounds",
                                                # ylim = [1000,10000],
                                                xlim = [0,60],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                size=(400,300),
                                                marker=(:circle, 3, 1.),
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"), legend=:right) # :right, :left, :top, :bottom, :inside, :best, :legend, :topright, :topleft, :bottomleft, :bottomright
@df load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddpResult-LC-$tightness.jld2")["sddpResult"][:solHistory] plot!(:Iter, :LB, marker=(:star, 3, 1.), label="LC")
@df load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/MagnantiWong/sddpResult-$ℓ-$tightness.jld2")["sddpResult"][:solHistory] plot!(:Iter, :LB, marker=(:hexagon, 3, 1.), label="ELC")
lbiter |> save("/Users/aaron/Downloads/figures/sddpResult-cut_comparsion.pdf")

lbiter = @df load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddlpResult-SMC-$tightness.jld2")["sddlpResult"][:solHistory] plot(:Iter, :LB, 
                                                label="SMC", 
                                                title = "Lower Bounds vs. Iteration", 
                                                xlab = "Iteration", 
                                                # ylab = "Lower Bounds",
                                                # ylim = [1000,10000],
                                                xlim = [0,60],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                size=(400,300),
                                                marker=(:star, 3, 1.),
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"), legend=:right) # :right, :left, :top, :bottom, :inside, :best, :legend, :topright, :topleft, :bottomleft, :bottomright
@df load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddlpResult-LC-$tightness.jld2")["sddlpResult"][:solHistory] plot!(:Iter, :LB, marker=(:hexagon, 3, 1.), label="LC")

@df load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/MagnantiWong/sddlpResult-$ℓ-$tightness.jld2")["sddlpResult"][:solHistory] plot!(:Iter, :LB, marker=(:circle, 3, 1.), label="PLC")
lbiter |> save("/Users/aaron/Downloads/figures/sddlp-cut_comparsion.pdf")


# LB vs Time: the same algorithm with different cuts
lbtime = @df load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddlpResult-SMC-$tightness.jld2")["sddlpResult"][:solHistory] plot(:Time, :LB, 
                                                label="SMC", 
                                                title = "Lower Bounds vs. Time", 
                                                xlab = "Time (sec.)", 
                                                # ylab = "Lower Bounds",
                                                # ylim = [1000,10000],
                                                # xlim = [0,100],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman")
                                    )


@df load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddlpResult-LC-$tightness.jld2")["sddlpResult"][:solHistory] plot!(:Time, :LB, label="LC")

@df load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/MagnantiWong/sddlpResult-$ℓ-$tightness.jld2")["sddlpResult"][:solHistory] plot!(:Time, :LB, label="ELC")



## ================================================================ different algorithms with the same cuts ================================================================ ##
# LB vs Iter: 
cutSelection = "SMC"; 
isddlpResult = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddpResult-$cutSelection-$tightness.jld2")["sddpResult"][:solHistory]
sddlpResult = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddlpResult-$cutSelection-$tightness.jld2")["sddlpResult"][:solHistory]
sddpResult = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/isddlpResult-$cutSelection-$tightness.jld2")["sddlpResult"][:solHistory]
lbiter = @df sddpResult plot(:Iter, :LB, label="SDDP", 
                                                title = "Lower bounds vs. Iteration", 
                                                xlab = "Iteration", 
                                                # ylab = "Lower Bounds",
                                                # ylim = [1000,10000],
                                                xlim = [0,60],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                size=(400,300),
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),marker=(:circle, 3, 1.),
                                                legendfont=font("Times New Roman"));

@df sddlpResult plot!(:Iter, :LB, marker=(:star, 3, 1.), label="SDDℓP");
@df isddlpResult plot!(:Iter, :LB, marker=(:hexagon, 3, 1.), label="iSDDℓP")

lbiter |> save("/Users/aaron/Downloads/figures/SMC-alg_comparsion_iter.pdf")

# LB vs Iter: different algorithms with the same cuts
lbtime = @df sddpResult plot(:Time, :LB, label="SDDP", 
                                                title = "Lower bounds vs. Time", 
                                                xlab = "Time (sec.)", 
                                                # ylab = "Lower Bounds",
                                                # ylim = [20000,23000],
                                                xlim = [0,600],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                size=(400,300),
                                                # yaxis=(formatter=y->string(round(Int, y / 10^3)) * "k"),
                                                # yformatter = :scientific,
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"), marker=(:circle, 2, 1.),
                                                legendfont=font("Times New Roman"))
@df sddlpResult plot!(:Time, :LB, marker=(:star, 2, 1.), label="SDDℓP")
@df isddlpResult plot!(:Time, :LB, marker=(:hexagon, 2, 1.), label="iSDDℓP")
lbtime |> save("/Users/aaron/Downloads/figures/SMC-alg_comparsion_time.pdf")


## ================================================================ accuracy of binarization ================================================================ ##
coef = 10; T = 6; num = 5;
sddipResult5 = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-5.jld2")["sddipResult"][:solHistory]
sddipResult6 = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-6.jld2")["sddipResult"][:solHistory]
sddipResult7 = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-7.jld2")["sddipResult"][:solHistory]
sddipResult8 = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-8.jld2")["sddipResult"][:solHistory]
sddipResult9 = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-9.jld2")["sddipResult"][:solHistory]
sddipResult10 = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-10.jld2")["sddipResult"][:solHistory]
# isddlpResult = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/isddlpResult-$cutSelection-$tightness.jld2")["sddlpResult"][:solHistory]

lbiter = @df sddipResult6 plot(:Iter, [:LB, :UB], label="1/64", 
                                                title = "Lower bounds vs. Iteration", 
                                                xlab = "Iteration", 
                                                # ylab = "Lower Bounds",
                                                # ylim = [7000,22000],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=10, 
                                                size=(400,300),
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"), marker=(:circle, 2, 1.), 
                                                linewidth=1,
                                                legend=:bottomright); # outertopright
@df sddipResult7 plot!(:Iter, [:LB, :UB], linewidth=1, marker=(:star, 2, 1.), label="1/128")
@df sddipResult8 plot!(:Iter, [:LB, :UB], linewidth=1, marker=(:hexagon, 2, 1.), label="1/256")
@df sddipResult9 plot!(:Iter, [:LB, :UB], linewidth=1, marker=(:star, 2, 1.), label="1/512")
@df sddipResult10 plot!(:Iter, [:LB, :UB], linewidth=1, marker=(:circle, 2, 1.), label="1/1024")


lbiter = @df sddipResult6 plot(:Iter, :LB, label="1/64", 
                                                title = "Upper and lower bounds vs. Iteration", 
                                                xlab = "Iteration", 
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=10, 
                                                # size=(400,300),
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"), 
                                                marker=(:circle, 2, 1.), 
                                                linewidth=1,
                                                legend=:bottomright,
                                                color=:blue)  # Use blue color for 1/64
@df sddipResult6 plot!(:Iter, :UB, label="", linewidth=1, marker=(:circle, 2, 1.), color=:blue)  # Plot UB for 1/64, no label

# Plot for 1/128 (LB and UB with the same color)
@df sddipResult7 plot!(:Iter, :LB, label="1/128", linewidth=1, marker=(:star, 2, 1.), color=:red)  # Use red color for 1/128
@df sddipResult7 plot!(:Iter, :UB, label="", linewidth=1, marker=(:star, 2, 1.), color=:red)  # Plot UB for 1/128, no label

@df sddipResult8 plot!(:Iter, :LB, linewidth=1, marker=(:hexagon, 2, 1.), label="1/256", color=:orange)
@df sddipResult8 plot!(:Iter, :UB, linewidth=1, marker=(:hexagon, 2, 1.), label="", color=:orange)

@df sddipResult9 plot!(:Iter, :LB, linewidth=1, marker=(:star, 2, 1.), label="1/512", color=:green)
@df sddipResult9 plot!(:Iter, :UB, linewidth=1, marker=(:hexagon, 2, 1.), label="", color=:green)

@df sddipResult10 plot!(:Iter, :LB, linewidth=1, marker=(:circle, 2, 1.), label="1/1024", color=:purple)
@df sddipResult10 plot!(:Iter, :UB, linewidth=1, marker=(:hexagon, 2, 1.), label="", color=:purple)


lbiter |> save("/Users/aaron/Downloads/figures/SMC-sddip_binarization_comparsion_iter.pdf")

# @df isddlpResult plot!(:Iter, :LB, linewidth=1, marker=(:hexagon, 2, 1.), label="Benchmark")


lbtime = @df sddipResult6 plot(:Time, :LB, label="1/64", 
                                                title = "Upper and lower bounds vs. Time", 
                                                xlab = "Time (sec.)", 
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=10, 
                                                # size=(400,300),
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"), 
                                                marker=(:circle, 2, 1.), 
                                                linewidth=1,
                                                legend=:bottomright,
                                                color=:blue)  # Use blue color for 1/64
@df sddipResult6 plot!(:Time, :UB, label="", linewidth=1, marker=(:circle, 2, 1.), color=:blue)  # Plot UB for 1/64, no label

# Plot for 1/128 (LB and UB with the same color)
@df sddipResult7 plot!(:Time, :LB, label="1/128", linewidth=1, marker=(:star, 2, 1.), color=:red)  # Use red color for 1/128
@df sddipResult7 plot!(:Time, :UB, label="", linewidth=1, marker=(:star, 2, 1.), color=:red)  # Plot UB for 1/128, no label

@df sddipResult8 plot!(:Time, :LB, linewidth=1, marker=(:hexagon, 2, 1.), label="1/256", color=:orange)
@df sddipResult8 plot!(:Time, :UB, linewidth=1, marker=(:hexagon, 2, 1.), label="", color=:orange)

@df sddipResult9 plot!(:Time, :LB, linewidth=1, marker=(:star, 2, 1.), label="1/512", color=:green)
@df sddipResult9 plot!(:Time, :UB, linewidth=1, marker=(:hexagon, 2, 1.), label="", color=:green)

@df sddipResult10 plot!(:Time, :LB, linewidth=1, marker=(:circle, 2, 1.), label="1/1024", color=:purple)
@df sddipResult10 plot!(:Time, :UB, linewidth=1, marker=(:hexagon, 2, 1.), label="", color=:purple)
lbtime |> save("/Users/aaron/Downloads/figures/SMC-sddip_binarization_comparsion_time.pdf")

## ================================================================ Cut difference ================================================================ ##
T = 6; num = 5; tightness = true; cutSelection = "SMC";
sddlpResultLC = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddlpResult-LC-$tightness.jld2")["sddlpResult"][:solHistory]
sddlpResultSMC = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddlpResult-SMC-$tightness.jld2")["sddlpResult"][:solHistory]
sddlpResultELC = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/MagnantiWong/sddlpResult-0.0-$tightness.jld2")["sddlpResult"][:solHistory]
sddpResultSMC = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddpResult-SMC-$tightness.jld2")["sddpResult"][:solHistory]
sddpResultLC = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddpResult-LC-$tightness.jld2")["sddpResult"][:solHistory]
isddlpResult = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/isddlpResult-SMC-$tightness.jld2")["sddlpResult"][:solHistory]

timeiter = @df sddlpResultLC plot(:Iter, :time, label="SDDℓP-LC", 
                                                title = "Iteration time (sec.) vs. Iteration", 
                                                xlab = "Iteration", 
                                                # ylab = "Iteration time",
                                                # ylim = [20000,23000],
                                                # xlim = [0,6000],
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
                                                legendfont=font("Times New Roman"), legend=:topleft);
@df sddlpResultSMC plot!(:Iter, :time, marker=(:star, 2, 1.), label="SDDℓP-SMC")
@df sddlpResultELC plot!(:Iter, :time, marker=(:hexagon, 2, 1.), label="SDDℓP-ELC")
@df sddpResultSMC plot!(:Iter, :time, marker=(:circle, 2, 1.), label="SDDP-SMC")
@df sddpResultLC plot!(:Iter, :time, marker=(:circle, 2, 1.), label="SDDP-LC")
@df isddlpResult plot!(:Iter, :time, marker=(:star, 2, 1.), label="iSDDℓP-SMC")
@df sddipResult7 plot!(:Iter, :time, marker=(:hexagon, 2, 1.), label="SDDiP-128")

timeiter |> save("/Users/aaron/Downloads/figures/cut_iteration_time.pdf")



## ================================================================ ELC Exploration ================================================================ ##
T = 6; num = 3; tightness = true; cutSelection = "SMC";
sddlpResultELC0 = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/MagnantiWong/sddlpResult-0.0-$tightness.jld2")["sddlpResult"][:solHistory]
sddlpResultELC5 = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/MagnantiWong/sddlpResult-0.5-$tightness.jld2")["sddlpResult"][:solHistory]
sddlpResultELC9 = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/MagnantiWong/sddlpResult-0.9-$tightness.jld2")["sddlpResult"][:solHistory]

lbiter = @df sddlpResultELC0[3:100,:] plot(:Iter, :LB, label="ELC-0.0", 
                                                title = "Lower bounds vs. Iteration", 
                                                xlab = "Iteration", 
                                                # ylab = "Iteration time",
                                                # ylim = [20000,23000],
                                                xlim = [3,100],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                size=(400,300),
                                                marker=(:circle, 2, 1.), 
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"));
@df sddlpResultELC5[3:100,:] plot!(:Iter, :LB, marker=(:star, 2, 1.), label="ELC-0.5")
@df sddlpResultELC9[3:100,:] plot!(:Iter, :LB, marker=(:hexagon, 2, 1.), label="ELC-0.9")
lbiter |> save("/Users/aaron/Downloads/figures/ELC_lb_iter.pdf")

lbtime = @df sddlpResultELC0 plot(:Time, :LB, label="ELC-0.0", 
                                                title = "Lower bounds vs. Time", 
                                                xlab = "Time (sec.)", 
                                                # ylab = "Iteration time",
                                                # ylim = [20000,23000],
                                                # xlim = [3,100],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                size=(400,300),
                                                marker=(:circle, 2, 1.), 
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"));
@df sddlpResultELC5 plot!(:Time, :LB, marker=(:star, 2, 1.), label="ELC-0.5")
@df sddlpResultELC9 plot!(:Time, :LB, marker=(:hexagon, 2, 1.), label="ELC-0.9")
lbtime |> save("/Users/aaron/Downloads/figures/ELC_lb_time.pdf")


timeiter = @df sddlpResultELC0[1:100,:] plot(:Iter, :time, label="ELC-0.0", 
                                                title = "Iteration time (sec.) vs. Iteration", 
                                                xlab = "Iteration", 
                                                # ylab = "Iteration time",
                                                # ylim = [20000,23000],
                                                xlim = [3,100],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                size=(400,300),
                                                marker=(:circle, 2, 1.), 
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"));
@df sddlpResultELC5[1:100,:] plot!(:Iter, :time, marker=(:star, 2, 1.), label="ELC-0.5")
@df sddlpResultELC9[1:100,:] plot!(:Iter, :time, marker=(:hexagon, 2, 1.), label="ELC-0.9")
timeiter |> save("/Users/aaron/Downloads/figures/ELC_iteration_time.pdf")





## ================================================================ Overall Performance ================================================================ ##
T = 6; num = 5; tightness = true; cutSelection = "SMC";
sddlpResultLC = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddlpResult-LC-$tightness.jld2")["sddlpResult"][:solHistory]
sddlpResultELC = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/MagnantiWong/sddlpResult-0.0-$tightness.jld2")["sddlpResult"][:solHistory]
sddlpResultSMC = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddlpResult-SMC-$tightness.jld2")["sddlpResult"][:solHistory]

sddpResultLC = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddpResult-LC-$tightness.jld2")["sddpResult"][:solHistory]
sddpResultELC = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/MagnantiWong/sddpResult-0.0-$tightness.jld2")["sddpResult"][:solHistory]
sddpResultSMC = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddpResult-SMC-$tightness.jld2")["sddpResult"][:solHistory]

sddipResultLC = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddipResult-LC-$tightness-7.jld2")["sddipResult"][:solHistory][1:157,:]
sddipResultELC = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/MagnantiWong/sddipResult-$ℓ-$tightness-7.jld2")["sddipResult"][:solHistory]
sddipResultSMC = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddipResult-SMC-$tightness-7.jld2")["sddipResult"][:solHistory]

timeiter = @df sddlpResultLC plot(:Iter, :time, label="SDDℓP-LC", 
                                                title = "Iteration time (sec.) vs. Iteration", 
                                                xlab = "Iteration", 
                                                # ylab = "Iteration time",
                                                ylim = [-10,350],
                                                xlim = [0,100],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                size=(700,430),
                                                marker=(:circle, 1, 1.), 
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"), legend=:outerright)
@df sddlpResultELC plot!(:Iter, :time, marker=(:hexagon, 1, 1.), label="SDDℓP-PLC")
@df sddlpResultSMC plot!(:Iter, :time, marker=(:star, 1, 1.), label="SDDℓP-SMC")
@df sddpResultLC plot!(:Iter, :time, marker=(:circle, 1, 1.), label="SDDP-LC")
@df sddpResultELC plot!(:Iter, :time, marker=(:hexagon, 1, 1.), label="SDDP-PLC")
@df sddpResultSMC plot!(:Iter, :time, marker=(:star, 1, 1.), label="SDDP-SMC")
@df sddipResultLC plot!(:Iter, :time, marker=(:circle, 1, 1.), label="SDDiP-LC")
@df sddipResultELC plot!(:Iter, :time, marker=(:hexagon, 1, 1.), label="SDDiP-PLC")
@df sddipResultSMC plot!(:Iter, :time, marker=(:star, 1, 1.), label="SDDiP-SMC")
timeiter |> save("/Users/aaron/Downloads/cut_iteration_time.pdf")



lbtime = @df sddlpResultLC plot(:Iter, :LB, label="SDDℓP-LC", 
                                                title = "Lower bounds vs. Time", 
                                                xlab = "Time", 
                                                # ylab = "Iteration time",
                                                # ylim = [20000,23000],
                                                xlim = [0,100],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                size=(700,430),
                                                marker=(:circle, 1, 1.), 
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"), legend=:right)
@df sddlpResultELC plot!(:Iter, :LB, marker=(:hexagon, 1, 1.), label="SDDℓP-PLC")
@df sddlpResultSMC plot!(:Iter, :LB, marker=(:star, 1, 1.), label="SDDℓP-SMC")
@df sddpResultLC plot!(:Iter, :LB, marker=(:circle, 1, 1.), label="SDDP-LC")
@df sddpResultELC plot!(:Iter, :LB, marker=(:hexagon, 1, 1.), label="SDDP-PLC")
@df sddpResultSMC plot!(:Iter, :LB, marker=(:star, 1, 1.), label="SDDP-SMC")

# plot!(legend=false)



lbiter = @df sddipResult6 plot(:Iter, :LB, label="1/64", 
                                                title = "Upper and lower bounds vs. Iteration", 
                                                xlab = "Iteration", 
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=10, 
                                                # size=(400,300),
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"), 
                                                marker=(:circle, 2, 1.), 
                                                linewidth=1,
                                                legend=:bottomright,
                                                color=:blue)  # Use blue color for 1/64
@df sddipResult6 plot!(:Iter, :UB, label="", linewidth=1, marker=(:circle, 2, 1.), color=:blue)  # Plot UB for 1/64, no label

# Plot for 1/128 (LB and UB with the same color)
@df sddipResult7 plot!(:Iter, :LB, label="1/128", linewidth=1, marker=(:star, 2, 1.), color=:red)  # Use red color for 1/128
@df sddipResult7 plot!(:Iter, :UB, label="", linewidth=1, marker=(:star, 2, 1.), color=:red)  # Plot UB for 1/128, no label

@df sddipResult8 plot!(:Iter, :LB, linewidth=1, marker=(:hexagon, 2, 1.), label="1/256", color=:orange)
@df sddipResult8 plot!(:Iter, :UB, linewidth=1, marker=(:hexagon, 2, 1.), label="", color=:orange)

@df sddipResult9 plot!(:Iter, :LB, linewidth=1, marker=(:star, 2, 1.), label="1/512", color=:green)
@df sddipResult9 plot!(:Iter, :UB, linewidth=1, marker=(:hexagon, 2, 1.), label="", color=:green)

@df sddipResult10 plot!(:Iter, :LB, linewidth=1, marker=(:circle, 2, 1.), label="1/1024", color=:purple)
@df sddipResult10 plot!(:Iter, :UB, linewidth=1, marker=(:hexagon, 2, 1.), label="", color=:purple)