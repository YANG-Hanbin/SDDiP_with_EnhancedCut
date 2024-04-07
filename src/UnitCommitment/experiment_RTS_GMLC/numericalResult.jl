using Pkg
Pkg.activate(".")
using JuMP, Gurobi, PowerModels
using Statistics, StatsBase, Random, Dates, Distributions
using Distributed, ParallelDataTransfer
using CSV, DataFrames, Printf
using JLD2, FileIO
using StatsPlots



include("src/UnitCommitment/SDDiP/def.jl")



cutSelection = "SMC" # for cutSelection in ["LC", "ELC", "SMC"]
T = 8; num = 3;
case = "RTS_GMLC";
# sddipResult = load("src/UnitCommitment/experiment/stage($T)real($num)/sddipResult_5hr_LC.jld2")["sddipResult"]
sddipResult = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/SsddipResult_5hr_$cutSelection.jld2")["sddipResult"]
# sddipResult[:solHistory][1:30,:]
sddipResult[:solHistory]
describe(sddipResult[:solHistory])


p = @df sddipResult = load("src/UnitCommitment/experiment/stage($T)real($num)/SsddipResult_5hr_LC.jld2")["sddipResult"][:solHistory] plot(:iter, :LB, label="ASLC", 
                                                title = "Convergence v.s. Iteration", 
                                                xlab = "Iterations", 
                                                ylab = "Lower and Upper Bounds",
                                                # ylim = [1000,10000],
                                                xlim = [0,200],
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
@df load("src/UnitCommitment/experiment/stage($T)real($num)/sddipResult_5hr_LC.jld2")["sddipResult"][:solHistory][:, :] plot!(:Time, :LB, label="LC")



p |> save("src/experiments/cutPerformance/convergence_time.pdf")