using Pkg
Pkg.activate(".")
using JuMP, Gurobi, PowerModels
using Statistics, StatsBase, Random, Dates, Distributions
using Distributed, ParallelDataTransfer
using CSV, DataFrames, Printf
using JLD2, FileIO
using StatsPlots



include("src/UnitCommitment/SDDiP/def.jl")



cutSelection = "ELC" # for cutSelection in ["LC", "ELC", "SMC"]
T = 6; num= 5;
case = "RTS_GMLC";
# sddipResult = load("src/UnitCommitment/experiment/stage($T)real($num)/sddipResult_5hr_LC.jld2")["sddipResult"]
sddipResult = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/SsddipResult_5hr_$cutSelection.jld2")["sddipResult"]
# sddipResult[:solHistory][1:30,:]
sddipResult[:solHistory]
describe(sddipResult[:solHistory])


p = @df sddipResult = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/SsddipResult_5hr_LC.jld2")["sddipResult"][:solHistory] plot(:iter, :LB, label="AS-LC", 
                                                title = "Convergence v.s. Iteration", 
                                                xlab = "Iterations", 
                                                ylab = "Lower Bounds",
                                                # ylim = [1000,10000],
                                                # xlim = [0,200],
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
@df load("src/UnitCommitment/experiment_$case/stage($T)real($num)/SsddipResult_5hr_ELC.jld2")["sddipResult"][:solHistory][:, :] plot!(:iter, :LB, label="AS-ELC")
@df load("src/UnitCommitment/experiment_$case/stage($T)real($num)/SsddipResult_5hr_SMC.jld2")["sddipResult"][:solHistory][:, :] plot!(:iter, :LB, label="AS-SMC")
@df load("src/UnitCommitment/experiment_$case/stage($T)real($num)/sddipResult_5hr_LC.jld2")["sddipResult"][:solHistory][:, :] plot!(:iter, :LB, label="LC")
p |> save("src/UnitCommitment/experiment_$case/stage($T)real($num)/convergence_iter.pdf")


p = @df sddipResult = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/SsddipResult_5hr_LC.jld2")["sddipResult"][:solHistory] plot(:Time, :LB, label="AS-LC", 
                                                title = "Convergence v.s. Time", 
                                                xlab = "Time", 
                                                ylab = "Lower Bounds",
                                                # ylim = [1000,10000],
                                                # xlim = [0,200],
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
@df load("src/UnitCommitment/experiment_$case/stage($T)real($num)/SsddipResult_5hr_ELC.jld2")["sddipResult"][:solHistory][:, :] plot!(:Time, :LB, label="AS-ELC")
@df load("src/UnitCommitment/experiment_$case/stage($T)real($num)/SsddipResult_5hr_SMC.jld2")["sddipResult"][:solHistory][:, :] plot!(:Time, :LB, label="AS-SMC")
@df load("src/UnitCommitment/experiment_$case/stage($T)real($num)/sddipResult_5hr_LC.jld2")["sddipResult"][:solHistory][:, :] plot!(:Time, :LB, label="LC")
p |> save("src/UnitCommitment/experiment_$case/stage($T)real($num)/convergence_time.pdf")