using Pkg
Pkg.activate(".")
using JuMP, Gurobi, PowerModels
using Statistics, StatsBase, Random, Dates, Distributions
using Distributed, ParallelDataTransfer
using CSV, DataFrames, Printf
using JLD2, FileIO
using StatsPlots



include("src/UnitCommitment/SDDiP/def.jl")



for T in [3, 6, 8]
    for num in [3, 5, 10]
        indexSets = load("src/UnitCommitment/experiment/stage($T)real($num)/indexSets.jld2")["indexSets"]
        paramOPF = load("src/UnitCommitment/experiment/stage($T)real($num)/paramOPF.jld2")["paramOPF"]
        paramDemand = load("src/UnitCommitment/experiment/stage($T)real($num)/paramDemand.jld2")["paramDemand"]
        scenarioTree = load("src/UnitCommitment/experiment/stage($T)real($num)/scenarioTree.jld2")["scenarioTree"]
        initialStageDecision = load("src/UnitCommitment/experiment/stage(3)real(3)/initialStageDecision.jld2")["initialStageDecision"]
        sddipResult = SDDiP_algorithm(scenarioTree = scenarioTree, 
                            indexSets = indexSets, 
                                paramDemand = paramDemand, 
                                    paramOPF = paramOPF, 
                                        initialStageDecision = initialStageDecision,
                                        Output_Gap = Output_Gap, max_iter = max_iter, ϵ = ϵ, δ = δ, cutSelection = cutSelection)
        save("src/UnitCommitment/experiment/stage($T)real($num)/sddipResult_5hr_$cutSelection.jld2", "sddipResult", sddipResult)
    
    end
end


cutSelection = "LC" # for cutSelection in ["LC", "ELC", "SMC"]
T = 6; num = 3;
# sddipResult = load("src/UnitCommitment/experiment/stage($T)real($num)/sddipResult_5hr_LC.jld2")["sddipResult"]
sddipResult = load("src/UnitCommitment/experiment/stage($T)real($num)/SsddipResult_5hr_$cutSelection.jld2")["sddipResult"]
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