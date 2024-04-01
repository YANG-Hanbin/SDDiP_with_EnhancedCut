#############################################################################################
####################################    Without Parallel   ##################################
#############################################################################################
using Pkg
Pkg.activate(".")
using JuMP, Gurobi, PowerModels
using Statistics, StatsBase, Random, Dates, Distributions
using Distributed, ParallelDataTransfer
using CSV, DataFrames, Printf
using JLD2, FileIO

const GRB_ENV = Gurobi.Env()


include("src/UnitCommitment/Surrogate_SDDiP/def.jl")
include("src/UnitCommitment/Surrogate_SDDiP/backwardModel.jl");
include("src/UnitCommitment/Surrogate_SDDiP/forwardModel.jl");
include("src/UnitCommitment/Surrogate_SDDiP/LevelSetMethod.jl");
include("src/UnitCommitment/Surrogate_SDDiP/sddip.jl");
include("src/UnitCommitment/SDDiP/extForm.jl");
include("src/UnitCommitment/Surrogate_SDDiP/readin.jl");


#############################################################################################
####################################### Run Experiment ######################################
#############################################################################################
Output_Gap = false; max_iter = 200; cutSelection = "SMC"; δ = 30.; numScenarios = 30; tightness = true; ϵ = 1e-4; T = 3; num = 3;
case = "case30"; # "RTS_GMLC"
for cutSelection in ["LC", "ELC", "SMC"]
    for T in [3, 6, 8]
        for num in [3, 5, 10]
            indexSets = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/indexSets.jld2")["indexSets"]
            paramOPF = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/paramOPF.jld2")["paramOPF"]
            paramDemand = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/paramDemand.jld2")["paramDemand"]
            scenarioTree = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/scenarioTree.jld2")["scenarioTree"]
            initialStageDecision = load("src/UnitCommitment/experiment_$case/initialStageDecision.jld2")["initialStageDecision"]
            sddipResult = SDDiP_algorithm(scenarioTree = scenarioTree, 
                                indexSets = indexSets, 
                                    paramDemand = paramDemand, 
                                        paramOPF = paramOPF, 
                                            initialStageDecision = initialStageDecision, numScenarios = numScenarios, TimeLimit = 60 * 60 * 3.;
                                            Output_Gap = Output_Gap, max_iter = max_iter, δ = δ, cutSelection = cutSelection)
            save("src/UnitCommitment/experiment_$case/stage($T)real($num)/SsddipResult_5hr_$cutSelection.jld2", "sddipResult", sddipResult)
        
        end
    end
end