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


include("src/UnitCommitment/SDDiP/def.jl")
include("src/UnitCommitment/SDDiP/backwardModel.jl");
include("src/UnitCommitment/SDDiP/forwardModel.jl");
include("src/UnitCommitment/SDDiP/LevelSetMethod.jl");
include("src/UnitCommitment/SDDiP/sddip.jl");
include("src/UnitCommitment/SDDiP/extForm.jl");
include("src/UnitCommitment/SDDiP/readin.jl");



# indexSets = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/indexSets.jld2")["indexSets"]
# paramOPF = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/paramOPF.jld2")["paramOPF"]
# paramDemand = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/paramDemand.jld2")["paramDemand"]
# scenarioTree = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/scenarioTree.jld2")["scenarioTree"]
# initialStageDecision = load("src/UnitCommitment/experiment_$case/initialStageDecision.jld2")["initialStageDecision"]
# Ξ = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/Ξ.jld2")["Ξ"]
#############################################################################################
####################################### Run Experiment ######################################
#############################################################################################
Output_Gap = false; max_iter = 150; cutSelection = "ELC"; δ = 10.; numScenarios = 30; tightness = true; ϵ = 1e-4; T = 3; num = 3; OPT = Inf; ε = 1/2^(3);
case = "RTS_GMLC"; # "RTS_GMLC"
for cutSelection in ["LC", "ELC", "SMC"]
    for T in [3, 6, 8]
        for num in [3, 5, 10]
            indexSets = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/indexSets.jld2")["indexSets"]
            paramOPF = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/paramOPF.jld2")["paramOPF"]
            paramDemand = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/paramDemand.jld2")["paramDemand"]
            scenarioTree = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/scenarioTree.jld2")["scenarioTree"]
            initialStageDecision = load("src/UnitCommitment/experiment_$case/initialStageDecision.jld2")["initialStageDecision"]
            # Ξ = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/Ξ.jld2")["Ξ"]
            # @time extResult = extensive_form(indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF, scenarioTree = scenarioTree, Ξ = Ξ, initialStageDecision = initialStageDecision); OPT = extResult.OPT;
            sddipResult = SDDiP_algorithm(scenarioTree = scenarioTree, 
                                indexSets = indexSets, 
                                    paramDemand = paramDemand, 
                                        paramOPF = paramOPF, 
                                            initialStageDecision = initialStageDecision, numScenarios = numScenarios, TimeLimit = 60 * 60 * 5., OPT = OPT,
                                            Output_Gap = Output_Gap, max_iter = max_iter, ε = ε, δ = δ, cutSelection = cutSelection)
            save("src/UnitCommitment/experiment_$case/stage($T)real($num)/sddipResult_5hr_$cutSelection.jld2", "sddipResult", sddipResult)
        
        end
    end
end
T = 8; num = 10; cutSelection = "SMC";
sddipResult = load("src/UnitCommitment/experiment_RTS_GMLC/stage($T)real($num)/sddipResult_5hr_$cutSelection.jld2")["sddipResult"]
# sddipResult = load("src/UnitCommitment/experiment_RTS_GMLC/stage($T)real($num)/SsddipResult_5hr_$cutSelection(_false).jld2")["sddipResult"]
sddipResult[:solHistory]
describe(sddipResult[:solHistory])