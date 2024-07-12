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


include("src/UnitCommitment/SDDP/def.jl")
include("src/UnitCommitment/SDDP/backwardModel.jl");
include("src/UnitCommitment/SDDP/forwardModel.jl");
include("src/UnitCommitment/SDDP/LevelSetMethod.jl");
include("src/UnitCommitment/SDDP/sddp.jl");
include("src/UnitCommitment/SDDP/extForm.jl");
include("src/UnitCommitment/SDDP/readin.jl");

#############################################################################################
####################################### Run Experiment ######################################
#############################################################################################
Output_Gap = false; max_iter = 150; cutSelection = "LC"; δ = 10.; numScenarios = 30; tightness = true; ϵ = 1e-4; T = 3; num = 3; OPT = Inf;
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
            sddpResult = SDDP_algorithm(scenarioTree = scenarioTree, 
                                indexSets = indexSets, 
                                    paramDemand = paramDemand, 
                                        paramOPF = paramOPF, 
                                            initialStageDecision = initialStageDecision, numScenarios = numScenarios, TimeLimit = 60 * 60 * 5., OPT = OPT,
                                            Output_Gap = Output_Gap, max_iter = max_iter, δ = δ, cutSelection = cutSelection)
            save("src/UnitCommitment/experiment_$case/stage($T)real($num)/sddpResult_5hr_$cutSelection.jld2", "sddpResult", sddpResult)
        
        end
    end
end