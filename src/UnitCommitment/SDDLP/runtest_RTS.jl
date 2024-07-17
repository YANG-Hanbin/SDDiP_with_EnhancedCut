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


include("src/UnitCommitment/SDDLP/def.jl")
include("src/UnitCommitment/SDDLP/backwardModel.jl");
include("src/UnitCommitment/SDDLP/forwardModel.jl");
include("src/UnitCommitment/SDDLP/LevelSetMethod.jl");
include("src/UnitCommitment/SDDLP/sddip.jl");
include("src/UnitCommitment/SDDiP/extForm.jl");
include("src/UnitCommitment/SDDLP/readin.jl");

#############################################################################################
####################################### Run Experiment ######################################
#############################################################################################
Output_Gap = false; TimeLimit = 18000; max_iter = 100; MaxIter = 20; cutSelection = "SMC"; δ = 1.; numScenarios = 30; tightness = false; ϵ = 1e-4; case = "RTS_GMLC"; # "RTS_GMLC", "case30"
T = 3; num = 3; TimeLimit = 60 * 60 * 5.; OPT = Inf; seed = nothing
for cutSelection in ["LC", "ELC", "SMC"]
    for T in [3, 6, 8]
        for num in [3, 5, 10]
            for seed in 1:10
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
                                                initialStageDecision = initialStageDecision, numScenarios = numScenarios, TimeLimit = TimeLimit, OPT = OPT, tightness = tightness,
                                                Output_Gap = Output_Gap, max_iter = max_iter, MaxIter = MaxIter, δ = δ, cutSelection = cutSelection, seed = seed)
                save("src/UnitCommitment/experiment_data_$case/stage($T)real($num)/sddlp_tight($tightness)_seed($seed)_$cutSelection.jld2", "sddipResult", sddipResult)
            end
        end
    end
end

cutSelection = "SMC"; tightness = false; case = "RTS_GMLC"; # "RTS_GMLC", "case30"
T = 3; num = 3; seed = 1;
sddipResult = load("src/UnitCommitment/experiment_data_$case/stage($T)real($num)/sddlp_tight($tightness)_seed($seed)_$cutSelection.jld2")["sddipResult"]
sddipResult[:solHistory]
