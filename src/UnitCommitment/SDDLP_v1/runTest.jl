#############################################################################################
####################################    Without Parallel   ##################################
#############################################################################################
using Pkg
Pkg.activate(".")
using Distributed; addprocs(5); 
@everywhere begin
    using JuMP, Gurobi, PowerModels;
    using Statistics, StatsBase, Random, Dates, Distributions;
    using Distributed, ParallelDataTransfer;
    using CSV, DataFrames, Printf;
    using JLD2, FileIO;

    const GRB_ENV = Gurobi.Env();


    include("src/UnitCommitment/SDDLP_v1/def.jl")
    include("src/UnitCommitment/SDDLP_v1/backwardModel.jl");
    include("src/UnitCommitment/SDDLP_v1/forwardModel.jl");
    include("src/UnitCommitment/SDDLP_v1/LevelSetMethod.jl");
    include("src/UnitCommitment/SDDLP_v1/sddip.jl");
    include("src/UnitCommitment/SDDLP_v1/extForm.jl");


    Output_Gap = false; max_iter = 150; MaxIter = 100; cutSelection = "LC"; δ = 1.; numScenarios = 100; tightness = true; case = "case30"; # "RTS_GMLC", "case30"
    T = 12; num = 3; TimeLimit = 60 * 60 * 2.; OPT = Inf;
end


for cutSelection in ["LC", "ELC", "SMC"]
    for num in [3, 5, 10]
        for T in [6, 8, 12] 
            @everywhere begin
                T = $T; num = $num; cutSelection = $cutSelection;
                indexSets = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/indexSets.jld2")["indexSets"]
                paramOPF = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/paramOPF.jld2")["paramOPF"]
                paramDemand = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/paramDemand.jld2")["paramDemand"]
                scenarioTree = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/scenarioTree.jld2")["scenarioTree"]
                initialStageDecision = load("src/UnitCommitment/experiment_$case/initialStageDecision.jld2")["initialStageDecision"]
            end
            # Ξ = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/Ξ.jld2")["Ξ"]
            # @time extResult = extensive_form(indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF, scenarioTree = scenarioTree, Ξ = Ξ, initialStageDecision = initialStageDecision); OPT = extResult.OPT;
            sddlpResult = SDDiP_algorithm(scenarioTree = scenarioTree, 
                                indexSets = indexSets, 
                                    paramDemand = paramDemand, 
                                        paramOPF = paramOPF, 
                                            initialStageDecision = initialStageDecision, numScenarios = numScenarios, TimeLimit = TimeLimit, OPT = OPT, tightness = tightness,
                                            Output_Gap = Output_Gap, max_iter = max_iter, MaxIter = MaxIter, δ = δ, cutSelection = cutSelection)
            save("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddlpResult-$cutSelection-$tightness.jld2", "sddlpResult", sddlpResult)
        end
    end
end

T = 8; num = 5; cutSelection = "SMC"; sddlpResult = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddlpResult-$cutSelection-$tightness.jld2")["sddlpResult"][:solHistory]
