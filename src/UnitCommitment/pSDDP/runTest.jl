using Pkg
Pkg.activate(".")
using Distributed; addprocs(5); 
@everywhere begin
    using JuMP, Gurobi, PowerModels
    using Statistics, StatsBase, Random, Dates, Distributions
    using Distributed, ParallelDataTransfer
    using CSV, DataFrames, Printf
    using JLD2, FileIO

    const GRB_ENV = Gurobi.Env()


    include("src/UnitCommitment/pSDDP/def.jl")
    include("src/UnitCommitment/pSDDP/backwardModel.jl");
    include("src/UnitCommitment/pSDDP/forwardModel.jl");
    include("src/UnitCommitment/pSDDP/LevelSetMethod.jl");
    include("src/UnitCommitment/pSDDP/sddp.jl");
    # include("src/UnitCommitment/SDDP/extForm.jl");

    Output_Gap = false; max_iter = 150; MaxIter = 100; cutSelection = "SMC"; δ = .1; numScenarios = 100; tightness = true; case = "case30"; # "RTS_GMLC", "case30"
    T = 3; num = 3; TimeLimit = 60 * 60 * 2.; OPT = Inf; 
end


for cut in ["LC", "ELC", "SMC"]
    for num in [3, 5, 10]
        for T in [6, 8, 12]
            indexSets = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/indexSets.jld2")["indexSets"];
            paramOPF = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/paramOPF.jld2")["paramOPF"];
            paramDemand = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/paramDemand.jld2")["paramDemand"];
            scenarioTree = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/scenarioTree.jld2")["scenarioTree"];
            initialStageDecision = load("src/UnitCommitment/experiment_$case/initialStageDecision.jld2")["initialStageDecision"];
            cutSelection = cut
            @everywhere begin
                indexSets = $indexSets; paramOPF = $paramOPF; paramDemand = $paramDemand; scenarioTree = $scenarioTree; initialStageDecision = $initialStageDecision;
            end
            
            # Ξ = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/Ξ.jld2")["Ξ"]
            # @time extResult = extensive_form(indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF, scenarioTree = scenarioTree, Ξ = Ξ, initialStageDecision = initialStageDecision); OPT = extResult.OPT;
            sddpResult = SDDP_algorithm(scenarioTree = scenarioTree, 
                                indexSets = indexSets, 
                                    paramDemand = paramDemand, 
                                        paramOPF = paramOPF, 
                                            initialStageDecision = initialStageDecision, numScenarios = numScenarios, TimeLimit = TimeLimit, OPT = OPT,
                                            Output_Gap = Output_Gap, max_iter = max_iter, MaxIter = MaxIter, δ = δ, cutSelection = cutSelection, tightness = tightness)
            save("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddpResult-$cutSelection-$tightness.jld2", "sddpResult", sddpResult)
        end
    end
end


# Load the results
T = 12; num = 10; cutSelection = "SMC";
sddpResult = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddpResult-$cutSelection-$tightness.jld2")["sddpResult"][:solHistory]