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
    
    project_root = @__DIR__;
    include(joinpath(project_root, "src", "UnitCommitment_case30", "SDDP", "def.jl"))
    include(joinpath(project_root, "src", "UnitCommitment_case30", "SDDP", "backwardModel.jl"))
    include(joinpath(project_root, "src", "UnitCommitment_case30", "SDDP", "forwardModel.jl"))
    include(joinpath(project_root, "src", "UnitCommitment_case30", "SDDP", "LevelSetMethod.jl"))
    include(joinpath(project_root, "src", "UnitCommitment_case30", "SDDP", "sddp.jl"))

    Output_Gap = false; max_iter = 150; MaxIter = 200; cutSelection = "ELC"; δ = .1; numScenarios = 100; tightness = true; case = "case30"; # "RTS_GMLC", "case30"
    T = 6; num = 5; TimeLimit = 60 * 60 * 2.; OPT = Inf; ℓ = 0.0;
end
# extForm_path = joinpath(project_root, "src", "UnitCommitment_case30", "extForm.jl")
# include(extForm_path)  

for cut in ["LC", "ELC", "SMC"]
    for num in [3, 5, 10]
        for T in [6, 8, 12]
            indexSets = load(joinpath(project_root, "src", "UnitCommitment_case30", "experiment_$case", "stage($T)real($num)", "indexSets.jld2"))["indexSets"];
            paramOPF = load(joinpath(project_root, "src", "UnitCommitment_case30", "experiment_$case", "stage($T)real($num)", "paramOPF.jld2"))["paramOPF"];
            paramDemand = load(joinpath(project_root, "src", "UnitCommitment_case30", "experiment_$case", "stage($T)real($num)", "paramDemand.jld2"))["paramDemand"];
            scenarioTree = load(joinpath(project_root, "src", "UnitCommitment_case30", "experiment_$case", "stage($T)real($num)", "scenarioTree.jld2"))["scenarioTree"];
            initialStageDecision = load(joinpath(project_root, "src", "UnitCommitment_case30", "experiment_$case", "initialStageDecision.jld2"))["initialStageDecision"];

            # cutSelection = cut
            @everywhere begin
                indexSets = $indexSets; paramOPF = $paramOPF; paramDemand = $paramDemand; scenarioTree = $scenarioTree; initialStageDecision = $initialStageDecision;
            end
            
            # Ξ = load("src/UnitCommitment_case30/experiment_$case/stage($T)real($num)/Ξ.jld2")["Ξ"]
            # @time extResult = extensive_form(indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF, scenarioTree = scenarioTree, Ξ = Ξ, initialStageDecision = initialStageDecision); OPT = extResult.OPT;
            sddpResult = SDDP_algorithm(scenarioTree = scenarioTree, 
                                indexSets = indexSets, 
                                    paramDemand = paramDemand, 
                                        paramOPF = paramOPF, 
                                            initialStageDecision = initialStageDecision, numScenarios = numScenarios, TimeLimit = TimeLimit, OPT = OPT,
                                            Output_Gap = Output_Gap, max_iter = max_iter, MaxIter = MaxIter, δ = δ, cutSelection = cutSelection, tightness = tightness)
            save("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddpResult-$cutSelection-$tightness.jld2", "sddpResult", sddpResult)
        end
    end
end


# Load the results
T = 12; num = 5; cutSelection = "SMC";
sddpResult = load("src/UnitCommitment_case30/numericalResults-$case/Periods$T-Real$num/sddpResult-$cutSelection-$tightness.jld2")["sddpResult"][:solHistory]
