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


    include("src/UnitCommitment_RTS/pSDDLP/def.jl")
    include("src/UnitCommitment_RTS/pSDDLP/backwardModel.jl");
    include("src/UnitCommitment_RTS/pSDDLP/forwardModel.jl");
    include("src/UnitCommitment_RTS/pSDDLP/LevelSetMethod.jl");
    include("src/UnitCommitment_RTS/pSDDLP/sddip.jl");
    # include("src/UnitCommitment/pSDDiP/extForm.jl");


    Output_Gap = false; max_iter = 150; MaxIter = 100; δ = 1.; numScenarios = 100; tightness = true;  TimeLimit = 60 * 60 * 2.; OPT = Inf; case = "case_RTS_GMLC"; # "case_RTS_GMLC", "case30"
    T = 6; num = 3; cutSelection = "ELC"; # "LC", "ELC", "SMC"
    ℓ = .0; core_point_strategy = "Eps"; # "Mid", "In-Out", "Eps", "Relint", "Conv"
end


for cut in ["LC", "ELC", "SMC"]
    for num in [3, 5, 10]
        for T in [6, 8, 12] 
            indexSets = load("src/UnitCommitment_RTS/experiment_$case/stage($T)real($num)/indexSets.jld2")["indexSets"];
            paramOPF = load("src/UnitCommitment_RTS/experiment_$case/stage($T)real($num)/paramOPF.jld2")["paramOPF"];
            paramDemand = load("src/UnitCommitment_RTS/experiment_$case/stage($T)real($num)/paramDemand.jld2")["paramDemand"];
            scenarioTree = load("src/UnitCommitment_RTS/experiment_$case/stage($T)real($num)/scenarioTree.jld2")["scenarioTree"];
            initialStageDecision = load("src/UnitCommitment_RTS/experiment_$case/initialStageDecision.jld2")["initialStageDecision"];
            # cutSelection = cut;
            @everywhere begin
                indexSets = $indexSets; paramOPF = $paramOPF; paramDemand = $paramDemand; scenarioTree = $scenarioTree; initialStageDecision = $initialStageDecision;
            end
            # Ξ = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/Ξ.jld2")["Ξ"]
            # @time extResult = extensive_form(indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF, scenarioTree = scenarioTree, Ξ = Ξ, initialStageDecision = initialStageDecision); OPT = extResult.OPT;
            sddlpResult = SDDiP_algorithm(scenarioTree = scenarioTree, 
                                indexSets = indexSets, 
                                    paramDemand = paramDemand, 
                                        paramOPF = paramOPF, 
                                            initialStageDecision = initialStageDecision, numScenarios = numScenarios, TimeLimit = TimeLimit, OPT = OPT, tightness = tightness,
                                            Output_Gap = Output_Gap, max_iter = max_iter, MaxIter = MaxIter, δ = δ, cutSelection = cutSelection, core_point_strategy = core_point_strategy)
            if cutSelection == "ELC" 
                save("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/MagnantiWong/sddlpResult-$ℓ-$tightness.jld2", "sddlpResult", sddlpResult)
            else
                save("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddlpResult-$cutSelection-$tightness.jld2", "sddlpResult", sddlpResult)
            end
            
            # remove all redundant variables for processes
            @everywhere begin
                forwardInfoList = nothing;
                backwardInfoList = nothing;
                StateVarList = nothing; 
                stageDecision = nothing; 
                λ₀ = nothing; λ₁ = nothing;
                solCollection = nothing;  # to store every iteration results
            end
            @everywhere GC.gc(); # garbage collection
        end
    end
end

T = 12; num = 10; cutSelection = "ELC"; ℓ = 0.0;
if cutSelection == "ELC" 
    sddlpResult = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/MagnantiWong/sddlpResult-$ℓ-$tightness.jld2")["sddlpResult"][:solHistory]
else
    sddlpResult = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddlpResult-$cutSelection-$tightness.jld2")["sddlpResult"][:solHistory]
end
