using Pkg
Pkg.activate(".")
using Distributed; addprocs(3); 
@everywhere begin
    using JuMP, Gurobi, PowerModels;
    using Statistics, StatsBase, Random, Dates, Distributions;
    using Distributed, ParallelDataTransfer;
    using CSV, DataFrames, Printf;
    using JLD2, FileIO;

    const GRB_ENV = Gurobi.Env();

    project_root = @__DIR__;

    include(joinpath(project_root, "src", "UnitCommitment", "SDDLP_improved", "def.jl"))
    include(joinpath(project_root, "src", "UnitCommitment", "SDDLP_improved", "backwardModel.jl"))
    include(joinpath(project_root, "src", "UnitCommitment", "SDDLP_improved", "forwardModel.jl"))
    include(joinpath(project_root, "src", "UnitCommitment", "SDDLP_improved", "LevelSetMethod.jl"))
    include(joinpath(project_root, "src", "UnitCommitment", "SDDLP_improved", "sddip.jl"))

    Output_Gap = false; max_iter = 150; MaxIter = 200; δ = .1; numScenarios = 20; tightness = true; TimeLimit = 60 * 60 * 2.; OPT = Inf; 
    forwardMipGap = 1e-3; backwardMipGap = 1e-3; forwardTimeLimit = 10; backwardTimeLimit = 10; terminate_threshold = 1e-3; branch_threshold = 1e-6; 
    med_method = "interval_mid"; 
    cutSelection = "ELC"; # "SMC"; "LC"; "ELC";
    ℓ = .0; core_point_strategy = "Eps"; # "Mid", "In-Out", "Eps", "Relint", "Conv"
    para = (forwardMipGap = forwardMipGap, backwardMipGap = backwardMipGap, forwardTimeLimit = forwardTimeLimit, backwardTimeLimit = backwardTimeLimit, 
            Output_Gap = Output_Gap, max_iter = max_iter, MaxIter = MaxIter, terminate_threshold = terminate_threshold, branch_threshold = branch_threshold, med_method = med_method, cutSelection = cutSelection, δ = δ, numScenarios = numScenarios, tightness = tightness, TimeLimit = TimeLimit, OPT = OPT, ℓ = ℓ, core_point_strategy = core_point_strategy)
    T = 6; num = 5; case = "case30pwl"; # "case_RTS_GMLC", "case30", "case30pwl", "case24_ieee_rts"
end
# include(joinpath(project_root, "src", "UnitCommitment", "extForm.jl"))  


for cut in ["LC", "ELC", "SMC"]
    for num in [3, 5, 10]
        for T in [6, 8, 12] 
            indexSets = load(joinpath(project_root, "src", "UnitCommitment", "experiment_$case", "stage($T)real($num)", "indexSets.jld2"))["indexSets"];
            paramOPF = load(joinpath(project_root, "src", "UnitCommitment", "experiment_$case", "stage($T)real($num)", "paramOPF.jld2"))["paramOPF"];
            paramDemand = load(joinpath(project_root, "src", "UnitCommitment", "experiment_$case", "stage($T)real($num)", "paramDemand.jld2"))["paramDemand"];
            scenarioTree = load(joinpath(project_root, "src", "UnitCommitment", "experiment_$case", "stage($T)real($num)", "scenarioTree.jld2"))["scenarioTree"];
            initialStageDecision = load(joinpath(project_root, "src", "UnitCommitment", "experiment_$case", "initialStageDecision.jld2"))["initialStageDecision"];

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
                                            initialStageDecision = initialStageDecision, para = para)
            if cutSelection == "ELC" 
                save("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/MagnantiWong/isddlpHCResult-$ℓ-$tightness.jld2", "sddlpResult", sddlpResult)
            else
                save("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/isddlpResult-$cutSelection-$tightness.jld2", "sddlpResult", sddlpResult)
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

T = 6; num = 3; cutSelection = "SMC"; ℓ = 0.0;
if cutSelection == "ELC" 
    sddlpResult = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/MagnantiWong/isddlpResult-$ℓ-$tightness.jld2")["sddlpResult"][:solHistory]
else
    sddlpResult = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/isddlpResult-$cutSelection-$tightness.jld2")["sddlpResult"][:solHistory]
end


sddlpResult = load("src/UnitCommitment/numericalResults-case30pwl/Periods6-Real5/isddlpResult-SMC-true.jld2")["sddlpResult"][:solHistory]

