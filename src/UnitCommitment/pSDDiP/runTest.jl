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


    include("src/UnitCommitment/pSDDiP/def.jl")
    include("src/UnitCommitment/pSDDiP/backwardModel.jl");
    include("src/UnitCommitment/pSDDiP/forwardModel.jl");
    include("src/UnitCommitment/pSDDiP/LevelSetMethod.jl");
    include("src/UnitCommitment/pSDDiP/sddip.jl");
    # include("src/UnitCommitment/pSDDiP/extForm.jl");


    Output_Gap = false; max_iter = 150; MaxIter = 100; cutSelection = "SMC"; δ = 1.; numScenarios = 100; tightness = true; case = "case30"; # "RTS_GMLC", "case30"
    T = 8; num = 5; TimeLimit = 60 * 60 * 2.; OPT = Inf; coef = 10; ε = 1/2^(coef); ℓ = .7; 
end


for cutSelection in ["LC", "ELC", "SMC"]
    for num in [3, 5, 10]
        for T in [6, 8, 12]
            indexSets = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/indexSets.jld2")["indexSets"];
            paramOPF = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/paramOPF.jld2")["paramOPF"];
            paramDemand = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/paramDemand.jld2")["paramDemand"];
            scenarioTree = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/scenarioTree.jld2")["scenarioTree"];
            initialStageDecision = load("src/UnitCommitment/experiment_$case/initialStageDecision.jld2")["initialStageDecision"];
            # cutSelection = cut;
            @everywhere begin
                indexSets = $indexSets; paramOPF = $paramOPF; paramDemand = $paramDemand; scenarioTree = $scenarioTree; initialStageDecision = $initialStageDecision;
            end
            # Ξ = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/Ξ.jld2")["Ξ"]
            # @time extResult = extensive_form(indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF, scenarioTree = scenarioTree, Ξ = Ξ, initialStageDecision = initialStageDecision); OPT = extResult.OPT;
            sddipResult = SDDiP_algorithm(scenarioTree = scenarioTree, 
                                indexSets = indexSets, 
                                    paramDemand = paramDemand, 
                                        paramOPF = paramOPF, MaxIter = MaxIter, ℓ = ℓ,
                                            initialStageDecision = initialStageDecision, numScenarios = numScenarios, TimeLimit = TimeLimit, OPT = OPT,
                                            Output_Gap = Output_Gap, max_iter = max_iter, ε = ε, δ = δ, cutSelection = cutSelection)
            save("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-$coef.jld2", "sddipResult", sddipResult)

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
@everywhere begin T = 6; num = 10; coef = 10; ε = 1/2^(coef); end
sddipResult = SDDiP_algorithm(scenarioTree = scenarioTree, 
                                indexSets = indexSets, 
                                    paramDemand = paramDemand, 
                                        paramOPF = paramOPF, MaxIter = MaxIter, ℓ = ℓ,
                                            initialStageDecision = initialStageDecision, numScenarios = numScenarios, TimeLimit = TimeLimit, OPT = OPT,
                                            Output_Gap = Output_Gap, max_iter = max_iter, ε = ε, δ = δ, cutSelection = cutSelection)
save("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-$coef.jld2", "sddipResult", sddipResult)

sddipResult = load("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-$coef.jld2")["sddipResult"][:solHistory]