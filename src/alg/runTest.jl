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

    project_root = @__DIR__;

    include(joinpath(project_root, "src", "alg", "utilities", "structs.jl"))
    include(joinpath(project_root, "src", "alg", "utilities", "auxiliary.jl"))
    include(joinpath(project_root, "src", "alg", "utilities", "level_set_method.jl"))
    include(joinpath(project_root, "src", "alg", "utilities", "cut_variants.jl"))
    include(joinpath(project_root, "src", "alg", "utils.jl"))
    include(joinpath(project_root, "src", "alg", "forward_pass.jl"))
    include(joinpath(project_root, "src", "alg", "backward_pass.jl"))
    include(joinpath(project_root, "src", "alg", "partition_tree.jl"))
    include(joinpath(project_root, "src", "alg", "sddp.jl"))

end

case = "case30"; # "case_RTS_GMLC", "case30", "case30pwl",
algorithm = :SDDPL; 
cut = :SMC; 
num = 5; T = 12;
for algorithm in [:SDDPL, :SDDP, :SDDiP]
    for cut in [:PLC, :SMC, :LC]
        for num in [3, 5, 10]
            for T in [6, 8, 12] 
                param = param_setup(
                    terminate_time = 7200,
                    TimeLimit = 10,
                    ε = 1/32;
                    numScenarios = 10,
                    LiftIterThreshold = 30,
                    cutSelection = cut, 
                    algorithm = algorithm,
                    terminate_threshold = 1e-3,
                    T = T,
                    num = num,
                    case = case,
                    med_method = :IntervalMed # :IntervalMed, :ExactPoint
                );
                param_cut = param_cut_setup(
                    core_point_strategy = "Eps", # "Mid", "Eps"
                    δ = 1e-2,
                    ℓ = .0,
                );
                param_levelsetmethod = param_levelsetmethod_setup(
                    μ = 0.9,
                    λ = 0.5,
                    threshold = 1e-4,
                    nxt_bound = 1e10,
                    MaxIter = 200,
                    verbose = false
                );
                indexSets        = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "indexSets.jld2"))["indexSets"];
                paramOPF         = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "paramOPF.jld2"))["paramOPF"];
                paramDemand      = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "paramDemand.jld2"))["paramDemand"];
                scenarioTree     = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "scenarioTree.jld2"))["scenarioTree"];
                initialStateInfo = load(joinpath(project_root, "src", "alg", "experiment_$case", "initialStateInfo.jld2"))["initialStateInfo"];

                
                @everywhere begin
                    indexSets = $indexSets; paramOPF = $paramOPF; paramDemand = $paramDemand; scenarioTree = $scenarioTree; initialStateInfo = $initialStateInfo;
                    param_cut = $param_cut; param_levelsetmethod = $param_levelsetmethod; param = $param;
                end


                sddpResults = stochastic_dual_dynamic_programming_algorithm(
                        scenarioTree,                   
                        indexSets,                        
                        paramDemand,
                        paramOPF;
                        initialStateInfo = initialStateInfo,
                        param_cut = param_cut, 
                        param_levelsetmethod = param_levelsetmethod, 
                        param = param
                );
                @everywhere GC.gc(); # garbage collection
            end
        end
    end
end
