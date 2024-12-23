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
    include(joinpath(project_root, "src", "alg", "utilities", "bundle_method.jl"))
    include(joinpath(project_root, "src", "alg", "utilities", "cut_variants.jl"))
    include(joinpath(project_root, "src", "alg", "utils.jl"))
    include(joinpath(project_root, "src", "alg", "forward_pass.jl"))
    include(joinpath(project_root, "src", "alg", "backward_pass.jl"))
    include(joinpath(project_root, "src", "alg", "partition_tree.jl"))
    include(joinpath(project_root, "src", "alg", "config.jl"))
    include(joinpath(project_root, "src", "alg", "stochastic_dual_dynamic_programming.jl"))

end

case = "case30pwl"
for cut in [:PLC, :SMC, :LC]
    for num in [3, 5, 10]
        for T in [6, 8, 12] 
            param = param_setup(
                numScenarios = 50,
                LiftIterThreshold = 10,
                cutSelection = cut, 
                algorithm = :SDDP,
                T = T,
                num = num,
                case = case
            );
            indexSets        = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "indexSets.jld2"))["indexSets"];
            paramOPF         = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "paramOPF.jld2"))["paramOPF"];
            paramDemand      = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "paramDemand.jld2"))["paramDemand"];
            scenarioTree     = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "scenarioTree.jld2"))["scenarioTree"];
            initialStateInfo = load(joinpath(project_root, "src", "alg", "experiment_$case", "initialStateInfo.jld2"))["initialStateInfo"];

            
            @everywhere begin
                indexSets = $indexSets; paramOPF = $paramOPF; paramDemand = $paramDemand; scenarioTree = $scenarioTree; initialStateInfo = $initialStateInfo;
            end


            sddpResults = stochastic_dual_dynamic_programming_algorithm(
                    scenarioTree,                   
                    indexSets,                        
                    paramDemand,
                    paramOPF;
                    initialStateInfo = initialStateInfo,
                    param_PLC = param_PLC, 
                    param_levelsetmethod = param_levelsetmethod, 
                    param = param
            );
            @everywhere GC.gc(); # garbage collection
        end
    end
end