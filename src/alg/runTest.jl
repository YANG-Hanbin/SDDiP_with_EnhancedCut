using Pkg
Pkg.activate(".")
using Distributed; # addprocs(5); 
@everywhere begin
    using JuMP, Gurobi, PowerModels;
    using Statistics, StatsBase, Random, Dates, Distributions;
    using Distributed, ParallelDataTransfer;
    using CSV, DataFrames, Printf;
    using JLD2, FileIO;

    const GRB_ENV = Gurobi.Env();

    project_root = @__DIR__;

    include(joinpath(project_root, "src", "alg", "utils", "structs.jl"))
    include(joinpath(project_root, "src", "alg", "utils", "auxiliary.jl"))
    include(joinpath(project_root, "src", "alg", "utils", "bundle_method.jl"))
    include(joinpath(project_root, "src", "alg", "utils", "cut_variants.jl"))
    include(joinpath(project_root, "src", "alg", "auxiliary.jl"))
    include(joinpath(project_root, "src", "alg", "forward_pass.jl"))
    include(joinpath(project_root, "src", "alg", "backward_pass.jl"))
    include(joinpath(project_root, "src", "alg", "partition_tree.jl"))
    include(joinpath(project_root, "src", "alg", "stochastic_dual_dynamic_programming.jl"))

    param = (
        verbose             = false,
        MIPGap              = 1e-4,
        TimeLimit           = 3600,
        terminate_threshold = 1e-3,
        MaxIter             = 100,
        θ̲                   = 0.0,
        OPT                 = 0.0,
        tightness           = true,
        numScenarios        = 3,
        LiftIterThreshold   = 10,
        branch_threshold    = 1e-3,
        ## "interval_mid", "exact_point"
        med_method          = "interval_mid",   
        ## :PLC, :SMC, :LC
        cutSelection        = :PLC,             
        ## :SDDPL, :SDDP, :SDDiP
        algorithm           = :SDDPL,          
    )


    param_levelsetmethod = (
        μ             = 0.9,
        λ             = 0.5,
        threshold     = 1e-4,
        nxt_bound     = 1e10,
        MaxIter       = 200,
        verbose       = false,
    )
    
    param_PLC = (
        core_point_strategy = "Eps", # "Mid", "Eps"
        δ                   = 1e-3,
        ℓ                   = .0,
    )

    T = 6; num = 3; case = "case30pwl"; # "case_RTS_GMLC", "case30", "case30pwl", "case24_ieee_rts"
    indexSets        = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "indexSets.jld2"))["indexSets"];
    paramOPF         = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "paramOPF.jld2"))["paramOPF"];
    paramDemand      = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "paramDemand.jld2"))["paramDemand"];
    scenarioTree     = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "scenarioTree.jld2"))["scenarioTree"];
    initialStateInfo = load(joinpath(project_root, "src", "alg", "experiment_$case", "initialStateInfo.jld2"))["initialStateInfo"];

end
# @everywhere begin
#     T = 12; num = 10; case = "case30pwl"; # "case_RTS_GMLC", "case30", "case30pwl", "case24_ieee_rts"
#     indexSets        = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "indexSets.jld2"))["indexSets"];
#     paramOPF         = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "paramOPF.jld2"))["paramOPF"];
#     paramDemand      = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "paramDemand.jld2"))["paramDemand"];
#     scenarioTree     = load(joinpath(project_root, "src", "alg", "experiment_$case", "stage($T)real($num)", "scenarioTree.jld2"))["scenarioTree"];
#     initialStateInfo = load(joinpath(project_root, "src", "alg", "experiment_$case", "initialStateInfo.jld2"))["initialStateInfo"];
# end
sddpResult = stochastic_dual_dynamic_programming_algorithm(
        scenarioTree,                   
        indexSets,                        
        paramDemand,
        paramOPF;
        initialStateInfo = initialStateInfo,
        param_PLC = param_PLC, 
        param_levelsetmethod = param_levelsetmethod, 
        param = param
)


