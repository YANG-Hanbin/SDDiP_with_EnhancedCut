using Pkg
Pkg.activate(".")

using Distributed
addprocs(5)

@everywhere begin
    using JuMP, Gurobi, ParallelDataTransfer
    using Distributions, Statistics, StatsBase, Distributed, Random
    using Test, Dates, Printf
    using CSV, DataFrames
    using JLD2, FileIO

    # One Gurobi environment per worker
    const GRB_ENV = Gurobi.Env()

    const PROJECT_ROOT = @__DIR__

    include(joinpath(PROJECT_ROOT, "src", "GenerationExpansion", "utilities", "structs.jl"))
    include(joinpath(PROJECT_ROOT, "src", "GenerationExpansion", "forwardPass.jl"))
    include(joinpath(PROJECT_ROOT, "src", "GenerationExpansion", "backwardPass.jl"))
    include(joinpath(PROJECT_ROOT, "src", "GenerationExpansion", "level_method.jl"))
    include(joinpath(PROJECT_ROOT, "src", "GenerationExpansion", "utilities", "setting.jl"))
    include(joinpath(PROJECT_ROOT, "src", "GenerationExpansion", "utilities", "utils.jl"))
    include(joinpath(PROJECT_ROOT, "src", "GenerationExpansion", "cut_variants.jl"))
    include(joinpath(PROJECT_ROOT, "src", "GenerationExpansion", "sddp.jl"))
end

timeSDDP        = 3600.0;
gapSDDP         = 1e-3;
iterSDDP        = 300;
sample_size_SDDP= 100;
solverGap = 1e-4; solverTime = 20.0;
ε               = 1e-4;
discreteZ       = true;
cutType         = :SMC;
cutSparsity     = true;
verbose         = false;
ℓ1              = 0.0;
ℓ2              = 0.0;
nxt_bound       = 1e8;
logger_save     = false;
algorithm       = :SDDPL;
T               = 10;
num             = 10;


stageDataList = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/stageDataList.jld2")["stageDataList"]
Ω = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/Ω.jld2")["Ω"]
binaryInfo = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/binaryInfo.jld2")["binaryInfo"]
probList = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/probList.jld2")["probList"]

param = param_setup(
    timeSDDP        = timeSDDP,
    gapSDDP         = gapSDDP,
    iterSDDP        = iterSDDP,
    sample_size_SDDP= sample_size_SDDP,
    solverGap       = solverGap,
    solverTime      = solverTime,
    ε               = ε,
    discreteZ       = discreteZ,
    cutType         = cutType,
    cutSparsity     = cutSparsity,
    T               = T,
    num             = num,
    verbose         = verbose,
    ℓ1              = ℓ1,
    ℓ2              = ℓ2,
    nxt_bound       = nxt_bound,
    logger_save     = logger_save,
    algorithm       = algorithm,
);
@everywhere begin
    stageDataList = $stageDataList; Ω = $Ω; binaryInfo = $binaryInfo; probList = $probList; 
    param = $param;
end

sddipResults = stochastic_dual_dynamic_programming_algorithm(
    Ω, 
    probList, 
    stageDataList; 
    binaryInfo = binaryInfo,
    param
);
@everywhere GC.gc(); # garbage collection