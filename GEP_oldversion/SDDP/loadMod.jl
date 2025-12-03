using Pkg
Pkg.activate(".")
using Distributed; addprocs(5); 
@everywhere begin
    using JuMP, Gurobi, ParallelDataTransfer;
    using Distributions, Statistics, StatsBase, Distributed, Random;
    using Test, Dates, Printf;
    using CSV, DataFrames;
    using JLD2, FileIO;


    const GRB_ENV = Gurobi.Env();

    project_root = @__DIR__;

    include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "def.jl"));
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "backwardPass.jl"));
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "forwardPass.jl"));
    # include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "extFormGurobi.jl"));
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "LevelSetMethod.jl"));
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "setting.jl"));
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "SDDiP.jl"));
end