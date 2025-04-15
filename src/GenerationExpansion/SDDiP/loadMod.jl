using Pkg
Pkg.activate(".")
using Distributed; addprocs(5); 
@everywhere begin
    using JuMP, Gurobi, ParallelDataTransfer
    using Distributions, Statistics, StatsBase, Distributed, Random
    using Test, Dates, Printf
    using CSV, DataFrames
    using JLD2, FileIO


    const GRB_ENV = Gurobi.Env();

    include(joinpath(project_root, "src", "GenerationExpansion", "SDDiP", "def.jl"))
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDiP", "backwardPass.jl"))
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDiP", "forwardPass.jl"))
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDiP", "LevelSetMethod.jl"))
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDiP", "setting.jl"))
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDiP", "SDDiP.jl"))
    # include(joinpath(project_root, "src", "GenerationExpansion", "SDDiP", "extFormGurobi.jl"))
end