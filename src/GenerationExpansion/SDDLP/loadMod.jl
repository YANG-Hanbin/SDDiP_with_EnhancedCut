using Pkg
Pkg.activate(".")
using Distributed; addprocs(5); 
@everywhere begin
    using JuMP, Gurobi, ParallelDataTransfer
    using Distributions, Statistics, StatsBase, Distributed, Random
    using Test, Dates, Printf
    using CSV, DataFrames
    using JLD2, FileIO


    const GRB_ENV = Gurobi.Env()

    include(joinpath(project_root, "src", "GenerationExpansion", "SDDLP", "def.jl"))
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDLP", "backwardPass.jl"))
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDLP", "forwardPass.jl"))
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDLP", "LevelSetMethod.jl"))
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDLP", "setting.jl"))
    include(joinpath(project_root, "src", "GenerationExpansion", "SDDLP", "SDDiP.jl"))
    # include(joinpath(project_root, "src", "GenerationExpansion", "SDDLP", "extFormGurobi.jl"))
end