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

    include(joinpath(project_root, "src", "multistage_stochastic_unit_commitment", "utilities", "structs.jl"))
    include(joinpath(project_root, "src", "multistage_stochastic_unit_commitment", "utilities", "auxiliary.jl"))
    include(joinpath(project_root, "src", "multistage_stochastic_unit_commitment", "utilities", "level_set_method.jl"))
    include(joinpath(project_root, "src", "multistage_stochastic_unit_commitment", "utilities", "cut_variants.jl"))
    include(joinpath(project_root, "src", "multistage_stochastic_unit_commitment", "utils.jl"))
    include(joinpath(project_root, "src", "multistage_stochastic_unit_commitment", "forward_pass.jl"))
    include(joinpath(project_root, "src", "multistage_stochastic_unit_commitment", "backward_pass.jl"))
    include(joinpath(project_root, "src", "multistage_stochastic_unit_commitment", "partition_tree.jl"))
    include(joinpath(project_root, "src", "multistage_stochastic_unit_commitment", "sddp.jl"))
    # include(joinpath(project_root, "src", "multistage_stochastic_unit_commitment", "utilities", "extForm.jl"))
end