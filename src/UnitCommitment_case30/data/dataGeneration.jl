## generate data, simulation
using Pkg;
Pkg.activate(".");
using JuMP, Gurobi, PowerModels;
using Statistics, StatsBase, Random, Dates, Distributions;
using CSV, DataFrames, Printf;
using JLD2, FileIO;

project_root = @__DIR__;
include(joinpath(project_root, "src", "UnitCommitment_case30", "data", "def.jl"))
include(joinpath(project_root, "src", "UnitCommitment_case30", "data", "readin.jl"))

case = "case30" # case_RTS_GMLC"; "case30";
network_data = PowerModels.parse_file("src/UnitCommitment_case30/data/$case/$case.m");
branchInfo = CSV.read("src/UnitCommitment_case30/data/case_RTS_GMLC/branch.csv", DataFrame);
T = 12; numRealization = 10;

for T in [6, 8, 12]
    for numRealization in [3, 5, 10]
        (indexSets, paramOPF, paramDemand) = prepareIndexSets(T = T, network_data = network_data);

        scenarioTree = scenario_tree_generation(T = T, numRealization = numRealization, indexSets = indexSets)
        # path = Dict(1 => 1); prob = 1.; Ξ = Dict(); build_scenarios(t = 2, path = path, prob = prob, scenarioTree = scenarioTree);

        save("src/UnitCommitment_case30/experiment_$case/stage($T)real($numRealization)/indexSets.jld2", "indexSets", indexSets)
        save("src/UnitCommitment_case30/experiment_$case/stage($T)real($numRealization)/paramOPF.jld2", "paramOPF", paramOPF)
        save("src/UnitCommitment_case30/experiment_$case/stage($T)real($numRealization)/paramDemand.jld2", "paramDemand", paramDemand)
        save("src/UnitCommitment_case30/experiment_$case/stage($T)real($numRealization)/scenarioTree.jld2", "scenarioTree", scenarioTree)
        # save("src/UnitCommitment/experiment_$case/stage($T)real($numRealization)/Ξ.jld2", "Ξ", Ξ)
        # save("src/UnitCommitment/experiment_$case/initialStageDecision.jld2", "initialStageDecision", initialStageDecision)
    end
end
