## generate data, simulation
using JuMP, Gurobi, PowerModels;
using Statistics, StatsBase, Random, Dates, Distributions;
using CSV, DataFrames, Printf;
using JLD2, FileIO;

include("src/UnitCommitment/def.jl");
include("src/UnitCommitment/readin.jl");

network_data = PowerModels.parse_file("src/UnitCommitment/data/RTS_GMLC/case_RTS_GMLC.m");
branchInfo = CSV.read("src/UnitCommitment/data/RTS_GMLC/branch.csv", DataFrame);
T = 2; numRealization = 2;
(indexSets, paramOPF, paramDemand) = prepareIndexSets(T = T, network_data = network_data, branchInfo = branchInfo);

scenarioTree = scenario_tree_generation(T = T, numRealization = numRealization, indexSets = indexSets)
path = Dict(1 => 1); prob = 1.; Ξ = Dict(); build_scenarios(t = 2, path = path, prob = prob, scenarioTree = scenarioTree);

save("src/UnitCommitment/experiment/indexSets.jld2", "indexSets", indexSets)
save("src/UnitCommitment/experiment/paramOPF.jld2", "paramOPF", paramOPF)
save("src/UnitCommitment/experiment/paramDemand.jld2", "paramDemand", paramDemand)
save("src/UnitCommitment/experiment/scenarioTree.jld2", "scenarioTree", scenarioTree)
save("src/UnitCommitment/experiment/Ξ.jld2", "Ξ", Ξ)

