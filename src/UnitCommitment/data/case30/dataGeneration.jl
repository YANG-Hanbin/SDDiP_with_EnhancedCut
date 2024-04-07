## generate data, simulation
using JuMP, Gurobi, PowerModels;
using Statistics, StatsBase, Random, Dates, Distributions;
using CSV, DataFrames, Printf;
using JLD2, FileIO;

include("src/UnitCommitment/SDDiP/def.jl");
include("src/UnitCommitment/SDDiP/readin.jl");

case = "case30"
network_data = PowerModels.parse_file("src/UnitCommitment/data/$case/$case.m");
# branchInfo = CSV.read("src/UnitCommitment/data/RTS_GMLC/branch.csv", DataFrame);
T = 8; numRealization = 10;
(indexSets, paramOPF, paramDemand) = prepareIndexSets(T = T, network_data = network_data, scale = 1);

scenarioTree = scenario_tree_generation(T = T, numRealization = numRealization, indexSets = indexSets)
# path = Dict(1 => 1); prob = 1.; Ξ = Dict(); build_scenarios(t = 2, path = path, prob = prob, scenarioTree = scenarioTree);

save("src/UnitCommitment/experiment_$case/stage($T)real($numRealization)/indexSets.jld2", "indexSets", indexSets)
save("src/UnitCommitment/experiment_$case/stage($T)real($numRealization)/paramOPF.jld2", "paramOPF", paramOPF)
save("src/UnitCommitment/experiment_$case/stage($T)real($numRealization)/paramDemand.jld2", "paramDemand", paramDemand)
save("src/UnitCommitment/experiment_$case/stage($T)real($numRealization)/scenarioTree.jld2", "scenarioTree", scenarioTree)
# save("src/UnitCommitment/experiment_$case/stage($T)real($numRealization)/Ξ.jld2", "Ξ", Ξ)
save("src/UnitCommitment/experiment_$case/initialStageDecision.jld2", "initialStageDecision", initialStageDecision)