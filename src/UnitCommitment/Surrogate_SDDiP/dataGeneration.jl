## generate data, simulation
using JuMP, Gurobi, PowerModels;
using Statistics, StatsBase, Random, Dates, Distributions;
using CSV, DataFrames, Printf;
using JLD2, FileIO;

include("src/def.jl");
include("src/readin.jl");

network_data = PowerModels.parse_file("data/RTS_GMLC/case_RTS_GMLC.m");
branchInfo = CSV.read("data/RTS_GMLC/branch.csv", DataFrame);
T = 3; numRealization = 3;
(indexSets, paramOPF, paramDemand) = prepareIndexSets(T = T, network_data = network_data, branchInfo = branchInfo);

scenarioTree = scenario_tree_generation(T = T, numRealization = numRealization, indexSets = indexSets)
path = Dict(1 => 1); prob = 1.; Ξ = Dict(); build_scenarios(t = 2, path = path, prob = prob, scenarioTree = scenarioTree);

save("src/experiment/indexSets.jld2", "indexSets", indexSets)
save("src/experiment/paramOPF.jld2", "paramOPF", paramOPF)
save("src/experiment/paramDemand.jld2", "paramDemand", paramDemand)
save("src/experiment/scenarioTree.jld2", "scenarioTree", scenarioTree)
save("src/experiment/Ξ.jld2", "Ξ", Ξ)

