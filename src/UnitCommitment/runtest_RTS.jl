#############################################################################################
####################################    Without Parallel   ##################################
#############################################################################################
using Pkg
Pkg.activate(".")
using JuMP, Gurobi, PowerModels
using Statistics, StatsBase, Random, Dates, Distributions
using Distributed, ParallelDataTransfer
using CSV, DataFrames, Printf
using JLD2, FileIO

const GRB_ENV = Gurobi.Env()


include("src/UnitCommitment/def.jl")
include("src/UnitCommitment/backwardModel.jl");
include("src/UnitCommitment/forwardModel.jl");
include("src/UnitCommitment/LevelSetMethod.jl");
include("src/UnitCommitment/sddip.jl");

indexSets = load("src/UnitCommitment/experiment/indexSets.jld2")["indexSets"]
paramOPF = load("src/UnitCommitment/experiment/paramOPF.jld2")["paramOPF"]
paramDemand = load("src/UnitCommitment/experiment/paramDemand.jld2")["paramDemand"]
scenarioTree = load("src/UnitCommitment/experiment/scenarioTree.jld2")["scenarioTree"]

#############################################################################################
########################################## Run Test #########################################
#############################################################################################
Output_Gap = false; max_iter::Int64 = 100; ϵ = 1e-4; cutSelection = "LC";

sddpResult = SDDiP_algorithm(scenarioTree = scenarioTree, 
                                indexSets = indexSets, 
                                    paramDemand = paramDemand, 
                                        paramOPF = paramOPF, 
                                            Output_Gap = Output_Gap, max_iter = max_iter, ϵ = ϵ, cutSelection = cutSelection)