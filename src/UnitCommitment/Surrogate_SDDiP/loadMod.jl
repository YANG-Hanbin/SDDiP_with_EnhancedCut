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



include("src/UnitCommitment/Surrogate_SDDiP/def.jl")
include("src/UnitCommitment/Surrogate_SDDiP/backwardModel.jl");
include("src/UnitCommitment/Surrogate_SDDiP/forwardModel.jl");
include("src/UnitCommitment/Surrogate_SDDiP/LevelSetMethod.jl");
include("src/UnitCommitment/Surrogate_SDDiP/sddip.jl");
include("src/UnitCommitment/SDDiP/extForm.jl");
include("src/UnitCommitment/Surrogate_SDDiP/readin.jl");