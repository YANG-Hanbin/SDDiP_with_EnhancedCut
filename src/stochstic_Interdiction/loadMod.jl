using JuMP, Gurobi, ParallelDataTransfer
using Distributions, Statistics, StatsBase, Distributed, Random
using Test, Dates, Printf
using CSV, DataFrames
using JLD2, FileIO


const GRB_ENV = Gurobi.Env()


include("src/stochstic_Interdiction/def.jl")
include("src/stochstic_Interdiction/backwardPass.jl")
include("src/stochstic_Interdiction/forwardPass.jl")
include("src/stochstic_Interdiction/extFormGurobi.jl")
include("src/stochstic_Interdiction/LevelSetMethod.jl")
include("src/stochstic_Interdiction/SDDiP.jl")
include("src/stochstic_Interdiction/generationTest.jl")

(stageData, indexSets, prob, Ω) = dataGeneration(nv = 10, ne = 30, S = 10)


#############################################################################################
####################################    main function   #####################################
#############################################################################################
max_iter = 200; ϵ = 1e-4; cutSelection = "ELC"
