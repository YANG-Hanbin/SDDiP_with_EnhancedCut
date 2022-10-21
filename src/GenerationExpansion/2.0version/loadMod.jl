using JuMP, Gurobi, ParallelDataTransfer
using Distributions, Statistics, StatsBase, Distributed, Random
using Test, Dates, Printf
using CSV, DataFrames
using JLD2, FileIO


const GRB_ENV = Gurobi.Env()


include("src/GenerationExpansion/2.0version/def.jl")
include("src/GenerationExpansion/2.0version/backwardPass.jl")
include("src/GenerationExpansion/2.0version/forwardPass.jl")
include("src/GenerationExpansion/2.0version/extFormGurobi.jl")
include("src/GenerationExpansion/2.0version/LevelSetMethod.jl")
include("src/GenerationExpansion/2.0version/setting.jl")

stageDataList = load("src/GenerationExpansion/2.0version/testData/stageDataList.jld2")["stageDataList"]
Ω = load("src/GenerationExpansion/2.0version/testData/Ω.jld2")["Ω"]
binaryInfo = load("src/GenerationExpansion/2.0version/testData/binaryInfo.jld2")["binaryInfo"]
scenario_sequence = load("src/GenerationExpansion/2.0version/testData/scenario_sequence.jld2")["scenario_sequence"]
probList = load("src/GenerationExpansion/2.0version/testData/probList.jld2")["probList"]




#############################################################################################
####################################    main function   #####################################
#############################################################################################
max_iter = 200; ϵ = 1e-2; Enhanced_Cut = true