using JuMP, Gurobi, ParallelDataTransfer
using Distributions, Statistics, StatsBase, Distributed, Random
using Test, Dates, Printf
using CSV, DataFrames
using JLD2, FileIO


const GRB_ENV = Gurobi.Env()


include("src/GenerationExpansion/3.0version/def.jl")
include("src/GenerationExpansion/3.0version/backwardPass.jl")
include("src/GenerationExpansion/3.0version/forwardPass.jl")
include("src/GenerationExpansion/3.0version/extFormGurobi.jl")
include("src/GenerationExpansion/3.0version/LevelSetMethod.jl")
include("src/GenerationExpansion/3.0version/setting.jl")
include("src/GenerationExpansion/3.0version/SDDiP.jl")

stageDataList = load("src/GenerationExpansion/3.0version/testData2/stageDataList.jld2")["stageDataList"]
Ω = load("src/GenerationExpansion/3.0version/testData2/Ω.jld2")["Ω"]
binaryInfo = load("src/GenerationExpansion/3.0version/testData2/binaryInfo.jld2")["binaryInfo"]
scenario_sequence = load("src/GenerationExpansion/3.0version/testData2/scenario_sequence.jld2")["scenario_sequence"]
probList = load("src/GenerationExpansion/3.0version/testData2/probList.jld2")["probList"]




#############################################################################################
####################################    main function   #####################################
#############################################################################################
max_iter = 200; ϵ = 1e-2; Enhanced_Cut = true