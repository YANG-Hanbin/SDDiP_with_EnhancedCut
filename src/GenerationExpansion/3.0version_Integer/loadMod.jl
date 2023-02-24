using JuMP, Gurobi, ParallelDataTransfer
using Distributions, Statistics, StatsBase, Distributed, Random
using Test, Dates, Printf
using CSV, DataFrames
using JLD2, FileIO


const GRB_ENV = Gurobi.Env()


include("src/GenerationExpansion/3.0version_Integer/def.jl")
include("src/GenerationExpansion/3.0version_Integer/backwardPass.jl")
include("src/GenerationExpansion/3.0version_Integer/forwardPass.jl")
include("src/GenerationExpansion/3.0version_Integer/extFormGurobi.jl")
include("src/GenerationExpansion/3.0version_Integer/LevelSetMethod.jl")
include("src/GenerationExpansion/3.0version_Integer/setting.jl")
include("src/GenerationExpansion/3.0version_Integer/SDDiP.jl")

stageDataList = load("src/GenerationExpansion/3.0version_Integer/testData_2/stageDataList.jld2")["stageDataList"]
Ω = load("src/GenerationExpansion/3.0version_Integer/testData_2/Ω.jld2")["Ω"]
binaryInfo = load("src/GenerationExpansion/3.0version_Integer/testData_2/binaryInfo.jld2")["binaryInfo"]
scenario_sequence = load("src/GenerationExpansion/3.0version_Integer/testData_2/scenario_sequence.jld2")["scenario_sequence"]
probList = load("src/GenerationExpansion/3.0version_Integer/testData_2/probList.jld2")["probList"]




#############################################################################################
####################################    main function   #####################################
#############################################################################################
max_iter = 200; ϵ = 1e-4; 