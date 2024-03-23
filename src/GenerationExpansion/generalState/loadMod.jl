using JuMP, Gurobi, ParallelDataTransfer
using Distributions, Statistics, StatsBase, Distributed, Random
using Test, Dates, Printf
using CSV, DataFrames
using JLD2, FileIO


const GRB_ENV = Gurobi.Env()


include("src/GenerationExpansion/generalState/def.jl")
include("src/GenerationExpansion/generalState/backwardPass.jl")
include("src/GenerationExpansion/generalState/forwardPass.jl")
include("src/GenerationExpansion/generalState/extFormGurobi.jl")
include("src/GenerationExpansion/generalState/LevelSetMethod.jl")
include("src/GenerationExpansion/generalState/setting.jl")
include("src/GenerationExpansion/generalState/SDDiP.jl")

T = 3; num = 5;
stageDataList = load("src/GenerationExpansion/data/testData_stage($T)_real($num)/stageDataList.jld2")["stageDataList"]
Ω = load("src/GenerationExpansion/data/testData_stage($T)_real($num)/Ω.jld2")["Ω"]
binaryInfo = load("src/GenerationExpansion/data/testData_stage($T)_real($num)/binaryInfo.jld2")["binaryInfo"]
scenario_sequence = load("src/GenerationExpansion/data/testData_stage($T)_real($num)/scenario_sequence.jld2")["scenario_sequence"]
probList = load("src/GenerationExpansion/data/testData_stage($T)_real($num)/probList.jld2")["probList"]