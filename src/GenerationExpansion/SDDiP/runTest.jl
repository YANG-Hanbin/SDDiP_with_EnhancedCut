#############################################################################################
####################################    Without Parallel   ##################################
#############################################################################################
using Pkg
Pkg.activate(".")
using JuMP, Gurobi, ParallelDataTransfer
using Distributions, Statistics, StatsBase, Distributed, Random
using Test, Dates, Printf
using CSV, DataFrames
using JLD2, FileIO


const GRB_ENV = Gurobi.Env()


include("src/GenerationExpansion/SDDiP/def.jl")
include("src/GenerationExpansion/SDDiP/backwardPass.jl")
include("src/GenerationExpansion/SDDiP/forwardPass.jl")
include("src/GenerationExpansion/SDDiP/extFormGurobi.jl")
include("src/GenerationExpansion/SDDiP/LevelSetMethod.jl")
include("src/GenerationExpansion/SDDiP/setting.jl")
include("src/GenerationExpansion/SDDiP/SDDiP.jl")


# T = 3; num = 5;
#############################################################################################
####################################    main function   #####################################
#############################################################################################
max_iter = 100; ϵ = 1e-4; M = 30; Output_Gap = false; tightness = false; 
cutSelection = "ELC"; # "LC", "ShrinkageLC", "ELC"
T = 5; # 3, 5
num = 10; # 5, 10
for cutSelection in ["ELC", "ShrinkageLC" ,"LC"]
    for T in [3, 5]
        for num in [5, 10]
            stageDataList = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/stageDataList.jld2")["stageDataList"]
            Ω = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/Ω.jld2")["Ω"]
            binaryInfo = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/binaryInfo.jld2")["binaryInfo"]
            scenario_sequence = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/scenario_sequence.jld2")["scenario_sequence"]
            probList = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/probList.jld2")["probList"]

            sddipResults = SDDiP_algorithm(Ω, probList, stageDataList, 
                                        scenario_sequence = scenario_sequence,
                                            ϵ = ϵ, M = M, max_iter = max_iter, Output_Gap = Output_Gap, tightness = tightness,
                                                cutSelection = cutSelection, binaryInfo = binaryInfo)
            save("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/sddip_tight($tightness)_$cutSelection.jld2", "sddipResults", sddipResults)
        end
    end
end

T = 3; num = 5; cutSelection = "LC"; tightness = true; 
sddipResults = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/sddip_tight($tightness)_$cutSelection.jld2")["sddipResults"]
sddipResults[:solHistory]
