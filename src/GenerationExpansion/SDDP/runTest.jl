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


include("src/GenerationExpansion/SDDP/def.jl")
include("src/GenerationExpansion/SDDP/backwardPass.jl")
include("src/GenerationExpansion/SDDP/forwardPass.jl")
include("src/GenerationExpansion/SDDP/extFormGurobi.jl")
include("src/GenerationExpansion/SDDP/LevelSetMethod.jl")
include("src/GenerationExpansion/SDDP/setting.jl")
include("src/GenerationExpansion/SDDP/SDDiP.jl")

#############################################################################################
####################################    main function   #####################################
#############################################################################################
max_iter = 100; ϵ = 1e-4; cutSelection = "ShrinkageLC"; M = 30; Output_Gap = false; tightness = false; TimeLimit = 60*60*2; MaxIter = 100;
T = 3; num = 5;
cutSelection = "ELC"; # "LC", "ShrinkageLC", "ELC"
for cutSelection in ["LC", "ELC", "ShrinkageLC"]
    for T in [3, 5]
        for num in [5, 10]
            stageDataList = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/stageDataList.jld2")["stageDataList"]
            Ω = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/Ω.jld2")["Ω"]
            binaryInfo = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/binaryInfo.jld2")["binaryInfo"]
            scenario_sequence = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/scenario_sequence.jld2")["scenario_sequence"]
            probList = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/probList.jld2")["probList"]
            sddipResults = SDDiP_algorithm(Ω, probList, stageDataList, 
                                        scenario_sequence = scenario_sequence,
                                            ϵ = ϵ, M = M, max_iter = max_iter, Output_Gap = Output_Gap, 
                                                cutSelection = cutSelection, binaryInfo = binaryInfo, TimeLimit = TimeLimit, MaxIter = MaxIter)
            save("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/sddp_tight($tightness)_$cutSelection.jld2", "sddipResults", sddipResults)  
        end
    end
end

T = 3; num = 5; cutSelection = "LC"; tightness = true; 
sddipResults = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/sddp_tight($tightness)_$cutSelection.jld2")["sddipResults"]
sddipResults[:solHistory]

