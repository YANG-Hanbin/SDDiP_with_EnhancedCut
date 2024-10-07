#############################################################################################
####################################    Without Parallel   ##################################
#############################################################################################
using Pkg;
Pkg.activate(".");
using JuMP, Gurobi, ParallelDataTransfer;
using Distributions, Statistics, StatsBase, Distributed, Random;
using Test, Dates, Printf;
using CSV, DataFrames;
using JLD2, FileIO;


const GRB_ENV = Gurobi.Env();

project_root = @__DIR__;

include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "def.jl"));
include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "backwardPass.jl"));
include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "forwardPass.jl"));
include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "extFormGurobi.jl"));
include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "LevelSetMethod.jl"));
include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "setting.jl"));
include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "SDDiP.jl"));


#############################################################################################
####################################    main function   #####################################
#############################################################################################
max_iter = 100; ϵ = 1e-4; cutSelection = "ShrinkageLC"; M = 30; Output_Gap = false; tightness = false; TimeLimit = 60*60*2; MaxIter = 100;
T = 10; num = 3;
cutSelection = "ShrinkageLC"; # "LC", "ShrinkageLC", "ELC"
for cutSelection in ["LC", "ELC", "ShrinkageLC"]
    for T in [5, 10]
        for num in [3, 5]
            stageDataList = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data", "testData_stage($T)_real($num)", "stageDataList.jld2"))["stageDataList"];
            binaryInfo = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data", "testData_stage($T)_real($num)", "binaryInfo.jld2"))["binaryInfo"];
            scenario_sequence = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data", "testData_stage($T)_real($num)", "scenario_sequence.jld2"))["scenario_sequence"];
            probList = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data", "testData_stage($T)_real($num)", "probList.jld2"))["probList"];
            Ω = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data", "testData_stage($T)_real($num)", "Ω.jld2"))["Ω"];

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

