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


const GRB_ENV = Gurobi.Env();


include("src/GenerationExpansion/SDDiP/def.jl")
include("src/GenerationExpansion/SDDiP/backwardPass.jl")
include("src/GenerationExpansion/SDDiP/forwardPass.jl")
# include("src/GenerationExpansion/SDDiP/extFormGurobi.jl")
include("src/GenerationExpansion/SDDiP/LevelSetMethod.jl")
include("src/GenerationExpansion/SDDiP/setting.jl")
include("src/GenerationExpansion/SDDiP/SDDiP.jl")


# T = 3; num = 5;
#############################################################################################
####################################    main function   #####################################
#############################################################################################
MaxIter = 200; ε = 1e-3; M = 500; Output_Gap = false; tightness = false; TimeLimit = 60*60.;
cutSelection = "ELC"; # "LC", "ShrinkageLC", "ELC"
T = 10; 
num = 5; 
logger_save = true;
ℓ1 = .0; ℓ2 = 0.0;
for cutSelection in ["ELC", "ShrinkageLC" ,"LC"]
    for T in [10, 15]
        for num in [5, 10]
            stageDataList = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/stageDataList.jld2")["stageDataList"]
            Ω = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/Ω.jld2")["Ω"]
            binaryInfo = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/binaryInfo.jld2")["binaryInfo"]
            # scenario_sequence = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/scenario_sequence.jld2")["scenario_sequence"]
            probList = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/probList.jld2")["probList"]

            param = param_setup(;
                terminate_time         = TimeLimit,
                terminate_threshold    = 1e-3,
                MaxIter                = MaxIter,
                M                      = M, 
                ε                      = ε,
                tightness              = tightness,
                cutSelection           = cutSelection,      ## :PLC, :LC, :SMC, :BC, :MDC
                T                      = T,
                num                    = num,
                ℓ1                     = ℓ1,
                ℓ2                     = ℓ2,
                logger_save            = logger_save,
                algorithm              = :SDDiP
            );

            sddipResults = SDDiP_algorithm(
                Ω, 
                probList, 
                stageDataList; 
                Output_Gap, 
                binaryInfo = binaryInfo,
                param
            );
        end
    end
end

T = 3; num = 5; cutSelection = "LC"; tightness = true; 
sddipResults = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/sddip_tight($tightness)_$cutSelection.jld2")["sddipResults"]
sddipResults[:solHistory]
