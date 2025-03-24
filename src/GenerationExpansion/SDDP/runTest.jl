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
# include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "extFormGurobi.jl"));
include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "LevelSetMethod.jl"));
include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "setting.jl"));
include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "SDDiP.jl"));


#############################################################################################
####################################    main function   #####################################
#############################################################################################
MaxIter = 200; ε = 1e-4; M = 500; Output_Gap = false; tightness = false; TimeLimit = 60*60; 
T = 15; num = 10;
cutSelection = "ELC"; # "LC", "ShrinkageLC", "ELC"
FeasibilityTol  = 1e-6;
MIPFocus        = 0;
NumericFocus    = 3;
MIPGap          = 1e-4;
improvement     = false;
logger_save = true
for cutSelection in ["ELC", "ShrinkageLC" ,"LC"]
    for T in [10, 15]
        for num in [5, 10]
            stageDataList = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data", "testData_stage($T)_real($num)", "stageDataList.jld2"))["stageDataList"];
            binaryInfo = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data", "testData_stage($T)_real($num)", "binaryInfo.jld2"))["binaryInfo"];
            # scenario_sequence = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data", "testData_stage($T)_real($num)", "scenario_sequence.jld2"))["scenario_sequence"];
            probList = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data", "testData_stage($T)_real($num)", "probList.jld2"))["probList"];
            Ω = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data", "testData_stage($T)_real($num)", "Ω.jld2"))["Ω"];

            param = param_setup(;
                terminate_time         = TimeLimit,
                terminate_threshold    = 1e-3,
                MaxIter                = MaxIter,
                verbose                = false,
                MIPGap                 = MIPGap,
                TimeLimit              = 20,
                MIPFocus               = MIPFocus,
                FeasibilityTol         = FeasibilityTol,
                NumericFocus           = NumericFocus,
                M                      = M, 
                ε                      = ε,
                tightness              = tightness,
                cutSelection           = cutSelection,      ## :PLC, :LC, :SMC, :BC, :MDC
                T                      = T,
                num                    = num,
                logger_save            = logger_save,
                algorithm              = :SDDP
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
sddipResults = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/sddp_tight($tightness)_$cutSelection.jld2")["sddipResults"]
sddipResults[:solHistory]

