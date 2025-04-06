using Pkg
Pkg.activate(".")
using Distributed; addprocs(5); 
@everywhere begin
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
end

#############################################################################################
####################################    main function   #####################################
#############################################################################################
MaxIter = 200; ε = 1e-4; M = 500; Output_Gap = false; tightness = false; TimeLimit = 60*60; 
T = 15; num = 10;
cutSelection = "LC"; # "LC", "ShrinkageLC", "ELC"
FeasibilityTol  = 1e-6;
MIPFocus        = 0;
NumericFocus    = 3;
MIPGap          = 1e-4;
improvement     = false;
logger_save     = true;
ℓ1              = 0.0; 
ℓ2              = 0.0;
nxt_bound = 1e6;

for cutSelection in ["ELC", "ShrinkageLC" ,"LC"]
    for T in [10, 15]
        for num in [5, 10]
            stageDataList = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data_1", "testData_stage($T)_real($num)", "stageDataList.jld2"))["stageDataList"];
            binaryInfo = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data_1", "testData_stage($T)_real($num)", "binaryInfo.jld2"))["binaryInfo"];
            # scenario_sequence = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data_1", "testData_stage($T)_real($num)", "scenario_sequence.jld2"))["scenario_sequence"];
            probList = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data_1", "testData_stage($T)_real($num)", "probList.jld2"))["probList"];
            Ω = load(joinpath(project_root, "src", "GenerationExpansion", "numerical_data_1", "testData_stage($T)_real($num)", "Ω.jld2"))["Ω"];

            if cutSelection == "LC"
                nxt_bound = 1e10
            elseif cutSelection == "ShrinkageLC"
                nxt_bound = 1e5
            elseif cutSelection == "ELC"
                nxt_bound = 1e5
            end
            param = param_setup(;
                terminate_time         = 60*60.,
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
                Output_Gap             = Output_Gap,
                nxt_bound              = nxt_bound,
                logger_save            = logger_save,
                algorithm              = :SDDP
            );
            @everywhere begin
                stageDataList = $stageDataList; Ω = $Ω; binaryInfo = $binaryInfo; probList = $probList; 
                param = $param;
            end

            sddipResults = SDDiP_algorithm(
                Ω, 
                probList, 
                stageDataList; 
                Output_Gap, 
                binaryInfo = binaryInfo,
                param
            );
            @everywhere GC.gc(); # garbage collection
        end
    end
end
