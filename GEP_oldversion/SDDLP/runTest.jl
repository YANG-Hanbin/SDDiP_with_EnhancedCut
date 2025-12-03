#############################################################################################
###################################### Parameter Setup ######################################
#############################################################################################
MaxIter         = 200; 
ε               = 1e-4; 
M               = 500; 
Output_Gap      = false; 
tightness       = true; 
cut             = "LC"; # "LC", "ShrinkageLC", "ELC", "SBC"
sparse_cut      = true; 
logger_save     = false;
nxt_bound       = 1e8;
ℓ1              = 0.0; 
ℓ2              = 0.0;
T = 10; num = 5; 

for cut in ["ELC", "ShrinkageLC" ,"LC", "SBC"]
    for T in [10, 15]
        for num in [5, 10]
            stageDataList = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/stageDataList.jld2")["stageDataList"]
            Ω = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/Ω.jld2")["Ω"]
            binaryInfo = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/binaryInfo.jld2")["binaryInfo"]
            # scenario_sequence = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/scenario_sequence.jld2")["scenario_sequence"]
            probList = load("src/GenerationExpansion/numerical_data/testData_stage($T)_real($num)/probList.jld2")["probList"]

            if cut == "LC"
                nxt_bound = 1e10
            elseif cut == "ShrinkageLC"
                nxt_bound = 1e6
            elseif cut == "ELC"
                nxt_bound = 1e6
            end
            
            param = param_setup(;
                terminate_time         = 60*60.,
                terminate_threshold    = 1e-3,
                MaxIter                = MaxIter,
                M                      = M, 
                ε                      = ε,
                tightness              = tightness,
                cutSelection           = cut,      ## :PLC, :LC, :SMC, :BC, :MDC
                sparse_cut             = sparse_cut,
                T                      = T,
                num                    = num,
                ℓ1                     = ℓ1,
                ℓ2                     = ℓ2,
                Output_Gap             = Output_Gap,
                nxt_bound              = nxt_bound,
                logger_save            = logger_save,
                algorithm              = :SDDLP
            );
            @everywhere begin
                stageDataList = $stageDataList; Ω = $Ω; binaryInfo = $binaryInfo; probList = $probList; 
                param = $param;
            end

            sddipResults = SDDiP_algorithm(
                Ω, 
                probList, 
                stageDataList; 
                binaryInfo = binaryInfo,
                param
            );
            @everywhere GC.gc(); # garbage collection
        end
    end
end