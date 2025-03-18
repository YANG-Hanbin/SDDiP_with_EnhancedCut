function setupLevelsetPara(forwardInfo::ForwardModelInfo, stageData::StageData, demand::Vector{Float64}, L̂::Vector{Float64};
                                    cutSelection::String = "ELC", 
                                    binaryInfo::BinaryInfo = binaryInfo, 
                                    Output_Gap::Bool = false, max_iter::Int64 = 100,
                                    λ::Union{Float64, Nothing} = .3, ℓ1::Real = 0.0, ℓ2::Real = 1.)
    if cutSelection == "ELC"
        forward_modify_constraints!(forwardInfo, 
                                        stageData, 
                                        demand, 
                                        L̂, 
                                        binaryInfo = binaryInfo
                                        )
        optimize!(forwardInfo.model); f_star_value = JuMP.objective_value(forwardInfo.model);
        # L̃ = L̂ .* ℓ2 .+ (1 - ℓ2)/2;
        L̃ = L̂ ./2;

        Output = 0; threshold = 1.0; 
        levelSetMethodParam = LevelSetMethodParam(0.95, λ, threshold, 
                                                            1e14, max_iter, Output, Output_Gap,
                                                                        L̂, cutSelection, L̃, f_star_value)


    elseif cutSelection == "ShrinkageLC" 
        forward_modify_constraints!(forwardInfo, 
                                        stageData, 
                                        demand, 
                                        L̂, 
                                        binaryInfo = binaryInfo
                                        )
        optimize!(forwardInfo.model); f_star_value = JuMP.objective_value(forwardInfo.model);
        Output = 0; threshold = 1.0; 
        levelSetMethodParam = LevelSetMethodParam(0.95, λ, threshold, 
                                                            1e14, max_iter, Output, Output_Gap,
                                                                        L̂,  cutSelection, nothing, f_star_value)

    elseif cutSelection == "ELCwithoutConstraint" 
        L̃ = L̂ .* ℓ2 .+ (1 - ℓ2)/2;
        Output = 0; threshold = 1.0; f_star_value = 0.0;
        levelSetMethodParam = LevelSetMethodParam(0.95, λ, threshold, 
                                                            1e14, max_iter, Output, Output_Gap,
                                                                        L̂,  cutSelection, L̃, f_star_value)

    elseif cutSelection == "LC" 
        L̃ = nothing; f_star_value = 0.0;
        Output = 0; threshold = 1.0; 
        levelSetMethodParam = LevelSetMethodParam(0.95, λ, threshold, 
                                                            1e14, max_iter, Output, Output_Gap,
                                                                        L̂,  cutSelection, L̃, f_star_value)
    end

    x₀ =  L̂ .* 0.0;

    return (levelSetMethodParam = levelSetMethodParam, x₀ = x₀)
end


function SDDiP_algorithm(   
    Ω::Dict{Int64,Dict{Int64,RandomVariables}}, 
    probList::Dict{Int64,Vector{Float64}}, 
    stageDataList::Dict{Int64, StageData}; 
    Output_Gap::Bool = false, 
    binaryInfo::BinaryInfo = binaryInfo,
    param::NamedTuple = param
)
    
    cutSelection = param.cutSelection; tightness = param.tightness; 
    OPT = Inf;
    # @time gurobiResult = gurobiOptimize!(Ω, 
    #                                 probList, 
    #                                 stageDataList;
    #                                 binaryInfo = binaryInfo, mipGap = 1e-2);
    # OPT = gurobiResult.OPT;
    i = 1; LB = - Inf; UB = Inf; solCollection = Dict(); u = 0;Scenarios = 0;

    col_names = [:iter, :LB, :OPT, :UB, :gap, :time, :Time]; # needs to be a vector Symbols
    col_types = [Int64, Float64, Float64, Float64, String, Float64, Float64];
    named_tuple = (; zip(col_names, type[] for type in col_types )...);
    sddipResult = DataFrame(named_tuple); # 0×7 DataFrame
    gapList = [];
    forwardInfoList = Dict{Int, ForwardModelInfo}();
    backwardInfoList = Dict{Int, BackwardModelInfo}();
    for t in 1:param.T 
        forwardInfoList[t] = forwardModel!(
            stageDataList[t],
            param, 
            binaryInfo = binaryInfo
        )
        backwardInfoList[t] = backwardModel!(
            stageDataList[t], 
            param, 
            binaryInfo = binaryInfo,                                      
            tightness = tightness
        );
    end 
    initial = now(); iter_time = 0.; total_Time = 0.; t0 = 0.0;

    while true
        t0 = now();
        solCollection = Dict();  # to store every iteration results
        u = Vector{Float64}(undef, param.M);  # to compute upper bound 
        Random.seed!(i);
        Scenarios = SampleScenarios(
            Ω, 
            probList, 
            M = param.M
        );
        
        ## Forward Step
        for k in 1:param.M
            Ŝ = [0.0 for i in 1:binaryInfo.d];  ## for the first-stage subproblem, we create a zero vector as 'x_ancestor'
            for t in 1:param.T
                forwardInfo = forwardInfoList[t];
                ## realization of k-th scenario at stage t
                ω = Scenarios[k][t];
                ## the following function is used to (1). change the problem coefficients for different node within the same stage t.
                forward_modify_constraints!(forwardInfo, 
                                                stageDataList[t], 
                                                Ω[t][ω].d, 
                                                Ŝ, 
                                                binaryInfo = binaryInfo
                                                );
                optimize!(forwardInfo.model);

                solCollection[t, k] = ( stageSolution = round.(JuMP.value.(forwardInfo.St), digits = 3), 
                                        stageValue = JuMP.objective_value(forwardInfo.model) - JuMP.value(forwardInfo.θ),
                                        OPT = JuMP.objective_value(forwardInfo.model)
                                        );

                Ŝ  = solCollection[t, k].stageSolution;
            end
            u[k] = sum(solCollection[t, k].stageValue for t in 1:param.T);
        end

        ## compute the upper bound
        LB = solCollection[1, 1].OPT;
        μ̄ = mean(u); UB = μ̄;
        σ̂² = var(u);
        UB = μ̄ + 1.96 * sqrt(σ̂²/param.M); # minimum([μ̄ + 1.96 * sqrt(σ̂²/M), UB]);
        gap = round((UB-LB)/UB * 100, digits = 2);
        gapString = string(gap,"%");
        push!(sddipResult, [i, LB, OPT, UB, gapString, iter_time, total_Time]); push!(gapList, gap);
        if i == 1
            print_iteration_info_bar();
        end
        print_iteration_info(i, LB, UB, gap, iter_time, 0, total_Time); 
        save_info(
            param, 
            Dict(
                :solHistory => sddipResult, 
                :gapHistory => gapList
            );
            logger_save = param.logger_save
        );
        if total_Time > param.terminate_time # || UB-LB ≤ param.terminate_threshold * UB || i >= param.MaxIter
            return Dict(
                :solHistory => sddipResult, 
                :solution => solCollection[1, 1].stageSolution, 
                :gapHistory => gapList
            ) 
        end

        for t = reverse(2:param.T)
            for k in [1] 
                c = [0, zeros(Float64,binaryInfo.d)]
                for j in keys(Ω[t])
                    # @info "$t, $k, $j"
                    backwardInfo = backwardInfoList[t]
                    backward_Constraint_Modification!(backwardInfo, 
                                                            Ω[t][j].d );
                    (levelSetMethodParam, x₀) = setupLevelsetPara(forwardInfoList[t], stageDataList[t], Ω[t][j].d, solCollection[t-1,k].stageSolution;
                                                                        cutSelection = cutSelection,  # "ShrinkageLC", "ELCwithoutConstraint", "LC", "ELC"
                                                                        binaryInfo = binaryInfo, Output_Gap = Output_Gap, max_iter = param.MaxIter,
                                                                        λ = .3, ℓ1 = 0.0, ℓ2 = 1.);
                    c = c + probList[t][j] .* LevelSetMethod_optimization!(
                        backwardInfo, 
                        x₀; 
                        levelSetMethodParam = levelSetMethodParam, 
                        stageData = stageDataList[t],   
                        ϵ = param.ε, 
                        binaryInfo = binaryInfo
                    ); 

                end
                # # add cut to both backward models and forward models
                @constraint(forwardInfoList[t-1].model, forwardInfoList[t-1].θ ≥ c[1] + c[2]' * forwardInfoList[t-1].St )
                @constraint(backwardInfoList[t-1].model, backwardInfoList[t-1].θ ≥ c[1] + c[2]' * backwardInfoList[t-1].St )              
            end
        end

        t1 = now();
        iter_time = (t1 - t0).value/1000;
        total_Time = (t1 - initial).value/1000;
        
        i = i + 1;

    end

end

# using DataFrames
# using Latexify
# df = DataFrame(A = 'x':'z', B = ["M", "F", "F"])
# latexify(df; env=:table, latex=false)

# save("src/GenerationExpansion/3.0version/testData_5/sddipResult_LCinterval.jld2", "sddipResult", sddipResult)
# save("src/GenerationExpansion/3.0version/testData_5/sddipResult_LCbinary.jld2", "sddipResult", sddipResult)
# save("src/GenerationExpansion/3.0version/testData_5/sddipResult_ELCinterval0.jld2", "sddipResult", sddipResult)
# save("src/GenerationExpansion/3.0version/testData_5/sddipResult_ELCbinary0.jld2", "sddipResult", sddipResult)



## ------------------------- review --------------------- ##
# 1. if we set the interior point is close the center, when need lambda to be close to 1
# 2. ...