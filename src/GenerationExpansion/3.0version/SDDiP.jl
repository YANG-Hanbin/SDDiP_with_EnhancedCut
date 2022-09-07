function setupLevelsetPara(t::Int64, j::Int64, k::Int64;
                                  cutSelection::String = "ELC", 
                                    L̂::Vector{Float64} = L̂, 
                                    binaryInfo::BinaryInfo = binaryInfo, 
                                    Output_Gap::Bool = false,
                                    λ::Union{Float64, Nothing} = .3, ℓ1::Real = 1.0, ℓ2::Real = 0.8)
    if cutSelection == "ELC"
        forward_modify_constraints!(forwardInfoList[t], 
                                        stageDataList[t], 
                                        Ω[t][j].d, 
                                        solCollection[t-1,k].stageSolution, 
                                        binaryInfo = binaryInfo
                                        )
        optimize!(forwardInfoList[t].model); f_star_value = JuMP.objective_value(forwardInfoList[t].model);
        L̃ = solCollection[t-1,k].stageSolution .* ℓ2 .+ (1 - ℓ2)/2;

        Output = 0; threshold = 1.0; 
        levelSetMethodParam = LevelSetMethodParam(0.95, λ, threshold, 
                                                            1e14, 1e2, Output, Output_Gap,
                                                                        L̂, cutSelection, L̃, f_star_value)


    elseif cutSelection == "ShrinkageLC" 
        forward_modify_constraints!(forwardInfoList[t], 
                                        stageDataList[t], 
                                        Ω[t][j].d, 
                                        solCollection[t-1,k].stageSolution, 
                                        binaryInfo = binaryInfo
                                        )
        optimize!(forwardInfoList[t].model); f_star_value = JuMP.objective_value(forwardInfoList[t].model);
        Output = 0; threshold = 1.0; 
        levelSetMethodParam = LevelSetMethodParam(0.95, λ, threshold, 
                                                            1e14, 1e2, Output, Output_Gap,
                                                                        L̂,  cutSelection, nothing, f_star_value)

    elseif cutSelection == "ELCwithoutConstraint" 
        L̃ = solCollection[t-1,k].stageSolution .* ℓ2 .+ (1 - ℓ2)/2;
        Output = 0; threshold = 1.0; f_star_value = 0.0;
        levelSetMethodParam = LevelSetMethodParam(0.95, λ, threshold, 
                                                            1e14, 1e2, Output, Output_Gap,
                                                                        L̂,  cutSelection, L̃, f_star_value)

    elseif cutSelection == "LC" 
        L̃ = nothing; f_star_value = 0.0;
        Output = 0; threshold = 1.0; 
        levelSetMethodParam = LevelSetMethodParam(0.95, λ, threshold, 
                                                            1e14, 1e2, Output, Output_Gap,
                                                                        L̂,  cutSelection, L̃, f_star_value)
    end

    x₀ =  L̂ .* f_star_value .* ℓ1 .- f_star_value * (ℓ1 ./ 2);

    return (levelSetMethodParam = levelSetMethodParam, x₀ = x₀)
end



function SDDiP_algorithm(   Ω::Dict{Int64,Dict{Int64,RandomVariables}}, 
                            probList::Dict{Int64,Vector{Float64}}, 
                            stageDataList::Dict{Int64, StageData}; 
                            scenario_sequence::Dict{Int64, Dict{Int64, Any}} = scenario_sequence, ϵ::Float64 = 0.001, M::Int64 = 30, max_iter::Int64 = 200, 
                            cutSelection::String = "ELC", binaryInfo::BinaryInfo = binaryInfo)
    ## d: dimension of x
    ## M: num of scenarios when doing one iteration
    initial = now(); iter_time = 0; total_Time = 0; t0 = 0.0;
    T = length(keys(Ω));
    i = 1; LB = - Inf; UB = Inf; solCollection = Dict(); u = 0;Scenarios = 0;

    col_names = [:iter, :LB, :OPT, :UB, :gap, :time, :Time]; # needs to be a vector Symbols
    col_types = [Int64, Float64, Float64, Float64, String, Float64, Float64];
    named_tuple = (; zip(col_names, type[] for type in col_types )...);
    sddipResult = DataFrame(named_tuple); # 0×7 DataFrame
    gapList = [];
    # @time gurobiResult = gurobiOptimize!(Ω, 
    #                                 probList, 
    #                                 stageDataList,
    #                                 binaryInfo = binaryInfo);
    # OPT = gurobiResult.OPT # 1.58e8, 3.8e8, 2.05e9    time 120s 720s 791s
    OPT = 1.58e8;
    forwardInfoList = Dict{Int, ForwardModelInfo}();
    backwardInfoList = Dict{Int, BackwardModelInfo}();
    for t in 1:T 
        forwardInfoList[t] = forwardModel!(stageDataList[t], binaryInfo = binaryInfo)
        backwardInfoList[t] = backwardModel!(stageDataList[t], binaryInfo = binaryInfo)
    end 
    
    println("---------------- print out iteration information -------------------")
    while true
        t0 = now();
        M = 4;
        solCollection = Dict();  # to store every iteration results
        u = Vector{Float64}(undef, M);  # to compute upper bound 
        Random.seed!(i*3)
        Scenarios = SampleScenarios(scenario_sequence, T = T, M = M);
        
        ## Forward Step
        for k in 1:M
             L̂= [0.0 for i in 1:binaryInfo.n];  ## for the first-stage subproblem, we create a zero vector as 'x_ancestor'
            for t in 1:T
                forwardInfo = forwardInfoList[t];
                ## realization of k-th scenario at stage t
                ω = scenario_sequence[Scenarios[k]][1][t];
                ## the following function is used to (1). change the problem coefficients for different node within the same stage t.
                forward_modify_constraints!(forwardInfo, 
                                                stageDataList[t], 
                                                Ω[t][ω].d, 
                                                L̂, 
                                                binaryInfo = binaryInfo
                                                );
                optimize!(forwardInfo.model);

                solCollection[t, k] = ( stageSolution = round.(JuMP.value.(forwardInfo.Lt)), 
                                        stageValue = JuMP.objective_value(forwardInfo.model) - JuMP.value(forwardInfo.θ),
                                        OPT = JuMP.objective_value(forwardInfo.model)
                                        );

                 L̂ = solCollection[t, k].stageSolution;
            end
            u[k] = sum(solCollection[t, k].stageValue for t in 1:T);
        end

        ## compute the upper bound
        LB = solCollection[1, 1].OPT;
        μ̄ = mean(u);
        σ̂² = var(u);
        UB = μ̄ + 1.96 * sqrt(σ̂²/M); # minimum([μ̄ + 1.96 * sqrt(σ̂²/M), UB]);
        gap = round((UB-LB)/UB * 100 ,digits = 2);
        gapString = string(gap,"%");
        push!(sddipResult, [i, LB, OPT, UB, gapString, iter_time, total_Time]); push!(gapList, gap);
        if i == 1
            println("---------------------------------- Iteration Info ------------------------------------")
            println("Iter |   LB                              UB                             gap")
        end
        @printf("%3d  |   %5.3g                         %5.3g                              %1.3f%s\n", i, LB, UB, gap, "%")
        if OPT-LB ≤ 1e-2 * OPT || i > max_iter
            return Dict(:solHistory => sddipResult, 
                            :solution => solCollection[1, 1].stageSolution, 
                            :gapHistory => gapList) 
        end

        for t = reverse(2:T)
            for k in 1:M 
                c = [0, zeros(Float64,binaryInfo.n)]
                for j in keys(Ω[t])
                    # @info "$t, $k, $j"
                    backwardInfo = backwardInfoList[t];
                    backward_Constraint_Modification!(backwardInfo, 
                                                            Ω[t][j].d );
                    # copyVariable_Constraint!(backwardInfo, 
                    #                                 stageDataList[t], 
                    #                                 Ω[t-1][j].d;
                    #                                 binaryInfo = binaryInfo             
                    #                                 );
                    (levelSetMethodParam, x₀) = setupLevelsetPara(t, j , k; cutSelection = "LC", # "ShrinkageLC", "ELCwithoutConstraint", "LC", "ELC"
                                                                L̂ = solCollection[t-1,k].stageSolution, 
                                                                        binaryInfo = binaryInfo, 
                                                                            Output_Gap = false,
                                                                                λ = 0.99, ℓ1 = .0, 
                                                                                            ℓ2 = .9);
                    c = c + probList[t][j] .* LevelSetMethod_optimization!(backwardInfo, x₀; 
                                                                                levelSetMethodParam = levelSetMethodParam, 
                                                                                        stageData = stageDataList[t],   
                                                                                                    ϵ = 1e-3,      
                                                                                                        binaryInfo = binaryInfo
                                                                        ) 
                end
                # # add cut to both backward models and forward models
                # cutCollection[t-1].v[i][k] = c[1]
                # cutCollection[t-1].π[i][k] = c[2]

                @constraint(forwardInfoList[t-1].model, forwardInfoList[t-1].θ ≥ c[1] + c[2]' * forwardInfoList[t-1].Lt )
                @constraint(forwardInfoList[t-1].model, forwardInfoList[t-1].θ ≥ c[1] + c[2]' * forwardInfoList[t-1].Lt ) 
                @constraint(backwardInfoList[t-1].model, backwardInfoList[t-1].θ ≥ c[1] + c[2]' * backwardInfoList[t-1].Lt )
                @constraint(backwardInfoList[t-1].model, backwardInfoList[t-1].θ ≥ c[1] + c[2]' * backwardInfoList[t-1].Lt )              
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




## ------------------------- review --------------------- ##
# 1. if we set the interior point is close the center, when need lambda to be close to 1
# 2. ...