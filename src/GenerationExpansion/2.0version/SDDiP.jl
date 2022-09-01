function SDDiP_algorithm(   Ω::Dict{Int64,Dict{Int64,RandomVariables}}, 
                            probList::Dict{Int64,Vector{Float64}}, 
                            stageDataList::Dict{Int64, StageData}; 
                            scenario_sequence::Dict{Int64, Dict{Int64, Any}} = scenario_sequence, ϵ::Float64 = 0.001, M::Int64 = 30, max_iter::Int64 = 200, 
                            Enhanced_Cut::Bool = true, binaryInfo::BinaryInfo = binaryInfo)
    ## d: dimension of x
    ## M: num of scenarios when doing one iteration
    initial = now(); iter_time = 0; total_Time = 0; t0 = 0.0;
    T = length(keys(Ω));
    i = 1; LB = - Inf; UB = Inf; 

    col_names = [:iter, :LB, :OPT, :UB, :gap, :time, :Time]; # needs to be a vector Symbols
    col_types = [Int64, Float64, Float64, Float64, String, Float64, Float64];
    named_tuple = (; zip(col_names, type[] for type in col_types )...);
    sddipResult = DataFrame(named_tuple); # 0×7 DataFrame
    gapList = [];
    @time gurobiResult = gurobiOptimize!(Ω, 
                                    probList, 
                                    stageDataList,
                                    binaryInfo = binaryInfo);
    OPT = gurobiResult.OPT
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
             L̂= [0.0 for i in 1:binaryInfo.n];
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
        UB = minimum([μ̄ + 1.96 * sqrt(σ̂²/M), UB]);
        gap = round((UB-LB)/UB * 100 ,digits = 2);
        gapString = string(gap,"%");
        push!(sddipResult, [i, LB, OPT, UB, gapString, iter_time, total_Time]); push!(gapList, gap);
        if i == 1
            println("---------------------------------- Iteration Info ------------------------------------")
            println("Iter |   LB                              UB                             gap")
        end
        @printf("%3d  |   %5.3g                         %5.3g                              %1.3f%s\n", i, LB, UB, gap, "%")
        if UB-LB ≤ 1e-2 * UB || i > max_iter
            return Dict(:solHistory => sddipResult, 
                            :solution => solCollection[1, 1].stageSolution, 
                            :gapHistory => gapList) 
        end

        for t = reverse(2:T)
            for k in 1:M 
                c = [0, zeros(Float64,binaryInfo.n)]
                for j in keys(Ω[t])
                    @info "$t, $k, $j"
                    backwardInfo = backwardInfoList[t]
                    backward_Constraint_Modification!(backwardInfo, 
                                                    stageDataList[t], 
                                                    Ω[t][j].d, 
                                                    solCollection[t-1,k].stageSolution;
                                                    binaryInfo = binaryInfo             
                                                    )
                    if Enhanced_Cut
                        forward_modify_constraints!(forwardInfoList[t], 
                                                        stageDataList[t], 
                                                        Ω[t][j].d, 
                                                        solCollection[t-1,k].stageSolution, 
                                                        binaryInfo = binaryInfo
                                                        )
                        optimize!(forwardInfoList[t].model); f_star_value = JuMP.objective_value(forwardInfoList[t].model);
                        # L̃ = [.8 for i in 1:n] - .5 *  L̂;
                        L̃ = solCollection[t-1,k].stageSolution .* .6 .+ .2;
                        # L̃ = [.5 for i in 1:binaryInfo.n]
                    else 
                        L̃ = nothing;
                        f_star_value = nothing;
                    end
                    λ_value = .3; Output = 0; Output_Gap = true; Enhanced_Cut = true; threshold = 1.0; 
                    levelSetMethodParam = LevelSetMethodParam(0.95, λ_value, threshold, 1e14, 3e3, Output, Output_Gap)
                    c = c + probList[t][j] * LevelSetMethod_optimization!(  backwardInfo, f_star_value; 
                                                                                levelSetMethodParam = levelSetMethodParam, 
                                                                                stageData = stageDataList[t], 
                                                                                         L̂ = solCollection[t-1,k].stageSolution, 
                                                                                            ϵ = 1e-4, 
                                                                                                Enhanced_Cut = Enhanced_Cut, 
                                                                                                    binaryInfo = binaryInfo, 
                                                                                                    L̃ = L̃) 
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



using JLD2, FileIO, DataFrames
result_enhanced = copy(Dict(:solHistory => sddipResult, :solution => _LB, :gapHistory => gapList) )
cut_enhanced = copy(cutCollection)

@save "runtests_small2_enhanced.jld2" result_enhanced cut_enhanced
# # @load "runtests_small2_enhanced.jld2" result_enhanced cut_enhanced

# result_LC = copy(Dict(:solHistory => sddipResult, :solution => _LB, :gapHistory => gapList) )
# cut_LC = copy(cutCollection)

# @save "runtests_small2_LC.jld2" result_LC cut_LC
# # @load "runtests_small2_LC.jld2" result_LC cut_LC



# using DataFrames
# using Latexify
# df = DataFrame(A = 'x':'z', B = ["M", "F", "F"])
# latexify(df; env=:table, latex=false)




