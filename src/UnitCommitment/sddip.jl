
function SDDiP_algorithm( ; scenarioTree::ScenarioTree = scenarioTree, 
                                indexSets::IndexSets = indexSets, 
                                    paramDemand::ParamDemand = paramDemand, 
                                        paramOPF::ParamOPF = paramOPF, 
                                            Output_Gap::Bool = true, max_iter::Int64 = 100, ϵ::Float64 = 1e-3, cutSelection::String = "LC"
                        )
    ## d: x dim
    initial = now(); i = 1; LB = - Inf; UB = Inf; OPT = Inf;
    iter_time = 0; total_Time = 0; t0 = 0.0; numScenarios = 5;

    col_names = [:iter, :LB, :OPT, :UB, :gap, :time, :Time]; # needs to be a vector Symbols
    col_types = [Int64, Float64, Union{Float64,Nothing}, Float64, String, Float64, Float64];
    named_tuple = (; zip(col_names, type[] for type in col_types )...);
    sddipResult = DataFrame(named_tuple); # 0×7 DataFrame
    gapList = [];
    # @time gurobiResult = gurobiOptimize!(indexSets, paramDemand, paramOPF, Ξ; β = β);  OPT = gurobiResult.OPT;

    forwardInfoList = Dict{Int, Model}();
    backwardInfoList = Dict{Int, Model}();
    for t in 1:indexSets.T 
        forwardInfoList[t] = forwardModel!(indexSets = indexSets, 
                                                paramDemand = paramDemand, 
                                                    paramOPF = paramOPF, 
                                                        stageRealization = scenarioTree.tree[t], 
                                                                outputFlag = 0, 
                                                                        mipGap = 1e-2
                                                );
        backwardInfoList[t] = backwardModel!(indexSets = indexSets, 
                                                paramDemand = paramDemand, 
                                                    paramOPF = paramOPF, 
                                                        stageRealization = scenarioTree.tree[t], 
                                                                outputFlag = 0, 
                                                                        mipGap = 1e-2, tightness = true
                                                );
    end 
    stageDecision = Dict{Symbol, Dict{Int64, Float64}}(); stageDecision[:s] = Dict{Int64, Float64}(); stageDecision[:y] = Dict{Int64, Float64}();
    λ₀ = 0.0; λ₁ = Dict{Symbol, Dict{Int64, Float64}}();
    solCollection = Dict();  # to store every iteration results
    while true
        t0 = now(); 
        Ξ = sample_scenarios(; scenarioTree = scenarioTree, numScenarios = numScenarios);
        u = Dict{Int64, Float64}();  # to store the value of each scenario
        ####################################################### Forward Steps ###########################################################
        for ω in keys(Ξ)
            stageDecision[:s] = Dict{Int64, Float64}(g => 0.0 for g in indexSets.G); stageDecision[:y] = Dict{Int64, Float64}(g => 0.0 for g in indexSets.G);  
            for t in 1:indexSets.T
                model = forwardInfoList[t];
                randomVariables = Ξ[ω][t];
                forwardModification!(model = model, randomVariables = randomVariables, paramOPF = paramOPF, indexSets = indexSets, stageDecision = stageDecision);
                optimize!(model);
                stageDecision[:s] = Dict{Int64, Float64}(g => JuMP.value(model[:s][g]) for g in indexSets.G);
                stageDecision[:y] = Dict{Int64, Float64}(g => round(JuMP.value(model[:y][g]), digits = 1) for g in indexSets.G);
                solCollection[i, t, ω] = ( stageSolution = stageDecision, 
                                        # stageValue = sum(paramOPF.slope[g] * JuMP.value(model[:s][g]) + paramOPF.intercept[g] * JuMP.value(model[:y][g]) + paramOPF.C_start[g] * JuMP.value(model[:v][g]) + paramOPF.C_down[g] * JuMP.value(model[:w][g]) for g in indexSets.G) + sum(paramDemand.w[d] * (1 - JuMP.value(model[:x][d])) for d in indexSets.D),
                                        stageValue = JuMP.objective_value(model) - sum(JuMP.value.(model[:θ])), 
                                        OPT = JuMP.objective_value(model));

            end
            u[ω] = sum(solCollection[i, t, ω].stageValue for t in 1:indexSets.T);
        end
        
        ####################################################### To Record Info ###########################################################
        LB = solCollection[i, 1, 1].OPT;
        μ̄ = mean(values(u));
        σ̂² = Statistics.var(values(u));
        UB = μ̄ + 1.96 * sqrt(σ̂²/numScenarios); # minimum([μ̄ + 1.96 * sqrt(σ̂²/numScenarios), UB]);
        gap = round((UB-LB)/UB * 100 ,digits = 2); gapString = string(gap,"%"); push!(sddipResult, [i, LB, OPT, UB, gapString, iter_time, total_Time]); push!(gapList, gap);
        if i == 1
            println("---------------------------------- Iteration Info ------------------------------------")
            println("Iter |   LB                              UB                             gap")
        end
        @printf("%3d  |   %5.3g                         %5.3g                              %1.3f%s\n", i, LB, UB, gap, "%")
        if UB-LB ≤ 1e-2 * OPT || i > max_iter
            return Dict(:solHistory => sddipResult, 
                            :solution => solCollection[i, 1, 1].stageSolution, 
                                :gapHistory => gapList) 
        end

        ####################################################### Backward Steps ###########################################################
        for t = reverse(2:indexSets.T)
            for ω in keys(Ξ)
                # λ₀ = 0.0; λ₁ = Dict{Symbol, Dict{Int64, Float64}}(); λ₁[:s] = Dict{Int64, Float64}(g => 0.0 for g in indexSets.G); λ₁[:y] = Dict{Int64, Float64}(g => 0.0 for g in indexSets.G);
                for n in keys(scenarioTree.tree[t].nodes)
                    # @info "$t, $k, $j"
                    forwardModification!(model = forwardInfoList[t], randomVariables = scenarioTree.tree[t].nodes[n], paramOPF = paramOPF, indexSets = indexSets, stageDecision = solCollection[i, t-1, ω].stageSolution);
                    optimize!(forwardInfoList[t]); f_star_value = JuMP.objective_value(forwardInfoList[t]);

                    backwardModification!(; model = backwardInfoList[t], randomVariables = scenarioTree.tree[t].nodes[n], paramOPF = paramOPF, indexSets = indexSets);

                    (x_interior, levelSetMethodParam, x₀) = setupLevelSetMethod(; stageDecision = solCollection[i, t-1, ω].stageSolution, 
                                                                        f_star_value = f_star_value, 
                                                                            cutSelection = cutSelection, 
                                                                                Output_Gap = Output_Gap, ℓ = .0, λ = .1 )
                   
                   
                    model = backwardInfoList[t];
                    stageDecision = solCollection[i, t-1, ω].stageSolution

                    (λ₀, λ₁) = LevelSetMethod_optimization!(levelSetMethodParam = levelSetMethodParam, model = model,
                                                            cutSelection = cutSelection,
                                                                stageDecision = stageDecision,
                                                                        x_interior = x_interior, x₀ = x₀, tightness = false, indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF, ϵ = ϵ)
                    # add cut to both backward models and forward models
                    @constraint(forwardInfoList[t-1], forwardInfoList[t-1][:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + sum(λ₁[:s][g] * forwardInfoList[t-1][:s][g] + λ₁[:y][g] * forwardInfoList[t-1][:y][g] for g in indexSets.G) )
                    @constraint(backwardInfoList[t-1], backwardInfoList[t-1][:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + sum(λ₁[:s][g] * backwardInfoList[t-1][:s][g] + λ₁[:y][g] * backwardInfoList[t-1][:y][g] for g in indexSets.G) )
                end            
            end
        end

        t1 = now(); iter_time = (t1 - t0).value/1000; total_Time = (t1 - initial).value/1000; i += 1;
    end

end
