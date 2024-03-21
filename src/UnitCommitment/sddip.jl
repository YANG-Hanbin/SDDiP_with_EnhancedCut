
function SDDiP_algorithm( ; scenarioTree::ScenarioTree = scenarioTree, 
                                indexSets::IndexSets = indexSets, 
                                    paramDemand::ParamDemand = paramDemand, 
                                        paramOPF::ParamOPF = paramOPF, 
                                            Output_Gap::Bool = true, max_iter::Int64 = 100, ϵ::Float64 = 1e-3, cutSelection::String = "LC", δ::Float64 = 50.
                        )
    ## d: x dim
    initial = now(); i = 1; LB = - Inf; UB = Inf; OPT = Inf;
    iter_time = 0; total_Time = 0; t0 = 0.0; numScenarios = 5;

    col_names = [:iter, :LB, :OPT, :UB, :gap, :time, :Time]; # needs to be a vector Symbols
    col_types = [Int64, Float64, Union{Float64,Nothing}, Float64, String, Float64, Float64];
    named_tuple = (; zip(col_names, type[] for type in col_types )...);
    sddipResult = DataFrame(named_tuple); # 0×7 DataFrame
    gapList = [];
    # @time extResult = extensive_form(indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF, scenarioTree = scenarioTree, Ξ = Ξ, initialStageDecision = initialStageDecision); OPT = extResult.OPT;
    OPT = 0.0;
    forwardInfoList = Dict{Int, Model}();
    backwardInfoList = Dict{Int, Model}();
    for t in 1:indexSets.T 
        forwardInfoList[t] = forwardModel!(indexSets = indexSets, 
                                                paramDemand = paramDemand, 
                                                    paramOPF = paramOPF, 
                                                        stageRealization = scenarioTree.tree[t], 
                                                                outputFlag = 0, 
                                                                        mipGap = 1e-2, θ_bound = 0.0, timelimit = 20
                                                );
        backwardInfoList[t] = backwardModel!(indexSets = indexSets, 
                                                paramDemand = paramDemand, 
                                                    paramOPF = paramOPF, 
                                                        stageRealization = scenarioTree.tree[t], 
                                                                outputFlag = 0, 
                                                                        mipGap = 1e-5, tightness = true, θ_bound = 0.0, timelimit = 20
                                                );
    end 
    stageDecision = Dict{Symbol, Dict{Int64, Float64}}(); stageDecision[:s] = Dict{Int64, Float64}(); stageDecision[:y] = Dict{Int64, Float64}();
    λ₀ = 0.0; λ₁ = Dict{Symbol, Dict{Int64, Float64}}();
    solCollection = Dict();  # to store every iteration results

    ####################################################### Main Loop ###########################################################
    while true
        t0 = now(); 
        Ξ̃ = sample_scenarios(; scenarioTree = scenarioTree, numScenarios = numScenarios);
        u = Dict{Int64, Float64}();  # to store the value of each scenario
        ####################################################### Forward Steps ###########################################################
        for ω in keys(Ξ̃)
            stageDecision = initialStageDecision
            stageDecision[:s] = Dict{Int64, Float64}(g => 0.0 for g in indexSets.G); stageDecision[:y] = Dict{Int64, Float64}(g => 0.0 for g in indexSets.G);  
            for t in 1:indexSets.T
                forwardModification!(model = forwardInfoList[t], randomVariables = Ξ̃[ω][t], paramOPF = paramOPF, indexSets = indexSets, stageDecision = stageDecision);
                optimize!(forwardInfoList[t]);
                stageDecision[:s] = Dict{Int64, Float64}(g => JuMP.value(forwardInfoList[t][:s][g]) for g in indexSets.G);
                stageDecision[:y] = Dict{Int64, Float64}(g => JuMP.value(forwardInfoList[t][:y][g]) for g in indexSets.G);
                solCollection[i, t, ω] = ( stageSolution = deepcopy(stageDecision), 
                                        stageValue = JuMP.objective_value(forwardInfoList[t]) - sum(JuMP.value.(forwardInfoList[t][:θ])), 
                                        OPT = JuMP.objective_value(forwardInfoList[t]));

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
        if UB-LB ≤ 1e-2 * UB || total_Time > 6000
            return Dict(:solHistory => sddipResult, 
                            :solution => solCollection[i, 1, 1].stageSolution, 
                                :gapHistory => gapList) 
        end

        ####################################################### Backward Steps ###########################################################
        for t = reverse(2:indexSets.T)
            for ω in keys(Ξ̃)
                # λ₀ = 0.0; λ₁ = Dict{Symbol, Dict{Int64, Float64}}(); λ₁[:s] = Dict{Int64, Float64}(g => 0.0 for g in indexSets.G); λ₁[:y] = Dict{Int64, Float64}(g => 0.0 for g in indexSets.G);
                for n in keys(scenarioTree.tree[t].nodes)
                    # @info "$t, $k, $j"
                    forwardModification!(model = forwardInfoList[t], randomVariables = scenarioTree.tree[t].nodes[n], stageDecision = solCollection[i, t-1, ω].stageSolution, paramOPF = paramOPF, indexSets = indexSets);
                    optimize!(forwardInfoList[t]); f_star_value = JuMP.objective_value(forwardInfoList[t]);

                    backwardModification!(model = backwardInfoList[t], randomVariables = scenarioTree.tree[t].nodes[n], paramOPF = paramOPF, indexSets = indexSets);

                    (x_interior, levelSetMethodParam, x₀) = setupLevelSetMethod(stageDecision = solCollection[i, t-1, ω].stageSolution, 
                                                                        f_star_value = f_star_value, 
                                                                            cutSelection = cutSelection, max_iter = max_iter,
                                                                                Output_Gap = Output_Gap, ℓ = .0, λ = .1 );
                    # model = backwardInfoList[t]; stageDecision = solCollection[i, t-1, ω].stageSolution;
                    (λ₀, λ₁) = LevelSetMethod_optimization!(levelSetMethodParam = levelSetMethodParam, model = backwardInfoList[t],
                                                            cutSelection = cutSelection,
                                                                stageDecision = solCollection[i, t-1, ω].stageSolution,
                                                                        x_interior = x_interior, x₀ = x₀, tightness = false, indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF, δ =δ);
                    # add cut to both backward models and forward models
                    @constraint(forwardInfoList[t-1], forwardInfoList[t-1][:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + sum(λ₁[:s][g] * forwardInfoList[t-1][:s][g] + λ₁[:y][g] * forwardInfoList[t-1][:y][g] for g in indexSets.G));
                    @constraint(backwardInfoList[t-1], backwardInfoList[t-1][:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + sum(λ₁[:s][g] * backwardInfoList[t-1][:s][g] + λ₁[:y][g] * backwardInfoList[t-1][:y][g] for g in indexSets.G));
                end            
            end
        end
        
        t1 = now(); iter_time = (t1 - t0).value/1000; total_Time = (t1 - initial).value/1000; i += 1;
    end
end

# λ₀ + sum(λ₁[:s][g] * solCollection[i, t-1, ω].stageSolution[:s][g] + λ₁[:y][g] * solCollection[i, t-1, ω].stageSolution[:y][g] for g in indexSets.G)

# t = 2; i = 1; ω = 3; j = 2;
# solCollection[i, t, ω].stageSolution[:s]
# solCollection[i, t, ω].stageSolution[:s] == solCollection[j, t, ω].stageSolution[:s]
# solCollection[i, t, ω].stageSolution[:y] == solCollection[j, t, ω].stageSolution[:y]