
function SDDiP_algorithm( ; scenarioTree::ScenarioTree = scenarioTree, 
                                indexSets::IndexSets = indexSets, 
                                    paramDemand::ParamDemand = paramDemand, 
                                        paramOPF::ParamOPF = paramOPF, 
                                            initialStageDecision::Dict{Symbol, Dict{Int64, Float64}} = initialStageDecision,
                                                Output_Gap::Bool = true, max_iter::Int64 = 100, TimeLimit::Float64 = 1e3, cutSelection::String = "LC", δ::Float64 = 50., numScenarios::Int64 = 2, OPT::Float64 = Inf,
                                                    ε::Real = 0.125, MaxIter::Int64 = 20, ℓ::Float64 = 0.0
                        )
    ## d: x dim
    initial = now(); i = 1; LB = - Inf; UB = Inf; 
    iter_time = 0; total_Time = 0; t0 = 0.0; LMiter = 0; LM_iter = 0; gap = 100.0; gapString = "100%";

    col_names = [:Iter, :LB, :OPT, :UB, :gap, :time, :LM_iter, :Time]; # needs to be a vector Symbols
    col_types = [Int64, Float64, Union{Float64,Nothing}, Float64, String, Float64, Int64, Float64];
    named_tuple = (; zip(col_names, type[] for type in col_types )...);
    sddipResult = DataFrame(named_tuple); # 0×7 DataFrame
    gapList = [];
    κ = Dict{Int64, Int64}()
    for g in indexSets.G
        if paramOPF.smax[g] ≥ ε
            κ[g] = floor(Int, log2(paramOPF.smax[g] / ε)) + 1
        else
            κ[g] = 1
        end
    end
   
    forwardInfoList = Dict{Int, Model}();
    backwardInfoList = Dict{Int, Model}();
    for t in 1:indexSets.T 
        forwardInfoList[t] = forwardModel!(indexSets = indexSets, 
                                                paramDemand = paramDemand, 
                                                    paramOPF = paramOPF, 
                                                        stageRealization = scenarioTree.tree[t], 
                                                                outputFlag = 0, κ = κ, ε = ε, 
                                                                        mipGap = 1e-4, θ_bound = 0.0, timelimit = 20
                                                );
        backwardInfoList[t] = backwardModel!(indexSets = indexSets, 
                                                paramDemand = paramDemand, 
                                                    paramOPF = paramOPF, 
                                                        stageRealization = scenarioTree.tree[t], 
                                                                outputFlag = 0, κ = κ, ε = ε, 
                                                                        mipGap = 1e-4, tightness = tightness, θ_bound = 0.0, timelimit = 5
                                                );
    end 
    stageDecision = Dict();
    λ₀ = 0.0; λ₁ = Dict{Symbol, Dict{Int64, Float64}}();
    solCollection = Dict();  # to store every iteration results

    ####################################################### Main Loop ###########################################################
    while true
        t0 = now(); 
        Ξ̃ = sample_scenarios(; scenarioTree = scenarioTree, numScenarios = numScenarios);
        u = Dict{Int64, Float64}();  # to store the value of each scenario
        ####################################################### Forward Steps ###########################################################
        for ω in keys(Ξ̃)
            stageDecision[:s] = Dict{Int64, Float64}(g => initialStageDecision[:s][g] for g in indexSets.G); stageDecision[:y] = Dict{Int64, Float64}(g => initialStageDecision[:y][g] for g in indexSets.G);  
            # stageDecision[:s] = Dict{Int64, Float64}(g => 0.0 for g in indexSets.G); stageDecision[:y] = Dict{Int64, Float64}(g => 0.0 for g in indexSets.G);  
            for t in 1:indexSets.T
                forwardModification!(model = forwardInfoList[t], randomVariables = Ξ̃[ω][t], paramOPF = paramOPF, indexSets = indexSets, stageDecision = stageDecision, paramDemand = paramDemand);
                optimize!(forwardInfoList[t]); st = termination_status(forwardInfoList[t]);
                stageDecision = Dict();
                stageDecision[:y] = Dict{Int64, Float64}(g => round(JuMP.value(forwardInfoList[t][:y][g]), digits = 6) for g in indexSets.G);
                stageDecision[:s] = Dict{Int64, Float64}(g => JuMP.value(forwardInfoList[t][:s][g]) for g in indexSets.G);
                stageDecision[:λ] = Dict{Int64, Dict{Int64, Float64}}(g => Dict(i => round(JuMP.value(forwardInfoList[t][:λ][g, i]), digits = 6) for i in 1:κ[g]) for g in indexSets.G);
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
        gap = round((UB-LB)/UB * 100 ,digits = 2); gapString = string(gap,"%"); push!(sddipResult, [i, LB, OPT, UB, gapString, iter_time, LM_iter, total_Time]); push!(gapList, gap);
        if i == 1
            println("---------------------------------- Iteration Info ------------------------------------")
            println("Iter |   LB                              UB                             gap")
        end
        @printf("%3d  |   %5.3g                         %5.3g                              %1.3f%s\n", i, LB, UB, gap, "%")
        if total_Time > TimeLimit || i > MaxIter # || UB-LB ≤ 1e-2 * UB
            return Dict(:solHistory => sddipResult, 
                            :solution => solCollection[i, 1, 1].stageSolution, 
                                :gapHistory => gapList) 
        end

        ####################################################### Backward Steps ###########################################################
        for t = reverse(2:indexSets.T)
            for ω in [1]#keys(Ξ̃)
                # λ₀ = 0.0; λ₁ = Dict{Symbol, Dict{Int64, Float64}}(); λ₁[:s] = Dict{Int64, Float64}(g => 0.0 for g in indexSets.G); λ₁[:y] = Dict{Int64, Float64}(g => 0.0 for g in indexSets.G);
                for n in keys(scenarioTree.tree[t].nodes)
                    # @info "$t, $k, $j"
                    forwardModification!(model = forwardInfoList[t], randomVariables = scenarioTree.tree[t].nodes[n], stageDecision = solCollection[i, t-1, ω].stageSolution, paramOPF = paramOPF, indexSets = indexSets, paramDemand = paramDemand);
                    optimize!(forwardInfoList[t]); f_star_value = JuMP.objective_value(forwardInfoList[t]);

                    backwardModification!(model = backwardInfoList[t], randomVariables = scenarioTree.tree[t].nodes[n], paramOPF = paramOPF, indexSets = indexSets, paramDemand = paramDemand);

                    (x_interior, levelSetMethodParam, x₀) = setupLevelSetMethod(stageDecision = solCollection[i, t-1, ω].stageSolution, 
                                                                        f_star_value = f_star_value, κ = κ,
                                                                            cutSelection = cutSelection, max_iter = max_iter,
                                                                                Output_Gap = Output_Gap, ℓ = ℓ, λ = .3);
                    # model = backwardInfoList[t]; stageDecision = solCollection[i, t-1, ω].stageSolution;
                    ((λ₀, λ₁), LMiter) = LevelSetMethod_optimization!(levelSetMethodParam = levelSetMethodParam, model = backwardInfoList[t],
                                                                    cutSelection = cutSelection, stageDecision = solCollection[i, t-1, ω].stageSolution, κ = κ,
                                                                        x_interior = x_interior, x₀ = x₀, indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF, δ = δ);
                    # λ₀ + sum(sum(λ₁[:λ][g][k] * solCollection[i, t-1, ω].stageSolution[:λ][g][k] for k in 1:κ[g]) + λ₁[:y][g] * solCollection[i, t-1, ω].stageSolution[:y][g] for g in indexSets.G)
                    # add cut to both backward models and forward models
                    @constraint(forwardInfoList[t-1], forwardInfoList[t-1][:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + sum(sum(λ₁[:λ][g][k] * forwardInfoList[t-1][:λ][g, k] for k in 1:κ[g]) + λ₁[:y][g] * forwardInfoList[t-1][:y][g] for g in indexSets.G));
                    @constraint(backwardInfoList[t-1], backwardInfoList[t-1][:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + sum(sum(λ₁[:λ][g][k] * backwardInfoList[t-1][:λ][g, k] for k in 1:κ[g]) + λ₁[:y][g] * backwardInfoList[t-1][:y][g] for g in indexSets.G));
                end            
                LM_iter += LMiter;     
            end
        end
        LM_iter = floor(Int64, LM_iter/sum(length(scenarioTree.tree[t].nodes) for t in 2:indexSets.T));
        
        t1 = now(); iter_time = (t1 - t0).value/1000; total_Time = (t1 - initial).value/1000; i += 1;
    end
end