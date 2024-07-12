
function SDDiP_algorithm( ; scenarioTree::ScenarioTree = scenarioTree, 
                                indexSets::IndexSets = indexSets, 
                                    paramDemand::ParamDemand = paramDemand, 
                                        paramOPF::ParamOPF = paramOPF, 
                                            initialStageDecision::Dict{Symbol, Dict{Int64, Float64}} = initialStageDecision, numScenarios::Real = 2,
                                            Output_Gap::Bool = true, max_iter::Int64 = 100, TimeLimit::Float64 = 1e3, cutSelection::String = "LC", δ::Float64 = 1., tightness::Bool = true, OPT::Float64 = Inf
                        )
    ## d: x dim
    initial = now(); i = 1; LB = - Inf; UB = Inf; 
    iter_time = 0; total_Time = 0; t0 = 0.0; 

    col_names = [:iter, :LB, :OPT, :UB, :gap, :time, :Time]; # needs to be a vector Symbols
    col_types = [Int64, Float64, Union{Float64,Nothing}, Float64, String, Float64, Float64];
    named_tuple = (; zip(col_names, type[] for type in col_types )...);
    sddipResult = DataFrame(named_tuple); # 0×7 DataFrame
    gapList = [];
    
    forwardInfoList = Dict{Int, Model}();
    backwardInfoList = Dict{Int, Model}();

    StateVarList = Dict(); 
    for t in 1:indexSets.T 
        forwardInfoList[t] = forwardModel!(indexSets = indexSets, 
                                                paramDemand = paramDemand, 
                                                    paramOPF = paramOPF, 
                                                        stageRealization = scenarioTree.tree[t], 
                                                                outputFlag = 0, timelimit = 20, max_sur = max_iter,
                                                                        mipGap = 1e-4, θ_bound = 0.0);
        var = Dict{Symbol, Dict{Int, VariableRef}}(:s => Dict(g => forwardInfoList[t][:s][g] for g in indexSets.G), :y => Dict(g => forwardInfoList[t][:y][g] for g in indexSets.G)); 
        sur = Dict{Int, Dict{Int, Dict{Symbol, Any}}}(g => Dict(1 => Dict(:lb => 0., :ub => paramOPF.smax[g], :var =>forwardInfoList[t][:sur][g, 1])) for g in indexSets.G); 
        leaf = Dict{Int, Vector{Int64}}(g => [1] for g in indexSets.G);
        StateVarList[t] = StateVar(var, sur, leaf);

        backwardInfoList[t] = backwardModel!(indexSets = indexSets, 
                                                paramDemand = paramDemand, 
                                                    paramOPF = paramOPF, 
                                                        stageRealization = scenarioTree.tree[t], 
                                                                outputFlag = 0, timelimit = 20, max_sur = max_iter,
                                                                        mipGap = 1e-4, tightness = tightness, θ_bound = 0.0);
        # var = Dict{Symbol, Dict{Int, VariableRef}}(:s => Dict(g => backwardInfoList[t][:s][g] for g in indexSets.G), :y => Dict(g => backwardInfoList[t][:y][g] for g in indexSets.G)); 
        # sur = Dict{Int, Dict{Int, Dict{Symbol, Any}}}(g => Dict(1 => Dict(:lb => 0., :ub => paramOPF.smax[g], :var =>backwardInfoList[t][:sur][g, 1])) for g in indexSets.G); 
        # backwardStateVarList[t] = StateVar(var, sur, leaf);
    end 

    stageDecision = Dict{Symbol, Dict{Int64, Any}}(); stageDecision[:s] = Dict{Int64, Float64}(); stageDecision[:y] = Dict{Int64, Float64}(); stageDecision[:sur] = Dict{Int64, Dict{Int64, Float64}}();
    λ₀ = 0.0; λ₁ = Dict{Symbol, Dict{Int64, Any}}();
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
            for g in indexSets.G stageDecision[:sur][g] = Dict(1 => 1.) end; # augmented state variables
            for t in 1:indexSets.T
                forwardModification!(model = forwardInfoList[t], randomVariables = Ξ̃[ω][t], paramOPF = paramOPF, indexSets = indexSets, stageDecision = stageDecision, paramDemand = paramDemand);
                optimize!(forwardInfoList[t]);
                stageDecision[:s] = Dict{Int64, Float64}(g => JuMP.value(forwardInfoList[t][:s][g]) for g in indexSets.G);
                stageDecision[:y] = Dict{Int64, Float64}(g => JuMP.value(forwardInfoList[t][:y][g]) for g in indexSets.G);
                for g in indexSets.G 
                    stageDecision[:sur][g] = Dict{Int64, Float64}()
                    for k in StateVarList[t].leaf[g]
                        stageDecision[:sur][g][k] = JuMP.value(forwardInfoList[t][:sur][g, k])
                    end
                end
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
        if UB-LB ≤ 1e-2 * UB || total_Time > TimeLimit || i > max_iter
            return Dict(:solHistory => sddipResult, 
                            :solution => solCollection[i, 1, 1].stageSolution, 
                                :gapHistory => gapList) 
        end

        ####################################################### Backward Steps ###########################################################
        for t = reverse(2:indexSets.T)
            for ω in [1]#keys(Ξ̃)
                # λ₀ = 0.0; λ₁ = Dict{Symbol, Dict{Int64, Float64}}(); λ₁[:s] = Dict{Int64, Float64}(g => 0.0 for g in indexSets.G); λ₁[:y] = Dict{Int64, Float64}(g => 0.0 for g in indexSets.G);
                for n in keys(scenarioTree.tree[t].nodes)
                    f_star_value = solCollection[i, t, ω].OPT; backwardModification!(model = backwardInfoList[t], randomVariables = scenarioTree.tree[t].nodes[n], paramOPF = paramOPF, indexSets = indexSets, paramDemand = paramDemand);
                    
                    (x_interior, levelSetMethodParam, x₀) = setupLevelSetMethod(stageDecision = solCollection[i, t-1, ω].stageSolution, f_star_value = f_star_value, cutSelection = "LC", max_iter = 150, paramOPF = paramOPF,
                                                                            Output_Gap = Output_Gap, ℓ = .0, λ = .1 );
                    (λ₀, λ₁) = LevelSetMethod_optimization!(levelSetMethodParam = levelSetMethodParam, model = backwardInfoList[t], cutSelection = "LC", indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF,
                                                            stageDecision = solCollection[i, t-1, ω].stageSolution,
                                                                    x_interior = nothing, x₀ = x₀, tightness = tightness, δ = δ);
                    f_star_value = λ₀ + sum(λ₁[:s][g] * solCollection[i, t-1, ω].stageSolution[:s][g] + 
                                                λ₁[:y][g] * solCollection[i, t-1, ω].stageSolution[:y][g] + 
                                                    sum(λ₁[:sur][g][k] * solCollection[i, t-1, ω].stageSolution[:sur][g][k] for k in keys(solCollection[i, t-1, ω].stageSolution[:sur][g])) for g in indexSets.G
                                                    );
                    # add cut to both backward models and forward models
                    @constraint(forwardInfoList[t-1], forwardInfoList[t-1][:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + 
                                                    sum(λ₁[:s][g] * forwardInfoList[t-1][:s][g] + λ₁[:y][g] * forwardInfoList[t-1][:y][g] for g in indexSets.G) + 
                                                        sum(sum(λ₁[:sur][g][k] * forwardInfoList[t-1][:sur][g, k] for k in keys(solCollection[i, t-1, ω].stageSolution[:sur][g])) for g in indexSets.G));
                    @constraint(backwardInfoList[t-1], backwardInfoList[t-1][:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + 
                                                    sum(λ₁[:s][g] * backwardInfoList[t-1][:s][g] + λ₁[:y][g] * backwardInfoList[t-1][:y][g] for g in indexSets.G) + 
                                                        sum(sum(λ₁[:sur][g][k] * backwardInfoList[t-1][:sur][g, k] for k in keys(solCollection[i, t-1, ω].stageSolution[:sur][g])) for g in indexSets.G));
                    ## using enhancement cuts under conditions
                    if cutSelection != "LC"  && gap ≥ 5.
                        (x_interior, levelSetMethodParam, x₀) = setupLevelSetMethod(stageDecision = solCollection[i, t-1, ω].stageSolution, f_star_value = f_star_value, cutSelection = cutSelection, max_iter = 150, paramOPF = paramOPF,
                                                                                    Output_Gap = Output_Gap, ℓ = .0, λ = .1 );
                        # model = backwardInfoList[t]; stageDecision = solCollection[i, t-1, ω].stageSolution;
                        (λ₀, λ₁) = LevelSetMethod_optimization!(levelSetMethodParam = levelSetMethodParam, model = backwardInfoList[t], cutSelection = cutSelection, indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF,
                                                                    stageDecision = solCollection[i, t-1, ω].stageSolution,
                                                                            x_interior = x_interior, x₀ = x₀, tightness = tightness, δ = δ);
                        # add cut to both backward models and forward models
                        @constraint(forwardInfoList[t-1], forwardInfoList[t-1][:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + 
                                                        sum(λ₁[:s][g] * forwardInfoList[t-1][:s][g] + λ₁[:y][g] * forwardInfoList[t-1][:y][g] for g in indexSets.G) + 
                                                            sum(sum(λ₁[:sur][g][k] * forwardInfoList[t-1][:sur][g, k] for k in keys(solCollection[i, t-1, ω].stageSolution[:sur][g])) for g in indexSets.G));
                        @constraint(backwardInfoList[t-1], backwardInfoList[t-1][:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + 
                                                        sum(λ₁[:s][g] * backwardInfoList[t-1][:s][g] + λ₁[:y][g] * backwardInfoList[t-1][:y][g] for g in indexSets.G) + 
                                                            sum(sum(λ₁[:sur][g][k] * backwardInfoList[t-1][:sur][g, k] for k in keys(solCollection[i, t-1, ω].stageSolution[:sur][g])) for g in indexSets.G));
                    end

                end            
            end
        end
        
        ####################################################### Backward Steps ###########################################################
        for t in 1:indexSets.T-1 
            for ω in  [1]#keys(Ξ̃)
                dev = Dict()
                for g in indexSets.G 
                    if solCollection[i, t, ω].stageSolution[:y][g] > .5
                        k = maximum([k for (k, v) in solCollection[i, t, ω].stageSolution[:sur][g] if v == maximum(values(solCollection[i, t, ω].stageSolution[:sur][g]))])
                        info = StateVarList[t].sur[g][k]
                        dev[g] = round(minimum([(info[:ub] - solCollection[i, t, ω].stageSolution[:s][g])/(info[:ub] - info[:lb] + 1e-6), (solCollection[i, t, ω].stageSolution[:s][g] - info[:lb])/(info[:ub] - info[:lb] + 1e-6)]), digits = 5)
                    end
                end
                g = [k for (k, v) in dev if v == maximum(values(dev))][1]
                if dev[g] >= 1e-6
                    # find the active leaf node 
                    keys_with_value_1 = maximum([k for (k, v) in solCollection[i, t, ω].stageSolution[:sur][g] if v == 1])
                    # find the lb and ub of this leaf node 
                    (lb, ub) = StateVarList[t].sur[g][keys_with_value_1][:lb], StateVarList[t].sur[g][keys_with_value_1][:ub]; med = solCollection[i, t, ω].stageSolution[:s][g]; # solCollection[i, t, ω].stageSolution[:s][g];# (lb + ub)/2; #round(solCollection[i, t, ω].stageSolution[:s][g], digits = 3); 
                    # create two new leaf nodes, and update their info (lb, ub)
                    left = length(StateVarList[t].sur[g]) + 1; right = length(StateVarList[t].sur[g]) + 2;
                    StateVarList[t].sur[g][left] = Dict(:lb => lb, :ub => med, :var => forwardInfoList[t][:sur][g, left])
                    StateVarList[t].sur[g][right] =  Dict(:lb => med, :ub => ub, :var => forwardInfoList[t][:sur][g, right])
                    # pop and push new leaf nodes
                    deleteat!(StateVarList[t].leaf[g], findall(x -> x == keys_with_value_1, StateVarList[t].leaf[g])); push!(StateVarList[t].leaf[g], left); push!(StateVarList[t].leaf[g], right);

                    # add logic constraints
                    ## for forward models
                    ### Parent-Child relationship
                    @constraint(forwardInfoList[t], forwardInfoList[t][:sur][g, left] + forwardInfoList[t][:sur][g, right] == forwardInfoList[t][:sur][g, keys_with_value_1])
                    ### bounding constraints
                    # if i > 1
                    #     for g in indexSets.G
                    #         delete(forwardInfoList[t], forwardInfoList[t][:bounding][g]);
                    #         delete(backwardInfoList[t], backwardInfoList[t][:bounding][g]);
                    #     end
                    #     unregister(forwardInfoList[t], :Nonanticipativity); unregister(backwardInfoList[t], :Nonanticipativity);
                    # end
                    @constraint(forwardInfoList[t], forwardInfoList[t][:s][g] ≥ sum(StateVarList[t].sur[g][k][:lb] * forwardInfoList[t][:sur][g, k] for k in StateVarList[t].leaf[g]))
                    @constraint(forwardInfoList[t], forwardInfoList[t][:s][g] ≤ sum(StateVarList[t].sur[g][k][:ub] * forwardInfoList[t][:sur][g, k] for k in StateVarList[t].leaf[g]))
                    # @constraint(forwardInfoList[t], forwardInfoList[t][:s][g] ≥ lb -  paramOPF.smax[g] * (1 - forwardInfoList[t][:sur][g, left]) );
                    # @constraint(forwardInfoList[t], forwardInfoList[t][:s][g] ≤ med + paramOPF.smax[g] * (1 - forwardInfoList[t][:sur][g, left]) );
                    # @constraint(forwardInfoList[t], forwardInfoList[t][:s][g] ≥ med -  paramOPF.smax[g] * (1 - forwardInfoList[t][:sur][g, right]));
                    # @constraint(forwardInfoList[t], forwardInfoList[t][:s][g] ≤ ub + paramOPF.smax[g] * (1 - forwardInfoList[t][:sur][g, right]) );
                    ## for backward models
                    ### Parent-Child relationship
                    @constraint(backwardInfoList[t+1], backwardInfoList[t+1][:sur_copy][g, left] + backwardInfoList[t+1][:sur_copy][g, right] == backwardInfoList[t+1][:sur_copy][g, keys_with_value_1]);
                    @constraint(backwardInfoList[t], backwardInfoList[t][:sur][g, left] + backwardInfoList[t][:sur][g, right] == backwardInfoList[t][:sur][g, keys_with_value_1]);
                    ### bounding constraints
                    @constraint(backwardInfoList[t+1], backwardInfoList[t+1][:s][g] ≥ sum(StateVarList[t].sur[g][k][:lb] * backwardInfoList[t+1][:sur_copy][g, k] for k in StateVarList[t].leaf[g]));
                    @constraint(backwardInfoList[t+1], backwardInfoList[t+1][:s][g] ≤ sum(StateVarList[t].sur[g][k][:ub] * backwardInfoList[t+1][:sur_copy][g, k] for k in StateVarList[t].leaf[g]));
                    @constraint(backwardInfoList[t], backwardInfoList[t][:s][g] ≥ sum(StateVarList[t].sur[g][k][:lb] * backwardInfoList[t][:sur][g, k] for k in StateVarList[t].leaf[g]));
                    @constraint(backwardInfoList[t], backwardInfoList[t][:s][g] ≤ sum(StateVarList[t].sur[g][k][:ub] * backwardInfoList[t][:sur][g, k] for k in StateVarList[t].leaf[g]));
                    # @constraint(backwardInfoList[t+1], backwardInfoList[t+1][:s][g] ≥ lb -  paramOPF.smax[g] * (1 - backwardInfoList[t+1][:sur][g, left]) );
                    # @constraint(backwardInfoList[t+1], backwardInfoList[t+1][:s][g] ≤ med + paramOPF.smax[g] * (1 - backwardInfoList[t+1][:sur][g, left]) );
                    # @constraint(backwardInfoList[t+1], backwardInfoList[t+1][:s][g] ≥ med -  paramOPF.smax[g] * (1 - backwardInfoList[t+1][:sur][g, right]));
                    # @constraint(backwardInfoList[t+1], backwardInfoList[t+1][:s][g] ≤ ub + paramOPF.smax[g] * (1 - backwardInfoList[t+1][:sur][g, right]) );
                end
            end
        end
        t1 = now(); iter_time = (t1 - t0).value/1000; total_Time = (t1 - initial).value/1000; i += 1; 

    end
end