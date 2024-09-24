"""
    SDDiP_algorithm

# Arguments

    Output_Gap   : whether output level method info
    MaxIter      : maximum number of SDDiP iterations
    max_iter     : maximum number of level method iterations
    Timelimit    : maximum time for SDDiP
    δ            : enhancement parameter for the Lagrangian cuts
    tightness    : whether use tightness cuts (binary or continuous in copy variables)
    OPT          : optimal value of the problem
    numScenarios : number of scenarios sampled in the forward pass
"""

function SDDiP_algorithm( ; scenarioTree::ScenarioTree = scenarioTree, 
                                indexSets::IndexSets = indexSets, 
                                    paramDemand::ParamDemand = paramDemand, 
                                        paramOPF::ParamOPF = paramOPF, 
                                            initialStageDecision::Dict{Symbol, Dict{Int64, Float64}} = initialStageDecision, numScenarios::Real = 2, ℓ::Float64 = .1,
                                            Output_Gap::Bool = true, MaxIter::Int64 = 100, max_iter::Int64 = 100, TimeLimit::Float64 = 1e3, cutSelection::String = "LC", 
                                            δ::Float64 = 1., tightness::Bool = true, OPT::Float64 = Inf, core_point_strategy::String = "Eps"
                        )
    ## d: x dim
    initial = now(); i = 1; LB = - Inf; UB = Inf; 
    iter_time = 0; total_Time = 0; t0 = 0.0; LMiter = 0; LM_iter = 0; gap = 100.0; gapString = "100%"; branchDecision = false;

    col_names = [:Iter, :LB, :OPT, :UB, :gap, :time, :LM_iter, :Time, :Branch]; # needs to be a vector Symbols
    col_types = [Int64, Float64, Union{Float64,Nothing}, Float64, String, Float64, Int64, Float64, Bool]; # needs to be a vector of types
    named_tuple = (; zip(col_names, type[] for type in col_types )...);
    sddipResult = DataFrame(named_tuple); # 0×7 DataFrame
    gapList = [];
    
    @everywhere begin
        forwardInfoList = Dict{Int, Model}();
        backwardInfoList = Dict{Int, Model}();

        StateVarList = Dict(); 
        for t in 1:indexSets.T 
            forwardInfoList[t] = forwardModel!(indexSets = indexSets, 
                                                    paramDemand = paramDemand, 
                                                        paramOPF = paramOPF, 
                                                            stageRealization = scenarioTree.tree[t], 
                                                                    outputFlag = 0, timelimit = 20,
                                                                            mipGap = 1e-3, θ_bound = 0.0);
            var = Dict{Symbol, Dict{Int, VariableRef}}(:s => Dict(g => forwardInfoList[t][:s][g] for g in indexSets.G), :y => Dict(g => forwardInfoList[t][:y][g] for g in indexSets.G)); 
            sur = Dict{Int, Dict{Int, Dict{Symbol, Any}}}(g => Dict(1 => Dict(:lb => 0., :ub => paramOPF.smax[g], :var =>forwardInfoList[t][:sur][g, 1])) for g in indexSets.G); 
            leaf = Dict{Int, Vector{Int64}}(g => [1] for g in indexSets.G);
            StateVarList[t] = StateVar(var, sur, leaf);

            backwardInfoList[t] = backwardModel!(indexSets = indexSets, 
                                                    paramDemand = paramDemand, 
                                                        paramOPF = paramOPF, 
                                                            stageRealization = scenarioTree.tree[t], 
                                                                    outputFlag = 0, timelimit = 20,
                                                                            mipGap = 1e-3, tightness = tightness, θ_bound = 0.0);
            # var = Dict{Symbol, Dict{Int, VariableRef}}(:s => Dict(g => backwardInfoList[t][:s][g] for g in indexSets.G), :y => Dict(g => backwardInfoList[t][:y][g] for g in indexSets.G)); 
            # sur = Dict{Int, Dict{Int, Dict{Symbol, Any}}}(g => Dict(1 => Dict(:lb => 0., :ub => paramOPF.smax[g], :var =>backwardInfoList[t][:sur][g, 1])) for g in indexSets.G); 
            # backwardStateVarList[t] = StateVar(var, sur, leaf);
        end 

        stageDecision = Dict{Symbol, Dict{Int64, Any}}(); stageDecision[:s] = Dict{Int64, Float64}(); stageDecision[:y] = Dict{Int64, Float64}(); stageDecision[:sur] = Dict{Int64, Dict{Int64, Float64}}();
        λ₀ = 0.0; λ₁ = Dict{Symbol, Dict{Int64, Any}}();
        solCollection = Dict();  # to store every iteration results
    end

    ####################################################### Main Loop ###########################################################
    while true
        t0 = now(); 
        Ξ̃ = sample_scenarios(; scenarioTree = scenarioTree, numScenarios = numScenarios);
        u = Dict{Int64, Float64}();  # to store the value of each scenario
        ####################################################### Forward Steps ###########################################################

        forwarPassResult = pmap(forwardPass, values(Ξ̃))

        for j in 1:numScenarios
            ω = collect(keys(Ξ̃))[j]
            for t in 1:indexSets.T
                solCollection[i, t, ω] = forwarPassResult[ω][t]
            end
            u[ω] = sum(solCollection[i, t, ω].stageValue for t in 1:indexSets.T);
        end
        @everywhere solCollection = $solCollection;
        ####################################################### To Record Info ###########################################################
        LB = solCollection[i, 1, 1].OPT;
        μ̄ = mean(values(u));
        σ̂² = Statistics.var(values(u));
        UB = μ̄ + 1.96 * sqrt(σ̂²/numScenarios); # minimum([μ̄ + 1.96 * sqrt(σ̂²/numScenarios), UB]);
        gap = round((UB-LB)/UB * 100 ,digits = 2); gapString = string(gap,"%"); push!(sddipResult, [i, LB, OPT, UB, gapString, iter_time, LM_iter, total_Time, branchDecision]); push!(gapList, gap); 
        if i == 1
            println("----------------------------------------- Iteration Info ------------------------------------------------")
            println("Iter |        LB        |        UB        |       Gap      |      i-time     |    #LM     |     T-Time")
            println("----------------------------------------------------------------------------------------------------------")
        end
        @printf("%4d | %12.2f     | %12.2f     | %9.2f%%     | %9.2f s     | %6d     | %10.2f s     \n", 
                i, LB, UB, gap, iter_time, LM_iter, total_Time); LM_iter = 0;
        if total_Time > TimeLimit || i ≥ MaxIter || UB-LB ≤ 1e-3 * UB  
            return Dict(:solHistory => sddipResult, 
                            :solution => solCollection[i, 1, 1].stageSolution, 
                                :gapHistory => gapList) 
        end

        ####################################################### Backward Steps ###########################################################
        if i ≥ 10                                      ## the first rule:: for branching: current convex envelope is good enough
            for t in reverse(1:indexSets.T-1) 
                for ω in  [1]#keys(Ξ̃)
                    dev = Dict()
                    for g in indexSets.G 
                        if solCollection[i, t, ω].stageSolution[:y][g] > .5
                            k = maximum([k for (k, v) in solCollection[i, t, ω].stageSolution[:sur][g] if v > 0.5]); # maximum(values(solCollection[i, t, ω].stageSolution[:sur][g]))
                            info = StateVarList[t].sur[g][k]
                            dev[g] = round(minimum([(info[:ub] - solCollection[i, t, ω].stageSolution[:s][g])/(info[:ub] - info[:lb] + 1e-6), (solCollection[i, t, ω].stageSolution[:s][g] - info[:lb])/(info[:ub] - info[:lb] + 1e-6)]), digits = 5)
                        end
                    end
                    g = [g for (g, v) in dev if v == maximum(values(dev))][1]
                    if dev[g] ≥ 1e-3                                        ## the second rule: current point is far from being an extreme point
                        branchDecision = true;
                        @everywhere begin
                            t = $t; ω = $ω; g = $g; i = $i; 
                            # find the active leaf node 
                            keys_with_value_1 = maximum([k for (k, v) in solCollection[i, t, ω].stageSolution[:sur][g] if v > 0.5]); ## find the active leaf node: maximum(values(solCollection[i, t, ω].stageSolution[:sur][g]))
                            # find the lb and ub of this leaf node 
                            (lb, ub) = StateVarList[t].sur[g][keys_with_value_1][:lb], StateVarList[t].sur[g][keys_with_value_1][:ub]; med = (lb + ub)/2; # solCollection[i, t, ω].stageSolution[:s][g];# (lb + ub)/2; #round(solCollection[i, t, ω].stageSolution[:s][g], digits = 3); 
                            # create two new leaf nodes, and update their info (lb, ub)
                            left = length(StateVarList[t].sur[g]) + 1; right = length(StateVarList[t].sur[g]) + 2; 
                            forwardInfoList[t][:sur][g, left] = @variable(forwardInfoList[t], base_name = "sur[$g, $left]", binary = true); 
                            forwardInfoList[t][:sur][g, right] = @variable(forwardInfoList[t], base_name = "sur[$g, $right]", binary = true);
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
                            if tightness
                                backwardInfoList[t+1][:sur_copy][g, left] = @variable(backwardInfoList[t+1], base_name = "sur_copy[$g, $left]", binary = true); 
                                backwardInfoList[t+1][:sur_copy][g, right] = @variable(backwardInfoList[t+1], base_name = "sur_copy[$g, $left]", binary = true); 
                                backwardInfoList[t][:sur][g, left] = @variable(backwardInfoList[t], base_name = "sur[$g, $left]", binary = true); 
                                backwardInfoList[t][:sur][g, right] = @variable(backwardInfoList[t], base_name = "sur[$g, $left]", binary = true); 
                            else
                                backwardInfoList[t+1][:sur_copy][g, left] = @variable(backwardInfoList[t+1], base_name = "sur_copy[$g, $left]", lower_bound = 0, upper_bound = 1); 
                                backwardInfoList[t+1][:sur_copy][g, right] = @variable(backwardInfoList[t+1], base_name = "sur_copy[$g, $left]", lower_bound = 0, upper_bound = 1); 
                                backwardInfoList[t][:sur][g, left] = @variable(backwardInfoList[t], base_name = "sur[$g, $left]", lower_bound = 0, upper_bound = 1); 
                                backwardInfoList[t][:sur][g, right] = @variable(backwardInfoList[t], base_name = "sur[$g, $left]", lower_bound = 0, upper_bound = 1); 
                            end
                            @constraint(backwardInfoList[t+1], backwardInfoList[t+1][:sur_copy][g, left] + backwardInfoList[t+1][:sur_copy][g, right] == backwardInfoList[t+1][:sur_copy][g, keys_with_value_1]);
                            @constraint(backwardInfoList[t], backwardInfoList[t][:sur][g, left] + backwardInfoList[t][:sur][g, right] == backwardInfoList[t][:sur][g, keys_with_value_1]);
                            ### bounding constraints
                            @constraint(backwardInfoList[t+1], backwardInfoList[t+1][:s_copy][g] ≥ sum(StateVarList[t].sur[g][k][:lb] * backwardInfoList[t+1][:sur_copy][g, k] for k in StateVarList[t].leaf[g]));
                            @constraint(backwardInfoList[t+1], backwardInfoList[t+1][:s_copy][g] ≤ sum(StateVarList[t].sur[g][k][:ub] * backwardInfoList[t+1][:sur_copy][g, k] for k in StateVarList[t].leaf[g]));
                            @constraint(backwardInfoList[t], backwardInfoList[t][:s][g] ≥ sum(StateVarList[t].sur[g][k][:lb] * backwardInfoList[t][:sur][g, k] for k in StateVarList[t].leaf[g]));
                            @constraint(backwardInfoList[t], backwardInfoList[t][:s][g] ≤ sum(StateVarList[t].sur[g][k][:ub] * backwardInfoList[t][:sur][g, k] for k in StateVarList[t].leaf[g]));
                            # @constraint(backwardInfoList[t+1], backwardInfoList[t+1][:s][g] ≥ lb -  paramOPF.smax[g] * (1 - backwardInfoList[t+1][:sur][g, left]) );
                            # @constraint(backwardInfoList[t+1], backwardInfoList[t+1][:s][g] ≤ med + paramOPF.smax[g] * (1 - backwardInfoList[t+1][:sur][g, left]) );
                            # @constraint(backwardInfoList[t+1], backwardInfoList[t+1][:s][g] ≥ med -  paramOPF.smax[g] * (1 - backwardInfoList[t+1][:sur][g, right]));
                            # @constraint(backwardInfoList[t+1], backwardInfoList[t+1][:s][g] ≤ ub + paramOPF.smax[g] * (1 - backwardInfoList[t+1][:sur][g, right]) );

                            ## update the stageDecision
                            stageDecision = deepcopy(solCollection[i, t, ω].stageSolution);
                            stageDecision[:sur][g] = Dict{Int64, Float64}()
                            for k in StateVarList[t].leaf[g]
                                stageDecision[:sur][g][k] = 0.0
                            end
                            stageDecision[:sur][g][left] = 1.0; 
                            solCollection[i, t, ω] = ( stageSolution = deepcopy(stageDecision), 
                                                        stageValue = solCollection[i, t, ω].stageValue, 
                                                        OPT = solCollection[i, t, ω].OPT
                                                    );
                        end
                    else 
                        branchDecision = false;
                    end            
                end
            end
        end

        ####################################################### Backward Steps ###########################################################
        for t = reverse(2:indexSets.T)
            for ω in [1]#keys(Ξ̃)
                backwardNodeInfoSet = Dict{Int64, Tuple}(); 
                for n in keys(scenarioTree.tree[t].nodes) backwardNodeInfoSet[n] = (i, t, n, ω, cutSelection, core_point_strategy) end
                backwardPassResult = pmap(backwardPass, values(backwardNodeInfoSet));

                for n in keys(scenarioTree.tree[t].nodes)
                    @everywhere begin
                        n = $n; t = $t; i = $i; ω = $ω; (λ₀, λ₁) = $backwardPassResult[n][1]; 
                        @constraint(forwardInfoList[t-1], forwardInfoList[t-1][:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + 
                                                            sum(λ₁[:s][g] * forwardInfoList[t-1][:s][g] + λ₁[:y][g] * forwardInfoList[t-1][:y][g] for g in indexSets.G) + 
                                                                sum(sum(λ₁[:sur][g][k] * forwardInfoList[t-1][:sur][g, k] for k in keys(solCollection[i, t-1, ω].stageSolution[:sur][g])) for g in indexSets.G));
                        @constraint(backwardInfoList[t-1], backwardInfoList[t-1][:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + 
                                                            sum(λ₁[:s][g] * backwardInfoList[t-1][:s][g] + λ₁[:y][g] * backwardInfoList[t-1][:y][g] for g in indexSets.G) + 
                                                                sum(sum(λ₁[:sur][g][k] * backwardInfoList[t-1][:sur][g, k] for k in keys(solCollection[i, t-1, ω].stageSolution[:sur][g])) for g in indexSets.G));   
                    end
                end

                LM_iter = LM_iter + sum(backwardPassResult[n][2] for n in keys(scenarioTree.tree[t].nodes))
            end
        end
        LM_iter = floor(Int64, LM_iter/sum(length(scenarioTree.tree[t].nodes) for t in 2:indexSets.T));
        
        t1 = now(); iter_time = (t1 - t0).value/1000; total_Time = (t1 - initial).value/1000; i += 1; 

    end
end