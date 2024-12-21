"""
    function stochastic_dual_dynamic_programming_algorithm(
            scenarioTree::ScenarioTree,                   
            indexSets::IndexSets,                        
            paramDemand::ParamDemand,
            paramOPF::ParamOPF;
            initialStateInfo::StateInfo = initialStateInfo,
            param_PLC::NamedTuple = param_PLC, 
            param_levelsetmethod::NamedTuple = param_levelsetmethod, 
            param::NamedTuple = param
    )

# Arguments
  1. `scenarioTree::ScenarioTree`: A scenario tree
  2. `indexSets::IndexSets`: A set of indices
  3. `paramDemand::ParamDemand`: Demand parameters
  4. `paramOPF::ParamOPF`: OPF parameters
  5. `initialStateInfo::StateInfo`: Initial state information
  6. `param_PLC::NamedTuple`: Named tuple of parameters for the PLC algorithm
  7. `param_levelsetmethod::NamedTuple`: Named tuple of parameters for the level set method
  8. `param::NamedTuple`: Named tuple of parameters for the SDDiP algorithm
"""
function stochastic_dual_dynamic_programming_algorithm(
        scenarioTree::ScenarioTree,                   
        indexSets::IndexSets,                        
        paramDemand::ParamDemand,
        paramOPF::ParamOPF;
        initialStateInfo::StateInfo = initialStateInfo,
        param_PLC::NamedTuple = param_PLC, 
        param_levelsetmethod::NamedTuple = param_levelsetmethod, 
        param::NamedTuple = param
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
        ModelList = Dict{Int, SDDPLModel}();
        for t in 1:indexSets.T 
            ModelList[t] = forwardModel!( 
                paramDemand, 
                paramOPF, 
                scenarioTree.tree[t];
                indexSets = indexSets, 
                param = param
            );
        end 
        stateInfoCollection = Dict();  # to store every iteration results
    end

    ####################################################### Main Loop ###########################################################
    while true
        t0 = now(); 
        Ξ̃ = sample_scenarios(; scenarioTree = scenarioTree, numScenarios = param.numScenarios);
        u = Dict{Int64, Float64}();  # to store the value of each scenario
        ####################################################### Forward Steps ###########################################################

        forwardPassResult = pmap(forwardPass, values(Ξ̃));

        for j in 1:param.numScenarios
            ω = collect(keys(Ξ̃))[j]
            for t in 1:indexSets.T
                stateInfoCollection[i, t, ω] = forwardPassResult[ω][t]
            end
            u[ω] = sum(stateInfoCollection[i, t, ω].StageValue for t in 1:indexSets.T);
        end
        @everywhere stateInfoCollection = $stateInfoCollection;
        ####################################################### To Record Info ###########################################################
        LB = maximum([stateInfoCollection[i, 1, 1].StateValue, LB]);
        μ̄ = mean(values(u));
        σ̂² = Statistics.var(values(u));
        UB = μ̄ + 1.96 * sqrt(σ̂²/param.numScenarios); # minimum([μ̄ + 1.96 * sqrt(σ̂²/numScenarios), UB]);
        gap = round((UB-LB)/UB * 100 ,digits = 2); 
        gapString = string(gap,"%"); 
        push!(sddipResult, [i, LB, param.OPT, UB, gapString, iter_time, LM_iter, total_Time, branchDecision]); 
        push!(gapList, gap); 
        branchDecision = false;
        if i == 1
            println("----------------------------------------- Iteration Info ------------------------------------------------")
            println("Iter |        LB        |        UB        |       Gap      |      i-time     |    #LM     |     T-Time")
            println("----------------------------------------------------------------------------------------------------------")
        end
        @printf("%4d | %12.2f     | %12.2f     | %9.2f%%     | %9.2f s     | %6d     | %10.2f s     \n", 
                i, LB, UB, gap, iter_time, LM_iter, total_Time); LM_iter = 0;
        if total_Time > param.TimeLimit || i ≥ param.MaxIter || UB-LB ≤ param.terminate_threshold * UB  
            return Dict(:solHistory => sddipResult, 
                            :solution => stateInfoCollection, 
                                :gapHistory => gapList) 
        end

        ####################################################### Update Partition Tree ###########################################################
        if i ≥ param.LiftingIter                                                                  ## the first rule:: for branching: current convex envelope is good enough
            for t in reverse(1:indexSets.T-1) 
                for ω in [1]#keys(Ξ̃)
                    dev = Dict()
                    for g in indexSets.G 
                        if stateInfoCollection[i, t, ω].BinVar[:y][g] > .5
                            k = maximum([k for (k, v) in stateInfoCollection[i, t, ω].ContVarLeaf[:s][g] if v[:var] > 0.5]);
                            info = ModelList[t].ContVarLeaf[:s][g][k]
                            dev[g] = round(minimum([(info[:ub] - stateInfoCollection[i, t, ω].ContVar[:s][g])/(info[:ub] - info[:lb] + 1e-6), (stateInfoCollection[i, t, ω].ContVar[:s][g] - info[:lb])/(info[:ub] - info[:lb] + 1e-6)]), digits = 5)
                        end
                    end
                    # g = [g for (g, v) in dev if v == maximum(values(dev))][1];
                    large_dev = [g for (g, g_dev) in dev if g_dev ≥ param.branch_threshold]
                    # if dev[g] ≥ param.branch_threshold                                        ## the second rule: current point is far from being an extreme point
                    for g in large_dev
                        branchDecision = true;
                        @everywhere begin
                            t = $t; ω = $ω; g = $g; i = $i; 
                            update_partition_tree!(
                                ModelList, 
                                stateInfoCollection[i, t, ω],
                                t, g; param = param
                            );
                        end
                    end            
                end
            end
        end

        ####################################################### Backward Steps ###########################################################
        for t = reverse(2:indexSets.T)
            for ω in [1]#keys(Ξ̃)
                backwardNodeInfoList = Dict{Int64, Tuple}(); 
                for n in keys(scenarioTree.tree[t].nodes) 
                    backwardNodeInfoList[n] = (i, t, n, ω, param.cutSelection, param_PLC.core_point_strategy) 
                end
                backwardPassResult = pmap(backwardPass, values(backwardNodeInfoList));

                for n in keys(scenarioTree.tree[t].nodes)
                    @everywhere begin
                        n = $n; t = $t; i = $i; ω = $ω; (λ₀, λ₁) = $backwardPassResult[n][1]; 
                        @constraint(
                            ModelList[t-1].model, 
                            ModelList[t-1].model[:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + 
                            sum(λ₁.ContVar[:s][g] * ModelList[t-1].model[:s][g] + λ₁.BinVar[:y][g] * ModelList[t-1].model[:y][g] for g in indexSets.G) + 
                            sum(sum(λ₁.ContAugState[:s][g][k] * ModelList[t-1].model[:augmentVar][g, k] for k in keys(stateInfoCollection[i, t-1, ω].ContAugState[:s][g]); init = 0.0) for g in indexSets.G)
                        );                 
                    end
                end

                LM_iter = LM_iter + sum(backwardPassResult[n][2] for n in keys(scenarioTree.tree[t].nodes))
            end
        end
        LM_iter = floor(Int64, LM_iter/sum(length(scenarioTree.tree[t].nodes) for t in 2:indexSets.T));
        
        t1 = now(); iter_time = (t1 - t0).value/1000; total_Time = (t1 - initial).value/1000; i += 1; 

    end
end