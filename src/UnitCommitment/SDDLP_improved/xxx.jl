for cut in ["LC", "ELC", "SMC"]
    for num in [3, 5, 10]
        for T in [6, 8, 12] 
            indexSets = load(joinpath(project_root, "src", "UnitCommitment", "experiment_$case", "stage($T)real($num)", "indexSets.jld2"))["indexSets"];
            paramOPF = load(joinpath(project_root, "src", "UnitCommitment", "experiment_$case", "stage($T)real($num)", "paramOPF.jld2"))["paramOPF"];
            paramDemand = load(joinpath(project_root, "src", "UnitCommitment", "experiment_$case", "stage($T)real($num)", "paramDemand.jld2"))["paramDemand"];
            scenarioTree = load(joinpath(project_root, "src", "UnitCommitment", "experiment_$case", "stage($T)real($num)", "scenarioTree.jld2"))["scenarioTree"];
            initialStageDecision = load(joinpath(project_root, "src", "UnitCommitment", "experiment_$case", "initialStageDecision.jld2"))["initialStageDecision"];

            # cutSelection = cut;
            @everywhere begin
                indexSets = $indexSets; paramOPF = $paramOPF; paramDemand = $paramDemand; scenarioTree = $scenarioTree; initialStageDecision = $initialStageDecision;
            end
            # Ξ = load("src/UnitCommitment/experiment_$case/stage($T)real($num)/Ξ.jld2")["Ξ"]
            # @time extResult = extensive_form(indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF, scenarioTree = scenarioTree, Ξ = Ξ, initialStageDecision = initialStageDecision); OPT = extResult.OPT;
            if true
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
                                                                                outputFlag = 0, timelimit = para.forwardTimeLimit,
                                                                                        mipGap = para.forwardMipGap, θ_bound = 0.0);
                        var = Dict{Symbol, Dict{Int, VariableRef}}(:s => Dict(g => forwardInfoList[t][:s][g] for g in indexSets.G), :y => Dict(g => forwardInfoList[t][:y][g] for g in indexSets.G)); 
                        sur = Dict{Int, Dict{Int, Dict{Symbol, Any}}}(g => Dict(1 => Dict(:lb => 0., :ub => paramOPF.smax[g], :var =>forwardInfoList[t][:sur][g, 1], :sibling => nothing, :parent => nothing)) for g in indexSets.G); 
                        leaf = Dict{Int, Vector{Int64}}(g => [1] for g in indexSets.G);
                        StateVarList[t] = StateVar(var, sur, leaf);

                        backwardInfoList[t] = backwardModel!(indexSets = indexSets, 
                                                                paramDemand = paramDemand, 
                                                                    paramOPF = paramOPF, 
                                                                        stageRealization = scenarioTree.tree[t], 
                                                                                outputFlag = 0, timelimit = para.backwardTimeLimit,
                                                                                        mipGap = para.backwardMipGap, tightness = para.tightness, θ_bound = 0.0);
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
                    Ξ̃ = sample_scenarios(; scenarioTree = scenarioTree, numScenarios = para.numScenarios);
                    u = Dict{Int64, Float64}();  # to store the value of each scenario
                    ####################################################### Forward Steps ###########################################################

                    forwarPassResult = pmap(forwardPass, values(Ξ̃));

                    for j in 1:para.numScenarios
                        ω = collect(keys(Ξ̃))[j]
                        for t in 1:indexSets.T
                            solCollection[i, t, ω] = forwarPassResult[ω][t]
                        end
                        u[ω] = sum(solCollection[i, t, ω].stageValue for t in 1:indexSets.T);
                    end
                    @everywhere solCollection = $solCollection;
                    ####################################################### To Record Info ###########################################################
                    LB = maximum([solCollection[i, 1, 1].OPT, LB]);
                    μ̄ = mean(values(u));
                    σ̂² = Statistics.var(values(u));
                    UB = μ̄ + 1.96 * sqrt(σ̂²/para.numScenarios); # minimum([μ̄ + 1.96 * sqrt(σ̂²/numScenarios), UB]);
                    gap = round((UB-LB)/UB * 100 ,digits = 2); gapString = string(gap,"%"); push!(sddipResult, [i, LB, para.OPT, UB, gapString, iter_time, LM_iter, total_Time, branchDecision]); push!(gapList, gap); branchDecision = false;
                    if i == 1
                        println("----------------------------------------- Iteration Info ------------------------------------------------")
                        println("Iter |        LB        |        UB        |       Gap      |      i-time     |    #LM     |     T-Time")
                        println("----------------------------------------------------------------------------------------------------------")
                    end
                    @printf("%4d | %12.2f     | %12.2f     | %9.2f%%     | %9.2f s     | %6d     | %10.2f s     \n", 
                            i, LB, UB, gap, iter_time, LM_iter, total_Time); LM_iter = 0;
                    if cutSelection == "ELC" 
                        save("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/MagnantiWong/isddlpHCResult-$ℓ-$tightness.jld2", "sddlpResult", Dict(:solHistory => sddipResult, 
                        :solution => solCollection[i, 1, 1].stageSolution, 
                            :gapHistory => gapList))
                    else
                        save("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/isddlpResult-$cutSelection-$tightness-5.jld2", "sddlpResult", Dict(:solHistory => sddipResult, 
                        :solution => solCollection[i, 1, 1].stageSolution, 
                            :gapHistory => gapList))
                    end
                    if total_Time > para.TimeLimit || i ≥ para.MaxIter || UB-LB ≤ 5e-2 * UB  
                        sddlpResult = Dict(:solHistory => sddipResult, 
                                        :solution => solCollection[i, 1, 1].stageSolution, 
                                            :gapHistory => gapList) 
                        break
                    end

                    ####################################################### Backward Steps ###########################################################
                    if i ≥ indexSets.T * 2                                      ## the first rule:: for branching: current convex envelope is good enough
                        for t in reverse(1:indexSets.T-1) 
                            for ω in [1]#keys(Ξ̃)
                                dev = Dict()
                                for g in indexSets.G 
                                    if solCollection[i, t, ω].stageSolution[:y][g] > .5
                                        k = maximum([k for (k, v) in solCollection[i, t, ω].stageSolution[:sur][g] if v > 0.5]); # maximum(values(solCollection[i, t, ω].stageSolution[:sur][g]))
                                        info = StateVarList[t].sur[g][k]
                                        dev[g] = round(minimum([(info[:ub] - solCollection[i, t, ω].stageSolution[:s][g])/(info[:ub] - info[:lb] + 1e-6), (solCollection[i, t, ω].stageSolution[:s][g] - info[:lb])/(info[:ub] - info[:lb] + 1e-6)]), digits = 5)
                                    end
                                end
                                # g = [g for (g, v) in dev if v == maximum(values(dev))][1];
                                large_dev = [g for (g, g_dev) in dev if g_dev ≥ para.branch_threshold]
                                # if dev[g] ≥ para.branch_threshold                                        ## the second rule: current point is far from being an extreme point
                                for g in large_dev
                                    branchDecision = true;
                                    @everywhere begin
                                        t = $t; ω = $ω; g = $g; i = $i; 
                                        # find the active leaf node 
                                        keys_with_value_1 = maximum([k for (k, v) in solCollection[i, t, ω].stageSolution[:sur][g] if v > 0.5]); ## find the active leaf node: maximum(values(solCollection[i, t, ω].stageSolution[:sur][g]))
                                        # find the lb and ub of this leaf node 
                                        (lb, ub) = StateVarList[t].sur[g][keys_with_value_1][:lb], StateVarList[t].sur[g][keys_with_value_1][:ub]; 
                                        if para.med_method == "interval_mid"
                                            med = (lb + ub)/2; 
                                        elseif para.med_method == "exact_point"
                                            med = solCollection[i, t, ω].stageSolution[:s][g]; # round(solCollection[i, t, ω].stageSolution[:s][g], digits = 3); 
                                        end
                                        # create two new leaf nodes, and update their info (lb, ub)
                                        left = length(StateVarList[t].sur[g]) + 1; right = length(StateVarList[t].sur[g]) + 2; 
                                        forwardInfoList[t][:sur][g, left] = @variable(forwardInfoList[t], base_name = "sur[$g, $left]", binary = true); 
                                        forwardInfoList[t][:sur][g, right] = @variable(forwardInfoList[t], base_name = "sur[$g, $right]", binary = true);
                                        StateVarList[t].sur[g][left] = Dict(:lb => lb, :ub => med, :var => forwardInfoList[t][:sur][g, left], :sibling => right, :parent => keys_with_value_1)
                                        StateVarList[t].sur[g][right] =  Dict(:lb => med, :ub => ub, :var => forwardInfoList[t][:sur][g, right], :sibling => left, :parent => keys_with_value_1)
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
                                        if solCollection[i, t, ω].stageSolution[:s][g] ≤ med
                                            stageDecision[:sur][g][left] = 1.0
                                            # stageDecision[:sur][g][right] = 0.0
                                        else
                                            stageDecision[:sur][g][right] = 1.0
                                            # stageDecision[:sur][g][left] = 0.0
                                        end
                                        solCollection[i, t, ω] = ( stageSolution = deepcopy(stageDecision), 
                                                                    stageValue = solCollection[i, t, ω].stageValue, 
                                                                    OPT = solCollection[i, t, ω].OPT
                                                                );
                                    end
                                end            
                            end
                        end
                    end

                    ####################################################### Backward Steps ###########################################################
                    for t = reverse(2:indexSets.T)
                        for ω in [1]#keys(Ξ̃)
                            backwardNodeInfoSet = Dict{Int64, Tuple}(); 
                            for n in keys(scenarioTree.tree[t].nodes) backwardNodeInfoSet[n] = (i, t, n, ω, para.cutSelection, para.core_point_strategy) end
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
            if cutSelection == "ELC" 
                save("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/MagnantiWong/isddlpHCResult-$ℓ-$tightness.jld2", "sddlpResult", sddlpResult)
            else
                save("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/isddlpResult-$cutSelection-$tightness.jld2", "sddlpResult", sddlpResult)
            end
            
            # remove all redundant variables for processes
            @everywhere begin
                forwardInfoList = nothing;
                backwardInfoList = nothing;
                StateVarList = nothing; 
                stageDecision = nothing; 
                λ₀ = nothing; λ₁ = nothing;
                solCollection = nothing;  # to store every iteration results
            end
            @everywhere GC.gc(); # garbage collection
        end
    end
end