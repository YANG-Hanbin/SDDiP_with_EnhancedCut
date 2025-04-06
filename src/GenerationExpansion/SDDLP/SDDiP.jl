function SDDiP_algorithm(
    Ω::Dict{Int64,Dict{Int64,RandomVariables}}, 
    probList::Dict{Int64,Vector{Float64}}, 
    stageDataList::Dict{Int64, StageData}; 
    binaryInfo::BinaryInfo = binaryInfo,
    param::NamedTuple = param
)
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
    @everywhere begin
        forwardInfoList = Dict{Int, ForwardModelInfo}();
        backwardInfoList = Dict{Int, BackwardModelInfo}();
        StateVarList = Dict(); 
        for t in 1:param.T
            forwardInfoList[t] = forwardModel!(stageDataList[t], binaryInfo = binaryInfo, timelimit = 10, mipGap = 1e-4);
            var = Dict{Symbol, Dict{Int, VariableRef}}(:St => Dict(g => forwardInfoList[t].model[:St][g] for g in 1:binaryInfo.d)); 
            sur = Dict{Int, Dict{Int, Dict{Symbol, Any}}}(g => Dict(1 => Dict(:lb => 0., :ub => stageDataList[t].ū[g], :var =>forwardInfoList[t].model[:sur][g, 1])) for g in 1:binaryInfo.d); 
            leaf = Dict{Int, Vector{Int64}}(g => [1] for g in 1:binaryInfo.d);
            StateVarList[t] = StateVar(var, sur, leaf);

            backwardInfoList[t] = backwardModel!(stageDataList[t], binaryInfo = binaryInfo, timelimit = 10, mipGap = 1e-4, tightness = param.tightness);
        end 
    end
    initial = now(); iter_time = 0.; total_Time = 0.; t0 = 0.0;

    while true
        t0 = now();
        solCollection = Dict();  # to store every iteration results
        u = Vector{Float64}(undef, param.M);  # to compute upper bound 
        Random.seed!(i)
        Scenarios = SampleScenarios(
            Ω, 
            probList, 
            M = param.M
        );
        
        ## Forward pass
        forwardPassResult = pmap(1:param.M) do k
            forwardPass(
                k, 
                Scenarios;
            )
        end;

        for k in 1:param.M
            for t in 1:param.T
                solCollection[t, k] = forwardPassResult[k][t, k]
            end
            u[k] = sum(solCollection[t, k].stageValue for t in 1:param.T);
        end
        
        ## compute the upper bound
        LB = solCollection[1, 1].OPT;
        μ̄ = mean(u);
        σ̂² = var(u);
        UB = μ̄ + 1.96 * sqrt(σ̂²/param.M); # minimum([μ̄ + 1.96 * sqrt(σ̂²/M), UB]);
        gap = round((UB-LB)/UB * 100 ,digits = 2);
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

        ####################################################### Adding Binary Variables ###########################################################
        for t in 1:param.T-1 
            for ω in [1]#keys(Ξ̃)
                dev = Dict()
                for g in 1:binaryInfo.d
                    k = maximum([k for (k, v) in solCollection[t, ω].stageSur[g] if v == maximum(values(solCollection[t, ω].stageSur[g]))])
                    info = StateVarList[t].sur[g][k]
                    dev[g] = minimum([(info[:ub] - solCollection[t, ω].stageSolution[g])/(info[:ub] - info[:lb] + 1e-6), (solCollection[t, ω].stageSolution[g] - info[:lb])/(info[:ub] - info[:lb] + 1e-6)])
                end
                g = [k for (k, v) in dev if v == maximum(values(dev))][1]
                if dev[g] ≥ 1e-8
                    # find the active leaf node 
                    keys_with_value_1 = maximum([k for (k, v) in solCollection[t, ω].stageSur[g] if v == 1])
                    # find the lb and ub of this leaf node 
                    (lb, ub) = StateVarList[t].sur[g][keys_with_value_1][:lb], StateVarList[t].sur[g][keys_with_value_1][:ub]; med = solCollection[t, ω].stageSolution[g]; # solCollection[i, t, ω].stageSolution[:s][g];# (lb + ub)/2; #round(solCollection[i, t, ω].stageSolution[:s][g], digits = 3); 
                    # create two new leaf nodes, and update their info (lb, ub)
                    left = length(StateVarList[t].sur[g]) + 1; right = length(StateVarList[t].sur[g]) + 2;
                    @everywhere begin
                        t = $t; left = $left; right = $right; g = $g; 
                        forwardInfoList[t].model[:sur][g, left] = @variable(forwardInfoList[t].model, base_name = "sur[$g, $left]", binary = true); 
                        forwardInfoList[t].model[:sur][g, right] = @variable(forwardInfoList[t].model, base_name = "sur[$g, $right]", binary = true);
                    end
                    StateVarList[t].sur[g][left] = Dict(:lb => lb, :ub => med, :var => forwardInfoList[t].model[:sur][g, left])
                    StateVarList[t].sur[g][right] =  Dict(:lb => med, :ub => ub, :var => forwardInfoList[t].model[:sur][g, right])
                    # pop and push new leaf nodes
                    deleteat!(StateVarList[t].leaf[g], findall(x -> x == keys_with_value_1, StateVarList[t].leaf[g])); push!(StateVarList[t].leaf[g], left); push!(StateVarList[t].leaf[g], right);

                    # add logic constraints
                    ## for forward models
                    ### Parent-Child relationship
                    @everywhere begin
                        t = $t; left = $left; right = $right; g = $g; keys_with_value_1 = $keys_with_value_1; StateVarList = $StateVarList;
                        @constraint(
                            forwardInfoList[t].model, 
                            forwardInfoList[t].model[:sur][g, left] + forwardInfoList[t].model[:sur][g, right] == forwardInfoList[t].model[:sur][g, keys_with_value_1]
                        );
                        ### bounding constraints
                        @constraint(
                            forwardInfoList[t].model, 
                            forwardInfoList[t].St[g] ≥ sum(StateVarList[t].sur[g][k][:lb] * forwardInfoList[t].model[:sur][g, k] for k in StateVarList[t].leaf[g])
                        );
                        @constraint(
                            forwardInfoList[t].model, 
                            forwardInfoList[t].St[g] ≤ sum(StateVarList[t].sur[g][k][:ub] * forwardInfoList[t].model[:sur][g, k] for k in StateVarList[t].leaf[g])
                        );
                    
                        ## for backward models
                        ### Parent-Child relationship
                        if param.tightness
                            backwardInfoList[t+1].model[:sur_copy][g, left] = @variable(backwardInfoList[t+1].model, base_name = "sur_copy[$g, $left]", binary = true); 
                            backwardInfoList[t+1].model[:sur_copy][g, right] = @variable(backwardInfoList[t+1].model, base_name = "sur_copy[$g, $left]", binary = true); 
                            backwardInfoList[t].model[:sur][g, left] = @variable(backwardInfoList[t].model, base_name = "sur[$g, $left]", binary = true); 
                            backwardInfoList[t].model[:sur][g, right] = @variable(backwardInfoList[t].model, base_name = "sur[$g, $left]", binary = true); 
                        else
                            backwardInfoList[t+1].model[:sur_copy][g, left] = @variable(backwardInfoList[t+1].model, base_name = "sur_copy[$g, $left]", lower_bound = 0, upper_bound = 1); 
                            backwardInfoList[t+1].model[:sur_copy][g, right] = @variable(backwardInfoList[t+1].model, base_name = "sur_copy[$g, $left]", lower_bound = 0, upper_bound = 1); 
                            backwardInfoList[t].model[:sur][g, left] = @variable(backwardInfoList[t].model, base_name = "sur[$g, $left]", lower_bound = 0, upper_bound = 1); 
                            backwardInfoList[t].model[:sur][g, right] = @variable(backwardInfoList[t].model, base_name = "sur[$g, $left]", lower_bound = 0, upper_bound = 1); 
                        end
                        @constraint(backwardInfoList[t+1].model, backwardInfoList[t+1].model[:sur_copy][g, left] + backwardInfoList[t+1].model[:sur_copy][g, right] == backwardInfoList[t+1].model[:sur_copy][g, keys_with_value_1]);
                        @constraint(backwardInfoList[t].model, backwardInfoList[t].model[:sur][g, left] + backwardInfoList[t].model[:sur][g, right] == backwardInfoList[t].model[:sur][g, keys_with_value_1]);
                        ### bounding constraints
                        @constraint(backwardInfoList[t+1].model, backwardInfoList[t+1].model[:Sc][g] ≥ sum(StateVarList[t].sur[g][k][:lb] * backwardInfoList[t+1].model[:sur_copy][g, k] for k in StateVarList[t].leaf[g]));
                        @constraint(backwardInfoList[t+1].model, backwardInfoList[t+1].model[:Sc][g] ≤ sum(StateVarList[t].sur[g][k][:ub] * backwardInfoList[t+1].model[:sur_copy][g, k] for k in StateVarList[t].leaf[g]));
                        @constraint(backwardInfoList[t].model, backwardInfoList[t].St[g] ≥ sum(StateVarList[t].sur[g][k][:lb] * backwardInfoList[t].model[:sur][g, k] for k in StateVarList[t].leaf[g]));
                        @constraint(backwardInfoList[t].model, backwardInfoList[t].St[g] ≤ sum(StateVarList[t].sur[g][k][:ub] * backwardInfoList[t].model[:sur][g, k] for k in StateVarList[t].leaf[g]));
                    end
                    stageDecision = deepcopy(solCollection[t, ω]);
                    stageDecision.stageSur[g] = Dict{Int64, Float64}()
                    for k in StateVarList[t].leaf[g]
                        stageDecision.stageSur[g][k] = 0.0
                    end
                    if solCollection[t, ω].stageSolution[g] ≤ med
                        stageDecision.stageSur[g][left] = 1.0
                        stageDecision.stageSur[g][right] = 0.0
                    else
                        stageDecision.stageSur[g][right] = 1.0
                        stageDecision.stageSur[g][left] = 0.0
                    end
                    solCollection[t, ω] = deepcopy(stageDecision);
                end
            end
        end

        ####################################################### Cut Generation Processes ###########################################################
        for t = reverse(2:param.T)
            for k in [1] 
                backwardNodeInfoList = Dict{Int64, Tuple}(); 
                for j in keys(Ω[t]) 
                    backwardNodeInfoList[j] = (t, j, k) 
                end

                backwardPassResult = pmap(values(backwardNodeInfoList)) do backwardNodeInfo
                    backwardPass(
                        backwardNodeInfo, 
                        solCollection; 
                    )
                end;

                L̂ = solCollection[t-1,k].stageSolution; 
                sur = solCollection[t-1,k].stageSur;
                (λ₀, λ₁) = (
                    0.0, 
                    Dict( 
                        :St => L̂ .* 0.0,                 
                        :sur => Dict(
                            g => Dict(
                                k => 0.0 for k in keys(sur[g])
                                ) for g in 1:binaryInfo.d
                        )
                    )
                );
                for j in keys(Ω[t])
                    λ₀ = λ₀ + probList[t][j] * backwardPassResult[j][1];
                    λ₁ = Dict(
                        :St => λ₁[:St] .+ probList[t][j] * backwardPassResult[j][2][:St],
                        :sur => Dict(
                            g => Dict(
                                k => λ₁[:sur][g][k] + probList[t][j] * backwardPassResult[j][2][:sur][g][k] for k in keys(sur[g]) 
                            ) for g in 1:binaryInfo.d 
                        )
                    );
                end
                # # add cut to both backward models and forward models
                @everywhere begin
                    t = $t; λ₁ = $λ₁; λ₀ = $λ₀;
                    @constraint(
                        forwardInfoList[t-1].model, 
                        forwardInfoList[t-1].θ ≥ λ₀ + 
                        λ₁[:St]' * forwardInfoList[t-1].St + 
                        sum(
                            sum(
                                λ₁[:sur][g][k] * forwardInfoList[t-1].model[:sur][g, k] for k in StateVarList[t-1].leaf[g]
                            ) for g in 1:binaryInfo.d
                        ) 
                    );
                    @constraint(
                        backwardInfoList[t-1].model, 
                        backwardInfoList[t-1].θ ≥ λ₀ + 
                        λ₁[:St]' * backwardInfoList[t-1].St + 
                        sum(
                            sum(
                                λ₁[:sur][g][k] * backwardInfoList[t-1].model[:sur][g, k] for k in StateVarList[t-1].leaf[g]
                            ) for g in 1:binaryInfo.d
                        ) 
                    );  
                end         
            end
        end

        t1 = now();
        iter_time = (t1 - t0).value/1000;
        total_Time = (t1 - initial).value/1000;
        i = i + 1;
    end

end