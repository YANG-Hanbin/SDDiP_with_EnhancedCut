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
                initial = now(); i = 1; LB = - Inf; UB = Inf; 
                iter_time = 0; total_Time = 0; t0 = 0.0; LMiter = 0; LM_iter = 0; gap = 100.0; gapString = "100%";
            
                col_names = [:Iter, :LB, :OPT, :UB, :gap, :time, :LM_iter, :Time]; # needs to be a vector Symbols
                col_types = [Int64, Float64, Union{Float64,Nothing}, Float64, String, Float64, Int64, Float64];
                named_tuple = (; zip(col_names, type[] for type in col_types )...);
                sddipResult = DataFrame(named_tuple); # 0×7 DataFrame
                gapList = [];
               
                @everywhere begin
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
                                                                                        mipGap = 1e-2, θ_bound = 0.0, timelimit = 20
                                                                );
                        backwardInfoList[t] = backwardModel!(indexSets = indexSets, 
                                                                paramDemand = paramDemand, 
                                                                    paramOPF = paramOPF, 
                                                                        stageRealization = scenarioTree.tree[t], 
                                                                                outputFlag = 0, κ = κ, ε = ε, 
                                                                                        mipGap = 1e-2, tightness = tightness, θ_bound = 0.0, timelimit = 20
                                                                );
                    end 
                    stageDecision = Dict{Symbol, Dict{Int64, Float64}}(); stageDecision[:s] = Dict{Int64, Float64}(); stageDecision[:y] = Dict{Int64, Float64}(); stageDecision[:λ] = Dict{Int64, Dict{Int64, Float64}}();
                    λ₀ = 0.0; λ₁ = Dict{Symbol, Dict{Int64, Float64}}();
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
                    LB = maximum([solCollection[i, 1, 1].OPT, LB]);
                    μ̄ = mean(values(u));
                    σ̂² = Statistics.var(values(u));
                    UB = μ̄ + 1.96 * sqrt(σ̂²/numScenarios); # minimum([μ̄ + 1.96 * sqrt(σ̂²/numScenarios), UB]);
                    gap = round((UB-LB)/UB * 100 ,digits = 2); gapString = string(gap,"%"); push!(sddipResult, [i, LB, OPT, UB, gapString, iter_time, LM_iter, total_Time]); push!(gapList, gap); 
                    if i == 1
                        println("----------------------------------------- Iteration Info ------------------------------------------------")
                        println("Iter |        LB        |        UB        |       Gap      |      i-time     |    #LM     |     T-Time")
                        println("----------------------------------------------------------------------------------------------------------")
                    end
                    @printf("%4d | %12.2f     | %12.2f     | %9.2f%%     | %9.2f s     | %6d     | %10.2f s     \n", 
                            i, LB, UB, gap, iter_time, LM_iter, total_Time); LM_iter = 0;
                    save("src/UnitCommitment/numericalResults-$case/Periods$T-Real$num/sddipResult-$cutSelection-$tightness-$coef.jld2", "sddipResult", Dict(:solHistory => sddipResult, 
                    :solution => solCollection[i, 1, 1].stageSolution, 
                        :gapHistory => gapList))
                    if total_Time > TimeLimit || i > MaxIter # || UB-LB ≤ 1e-2 * UB
                        sddipResult = Dict(:solHistory => sddipResult, 
                                        :solution => solCollection[i, 1, 1].stageSolution, 
                                            :gapHistory => gapList) 
                        break;
                    end
            
                    ####################################################### Backward Steps ###########################################################
                    for t = reverse(2:indexSets.T)
                        for ω in [1]#keys(Ξ̃)
                            backwardNodeInfoSet = Dict{Int64, Tuple}();
                            for n in keys(scenarioTree.tree[t].nodes) backwardNodeInfoSet[n] = (i, t, n, ω, cutSelection) end
                            backwardPassResult = pmap(backwardPass, values(backwardNodeInfoSet));
            
                            for n in keys(scenarioTree.tree[t].nodes)
                                @everywhere begin
                                    n = $n; t = $t; i = $i; ω = $ω; (λ₀, λ₁) = $backwardPassResult[n][1]; 
                                    @constraint(forwardInfoList[t-1], forwardInfoList[t-1][:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + sum(sum(λ₁[:λ][g][k] * forwardInfoList[t-1][:λ][g, k] for k in 1:κ[g]) + λ₁[:y][g] * forwardInfoList[t-1][:y][g] for g in indexSets.G));
                                    @constraint(backwardInfoList[t-1], backwardInfoList[t-1][:θ][n]/scenarioTree.tree[t-1].prob[n] ≥ λ₀ + sum(sum(λ₁[:λ][g][k] * backwardInfoList[t-1][:λ][g, k] for k in 1:κ[g]) + λ₁[:y][g] * backwardInfoList[t-1][:y][g] for g in indexSets.G));
                                end
                            end
            
                            LM_iter = LM_iter + sum(backwardPassResult[n][2] for n in keys(scenarioTree.tree[t].nodes))
                        end
                    end
                    LM_iter = floor(Int64, LM_iter/sum(length(scenarioTree.tree[t].nodes) for t in 2:indexSets.T));
                    
                    t1 = now(); iter_time = (t1 - t0).value/1000; total_Time = (t1 - initial).value/1000; i += 1;
                end

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