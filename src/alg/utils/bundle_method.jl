"""
    function SetupLevelSetMethodOracleParam(
        param::NamedTuple
    )::LevelSetMethodOracleParam

# Arguments

    1. `param::NamedTuple` : the parameters of the level set method

# Returns

    1. `LevelSetMethodOracleParam` : the parameters of the level set method

"""
function SetupLevelSetMethodOracleParam(
    stateInfo::StateInfo;
    indexSets::IndexSets = indexSets,
    param_levelsetmethod::NamedTuple = param_levelsetmethod
)::LevelSetMethodOracleParam
    
    μ             = get(param_levelsetmethod, :μ, 0.9);
    λ             = get(param_levelsetmethod, :λ, 0.5);
    threshold     = get(param_levelsetmethod, :threshold, nothing);
    nxt_bound     = get(param_levelsetmethod, :nxt_bound, 1e10);
    MaxIter       = get(param_levelsetmethod, :MaxIter, 200);
    verbose       = get(param_levelsetmethod, :levelsetmethod_verbose, true);

    return LevelSetMethodOracleParam(
        μ, 
        λ, 
        threshold, 
        nxt_bound, 
        MaxIter, 
        verbose, 
        setup_initial_point(
            stateInfo;
            indexSets = indexSets 
        )
    )
end

"""
   function Δ_model_formulation(
        functionHistory::FunctionHistory, 
        f_star::Float64, 
        iter::Int64; 
        para::NamedTuple = para
    )

# Arguments

    1. `functionHistory::FunctionHistory` : the history of the function
    2. `f_star::Float64` : the optimal value of the current approximate value function
    3. `iter::Int64` : the iteration number
    4. `para::NamedTuple` : the parameters of the level set method

# Returns

    1. `Dict{Int64, Float64}` : the value of the gap and the bounds of alpha
"""
function Δ_model_formulation(
    functionHistory::FunctionHistory, 
    f_star::Float64, 
    iter::Int64; 
    param::NamedTuple = param
)
    
    alphaModel = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                                    "Threads" => 0)); 
    MOI.set(alphaModel, MOI.Silent(), !param.verbose);
    set_optimizer_attribute(alphaModel, "MIPGap", param.MIPGap);
    set_optimizer_attribute(alphaModel, "TimeLimit", param.TimeLimit);

    @variable(alphaModel, z);
    @variable(alphaModel, 0 ≤ α ≤ 1);
    @constraint(alphaModel, con[j = 1:iter], z ≤  α * ( functionHistory.f_his[j] - f_star) + (1 - α) * functionHistory.G_max_his[j] );
    
    # we first compute gap Δ
    @objective(alphaModel, Max, z);
    optimize!(alphaModel);
    st = termination_status(alphaModel);
    Δ = JuMP.objective_value(alphaModel);

    
    ## then we modify above model to compute alpha
    # alpha_min
    @constraint(alphaModel, z .≥ 0);
    @objective(alphaModel, Min, α);
    optimize!(alphaModel);
    a_min = JuMP.value(α);

    # alpha_max
    @objective(alphaModel, Max, α);
    optimize!(alphaModel);
    a_max = JuMP.value(α);

    return Dict(1 => Δ, 2 => a_min, 3 => a_max)

end


"""
    function add_constraint(
        currentInfo::CurrentInfo,
        modelInfo::ModelInfo;
        indexSets::IndexSets = indexSets
    )::Nothing

# Arguments

    1. `currentInfo::CurrentInfo` : the current information
    2. `modelInfo::ModelInfo` : the model information
    3. `indexSets::IndexSets` : the index sets
    
# Returns
    This function is to add constraints for the model f_star and nxt pt.
"""
function add_constraint(
    currentInfo::CurrentInfo,
    modelInfo::ModelInfo;
    indexSets::IndexSets = indexSets
)::Nothing
    m = length(currentInfo.G)

    # add constraints     
    @constraint(
        modelInfo.model, 
        modelInfo.z .≥ currentInfo.f + 
        sum(currentInfo.df[:s][g] * (modelInfo.xs[g] - currentInfo.x.ContVar[:s][g]) for g in indexSets.G) + 
        sum(currentInfo.df[:y][g] * (modelInfo.xy[g] - currentInfo.x.BinVar[:y][g]) for g in indexSets.G) + 
        sum(sum(currentInfo.df[:sur][g][k] * (modelInfo.sur[g, k] - currentInfo.x.ContAugState[:s][g][k]) for k in keys(currentInfo.df[:sur][g]); init = 0.0) for g in indexSets.G)
    );

    @constraint(
        modelInfo.model, 
        [k = 1:m], 
        modelInfo.y .≥ currentInfo.G[k] + 
        sum(currentInfo.dG[k].ContVar[:s][g] * (modelInfo.xs[g] .- currentInfo.x.ContVar[:s][g]) for g in keys(currentInfo.df[:s])) + 
        sum(currentInfo.dG[k].BinVar[:y][g] * (modelInfo.xy[g] .- currentInfo.x.BinVar[:y][g]) for g in keys(currentInfo.df[:y])) + 
        sum(sum(currentInfo.dG[k].ContAugState[:s][g][j] * (modelInfo.sur[g, j] .- currentInfo.x.ContAugState[:s][g][j]) for j in keys(currentInfo.df[:sur][g]); init = 0.0) for g in indexSets.G)
    );

    return                                                                              
end


"""
    function solve_inner_minimization_problem(
        PLCGeneration::ParetoLagrangianCutGeneration,
        model::Model, 
        x₀::StateInfo, 
        stateInfo::StateInfo
    )

# Arguments

    1. `PLCGeneration::ParetoLagrangianCutGeneration` : the information of the cut that will be generated information
    2. `model::Model` : the backward model
    3. `x₀::StateInfo` : the dual information
    4. `stateInfo::StateInfo` : the last stage decision
  
# Returns
    1. `currentInfo::CurrentInfo` : the current information
  
"""
function solve_inner_minimization_problem(
    PLCGeneration::ParetoLagrangianCutGeneration,
    model::Model, 
    x₀::StateInfo, 
    stateInfo::StateInfo
)
    @objective(
        model, 
        Min,  
        model[:primal_objective_expression] +
        sum(
            x₀.ContVar[:s][g] * (stateInfo.ContVar[:s][g] - model[:s_copy][g]) + 
            x₀.BinVar[:y][g] * (stateInfo.BinVar[:y][g] - model[:y_copy][g]) + 
            sum(x₀.ContAugState[:s][g][k] * (stateInfo.ContAugState[:s][g][k] - model[:augmentVar_copy][g, k]) for k in keys(stateInfo.ContAugState[:s][g])) 
            for g in indexSets.G
        )
    );
    ## ==================================================== solve the model and display the result ==================================================== ##
    optimize!(model);
    F  = JuMP.objective_value(model);
    
    negative_∇F = StateInfo(
        Dict{Any, Dict{Any, Any}}(:y => Dict{Any, Any}(
            g => JuMP.value(model[:y_copy][g]) - stateInfo.BinVar[:y][g] for g in indexSets.G)
        ), 
        nothing, 
        Dict{Any, Dict{Any, Any}}(:s => Dict{Any, Any}(
            g => JuMP.value(model[:s_copy][g]) - stateInfo.ContVar[:s][g] for g in indexSets.G)
        ), 
        nothing, 
        nothing, 
        nothing, 
        nothing, 
        nothing, 
        Dict{Any, Dict{Any, Dict{Any, Any}}}(
        :s => Dict{Any, Dict{Any, Any}}(
            g => Dict{Any, Any}(
                k => JuMP.value(model[:augmentVar_copy][g, k]) - stateInfo.ContAugState[:s][g][k] for k in keys(stateInfo.ContAugState[:s][g])
            ) for g in indexSets.G)
        )
    );
    currentInfo = CurrentInfo(  
        x₀, 
        - F - sum(x₀.ContVar[:s][g] * (PLCGeneration.core_point.ContVar[:s][g] .- stateInfo.ContVar[:s][g]) + x₀.BinVar[:y][g] * (PLCGeneration.core_point.BinVar[:y][g] - stateInfo.BinVar[:y][g]) + 
            sum(x₀.ContAugState[:s][g][k] * (PLCGeneration.core_point.ContAugState[:s][g][k] - stateInfo.ContAugState[:s][g][k]) for k in keys(stateInfo.ContAugState[:s][g]); init = 0.0) 
                for g in indexSets.G
        ),                                                                                                                                                              ## obj function value
        Dict(
            1 => PLCGeneration.primal_bound - F - PLCGeneration.δ
        ),                                                                                                                                                              ## constraint value
        Dict{Symbol, Dict{Int64, Any}}(
            :s => Dict(g => negative_∇F.ContVar[:s][g] + stateInfo.ContVar[:s][g] - PLCGeneration.core_point.ContVar[:s][g] for g in indexSets.G),
            :y => Dict(g => negative_∇F.BinVar[:y][g] + stateInfo.BinVar[:y][g] - PLCGeneration.core_point.BinVar[:y][g] for g in indexSets.G), 
            :sur => Dict(g => Dict(k => negative_∇F.ContAugState[:s][g][k] + stateInfo.ContAugState[:s][g][k] - PLCGeneration.core_point.ContAugState[:s][g][k] for k in keys(stateInfo.ContAugState[:s][g])) for g in indexSets.G)
        ),                                                                                                                                                              ## obj gradient
        Dict(1 => negative_∇F )                                                                                                                                         ## constraint gradient
    );
    return (currentInfo = currentInfo, currentInfo_f = F)
end    


"""
LevelSetMethod_optimization!(; stageDecision::Dict{Symbol, Dict{Int64, Float64}} = stageDecision,
                                    f_star_value::Float64 = f_star_value,
                                        x_interior::Union{Dict{Symbol, Dict{Int64, Float64}}, Nothing} = nothing,
                                            x₀::Dict{Symbol, Dict{Int64, Float64}} = x₀
                            )

# Arguments

    1. `stageDecision::Dict{Symbol, Dict{Int64, Float64}}` : the decision of the last stage
    2. `f_star_value::Float64` : the optimal value of the current approximate value function
    3. `x_interior::Union{Dict{Symbol, Dict{Int64, Float64}}, Nothing}` : an interior point
    4. `x₀::Dict{Symbol, Dict{Int64, Float64}}` : the initial point of the lagrangian dual variables
    5. `model::Model` : backward model
  
# Returns
    1. `cutInfo::Array{Any,1}` : the cut information
  
"""
function LevelSetMethod_optimization!(
    model::Model, 
    levelsetmethodOracleParam::LevelSetMethodOracleParam, 
    stateInfo::StateInfo,
    PLCGeneration::ParetoLagrangianCutGeneration;
    indexSets::IndexSets = indexSets, 
    paramDemand::ParamDemand = paramDemand, 
    paramOPF::ParamOPF = paramOPF, 
    param::NamedTuple = param, param_levelsetmethod::NamedTuple = param_levelsetmethod
)

    ## ==================================================== auxiliary function for function information ==================================================== ##
    # μ larger is better
    # (μ, λ, threshold, nxt_bound, max_iter, verbose, x₀) = (
    #     levelsetmethodOracleParam.μ, 
    #     levelsetmethodOracleParam.λ, 
    #     levelsetmethodOracleParam.threshold, 
    #     levelsetmethodOracleParam.nxt_bound, 
    #     levelsetmethodOracleParam.MaxIter, 
    #     levelsetmethodOracleParam.verbose, 
    #     levelsetmethodOracleParam.x₀
    # );
    (core_point, primal_bound) = (PLCGeneration.core_point, PLCGeneration.primal_bound);
    (D, G, L, B) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B);
    
    ## ==================================================== Levelset Method ============================================== ##    
    iter = 1;
    α = 1/2;
    Δ = Inf; 

    # trajectory
    currentInfo, currentInfo_f = solve_inner_minimization_problem(
        PLCGeneration,
        model, 
        levelsetmethodOracleParam.x₀, 
        stateInfo
    );

    functionHistory = FunctionHistory(  
        Dict(1 => currentInfo.f), 
        Dict(1 => maximum(currentInfo.G[k] for k in keys(currentInfo.G)) )
    );

    # model for oracle
    oracleModel = Model(optimizer_with_attributes(
        ()->Gurobi.Optimizer(GRB_ENV), 
        "Threads" => 0)
    ); 
    MOI.set(oracleModel, MOI.Silent(), true);
    set_optimizer_attribute(oracleModel, "MIPGap", param.MIPGap);
    set_optimizer_attribute(oracleModel, "TimeLimit", param.TimeLimit);

    ## ==================================================== Levelset Method ============================================== ##
    para_oracle_bound = abs(currentInfo.f);
    z_rhs = 5 * 10^(ceil(log10(para_oracle_bound)));
    @variable(oracleModel, z ≥ - z_rhs);
    @variable(oracleModel, xs_oracle[G]);
    @variable(oracleModel, xy_oracle[G]);
    @variable(oracleModel, sur_oracle[g in G, k in keys(stateInfo.ContAugState[:s][g])]);
    @variable(oracleModel, y ≤ 0);

    @objective(oracleModel, Min, z);
    oracleInfo = ModelInfo(oracleModel, xs_oracle, xy_oracle, sur_oracle, y, z);


    nxtModel = Model(optimizer_with_attributes(
        ()->Gurobi.Optimizer(GRB_ENV), 
        "Threads" => 0)
    ); 
    MOI.set(nxtModel, MOI.Silent(), true);
    set_optimizer_attribute(nxtModel, "MIPGap", param.MIPGap);
    set_optimizer_attribute(nxtModel, "TimeLimit", param.TimeLimit);

    @variable(nxtModel, xs[G]);
    @variable(nxtModel, xy[G]);
    @variable(nxtModel, sur[g in G, k in keys(stateInfo.ContAugState[:s][g])]);
    @variable(nxtModel, z1);
    @variable(nxtModel, y1);
    nxtInfo = ModelInfo(nxtModel, xs, xy, sur, y1, z1);

    cutInfo =  [ 
        - currentInfo.f - 
        sum(currentInfo.x.ContVar[:s][g] * core_point.ContVar[:s][g] + 
        currentInfo.x.BinVar[:y][g] * core_point.BinVar[:y][g] for g in G) - 
        sum(
            sum(
                currentInfo.x.ContAugState[:s][g][k] * core_point.ContAugState[:s][g][k] for k in keys(stateInfo.ContAugState[:s][g]); init = 0.0
                ) for g in G
        ), 
        currentInfo.x
    ]; 

    while true
        add_constraint(currentInfo, oracleInfo);
        optimize!(oracleModel);
        st = termination_status(oracleModel);
        if termination_status(oracleModel) == MOI.OPTIMAL
            f_star = JuMP.objective_value(oracleModel);
        else 
            # @info "Oracle Model is $(st)!"
            return (cutInfo = cutInfo, iter = iter)
        end

        # formulate alpha model
        result = Δ_model_formulation(functionHistory, f_star, iter);
        previousΔ = copy.(Δ);
        Δ, a_min, a_max = result[1], result[2], result[3];

        if param_levelsetmethod.verbose # && (iter % 30 == 0)
            if iter == 1
                println("------------------------------------ Iteration Info --------------------------------------")
                println("Iter |   Gap                              Objective                             Constraint")
            end
            @printf("%3d  |   %5.3g                         %5.3g                              %5.3g\n", iter, Δ, - currentInfo.f, currentInfo.G[1])
        end

        x₀ = currentInfo.x;
        if (round(previousΔ) > round(Δ)) || ((currentInfo.G[1] ≤ 0.0))
            cutInfo =  [ 
                - currentInfo.f - 
                sum(currentInfo.x.ContVar[:s][g] * core_point.ContVar[:s][g] + 
                currentInfo.x.BinVar[:y][g] * core_point.BinVar[:y][g] for g in G) - 
                sum(
                    sum(
                        currentInfo.x.ContAugState[:s][g][k] * core_point.ContAugState[:s][g][k] for k in keys(stateInfo.ContAugState[:s][g]); init = 0.0
                        ) for g in G
                ), 
                currentInfo.x
            ]; 
        end

        # update α
        if param_levelsetmethod.μ/2 ≤ (α-a_min)/(a_max-a_min) .≤ 1 - param_levelsetmethod.μ/2
            α = α;
        else
            α = round.((a_min+a_max)/2, digits = 6);
        end

        # update level
        w = α * f_star;
        W = minimum( α * functionHistory.f_his[j] + (1-α) * functionHistory.G_max_his[j] for j in 1:iter);
        
        level = w + param_levelsetmethod.λ * (W - w);
        
        ## ==================================================== next iteration point ============================================== ##
        # obtain the next iteration point
        if iter == 1
            @constraint(nxtModel, levelConstraint, α * z1 + (1 - α) * y1 ≤ level);
        else 
            delete(nxtModel, nxtModel[:levelConstraint]);
            unregister(nxtModel, :levelConstraint);
            @constraint(nxtModel, levelConstraint, α * z1 + (1 - α) * y1 ≤ level);
        end
        add_constraint(currentInfo, nxtInfo);
        @objective(
            nxtModel, 
            Min, 
            sum(
                (xs[g] - x₀.ContVar[:s][g]) * (xs[g] - x₀.ContVar[:s][g]) + 
                (xy[g] - x₀.BinVar[:y][g]) * (xy[g] - x₀.BinVar[:y][g]) + 
                sum((sur[g, k] - x₀.ContAugState[:s][g][k]) * (sur[g, k] - x₀.ContAugState[:s][g][k]) for k in keys(x₀.ContAugState[:s][g]); init = 0.0) for g in G) 
        );
        optimize!(nxtModel);
        st = termination_status(nxtModel);
        if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
            BinVar = Dict{Any, Dict{Any, Any}}(:y => Dict{Any, Any}(
                g => JuMP.value(xy[g]) for g in indexSets.G)
            );
            ContVar = Dict{Any, Dict{Any, Any}}(:s => Dict{Any, Any}(
                g => JuMP.value(xs[g]) for g in indexSets.G)
            );
            ContAugState = Dict{Any, Dict{Any, Dict{Any, Any}}}(
                :s => Dict{Any, Dict{Any, Any}}(
                    g => Dict{Any, Any}(
                        k => JuMP.value(sur[g, k]) for k in keys(stateInfo.ContAugState[:s][g])
                    ) for g in indexSets.G
                )
            );

            x_nxt = StateInfo(
                BinVar, 
                nothing, 
                ContVar, 
                nothing, 
                nothing, 
                nothing, 
                nothing, 
                nothing, 
                ContAugState
            );
        elseif st == MOI.NUMERICAL_ERROR ## need to figure out why this case happened and fix it
            # @info "Numerical Error occurs -- Build a new nxtModel"
            nxtModel = Model(optimizer_with_attributes(
                ()->Gurobi.Optimizer(GRB_ENV), 
                "Threads" => 0)
            ); 
            MOI.set(nxtModel, MOI.Silent(), true);
            set_optimizer_attribute(nxtModel, "MIPGap", param.MIPGap);
            set_optimizer_attribute(nxtModel, "TimeLimit", param.TimeLimit);
            @variable(nxtModel, xs[G]);
            @variable(nxtModel, xy[G]);
            @variable(nxtModel, sur[g in G, k in keys(stateInfo.ContAugState[:s][g])]);
            @variable(nxtModel, z1);
            @variable(nxtModel, y1);
            nxtInfo = ModelInfo(nxtModel, xs, xy, sur, y1, z1);
            @constraint(nxtModel, levelConstraint, α * z1 + (1 - α) * y1 ≤ level);
            add_constraint(currentInfo, nxtInfo);
            @objective(
                nxtModel, 
                Min, 
                sum(
                    (xs[g] - x₀.ContVar[:s][g]) * (xs[g] - x₀.ContVar[:s][g]) + 
                    (xy[g] - x₀.BinVar[:y][g]) * (xy[g] - x₀.BinVar[:y][g]) + 
                    sum((sur[g, k] - x₀.ContAugState[:s][g][k]) * (sur[g, k] - x₀.ContAugState[:s][g][k]) for k in keys(x₀.ContAugState[:s][g]); init = 0.0) for g in G) 
            );
            optimize!(nxtModel);
            st = termination_status(nxtModel);
            if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
                BinVar = Dict{Any, Dict{Any, Any}}(:y => Dict{Any, Any}(
                    g => JuMP.value(xy[g]) for g in indexSets.G)
                );
                ContVar = Dict{Any, Dict{Any, Any}}(:s => Dict{Any, Any}(
                    g => JuMP.value(xs[g]) for g in indexSets.G)
                );
                ContAugState = Dict{Any, Dict{Any, Dict{Any, Any}}}(
                    :s => Dict{Any, Dict{Any, Any}}(
                        g => Dict{Any, Any}(
                            k => JuMP.value(sur[g, k]) for k in keys(stateInfo.ContAugState[:s][g])
                        ) for g in indexSets.G
                    )
                );

                x_nxt = StateInfo(
                    BinVar, 
                    nothing, 
                    ContVar, 
                    nothing, 
                    nothing, 
                    nothing, 
                    nothing, 
                    nothing, 
                    ContAugState
                );
            else
                return (cutInfo = cutInfo, iter = iter)
            end
        else
            # @info "Re-compute Next Iteration Point -- change to a safe level!"
            set_normalized_rhs( levelConstraint, w + .99 * (W - w))
            optimize!(nxtModel)
            st = termination_status(nxtModel);
            if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
                BinVar = Dict{Any, Dict{Any, Any}}(:y => Dict{Any, Any}(
                    g => JuMP.value(xy[g]) for g in indexSets.G)
                );
                ContVar = Dict{Any, Dict{Any, Any}}(:s => Dict{Any, Any}(
                    g => JuMP.value(xs[g]) for g in indexSets.G)
                );
                ContAugState = Dict{Any, Dict{Any, Dict{Any, Any}}}(
                    :s => Dict{Any, Dict{Any, Any}}(
                        g => Dict{Any, Any}(
                            k => JuMP.value(sur[g, k]) for k in keys(stateInfo.ContAugState[:s][g])
                        ) for g in indexSets.G
                    )
                );

                x_nxt = StateInfo(
                    BinVar, 
                    nothing, 
                    ContVar, 
                    nothing, 
                    nothing, 
                    nothing, 
                    nothing, 
                    nothing, 
                    ContAugState
                );
            else
                return (cutInfo = cutInfo, iter = iter)
            end
        end

        ## stop rules
        if Δ ≤ param_levelsetmethod.threshold * primal_bound || iter > param_levelsetmethod.MaxIter
            return (cutInfo = cutInfo, iter = iter)
        end
        ## ==================================================== end ============================================== ##
        ## save the trajectory
        currentInfo, currentInfo_f = solve_inner_minimization_problem(
            PLCGeneration,
            model, 
            x_nxt, 
            stateInfo
        );
        iter = iter + 1;
        functionHistory.f_his[iter] = currentInfo.f;
        functionHistory.G_max_his[iter] = maximum(currentInfo.G[k] for k in keys(currentInfo.G));
    end

end
