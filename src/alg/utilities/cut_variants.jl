"""
    function solve_inner_minimization_problem(
        CutGenerationInfo::ParetoLagrangianCutGeneration,
        model::Model, 
        x₀::StateInfo, 
        stateInfo::StateInfo
    )

# Arguments

    1. `CutGenerationInfo::ParetoLagrangianCutGeneration` : the information of the cut that will be generated information
    2. `model::Model` : the backward model
    3. `x₀::StateInfo` : the dual information
    4. `stateInfo::StateInfo` : the last stage decision
  
# Returns
    1. `currentInfo::CurrentInfo` : the current information
  
"""
function solve_inner_minimization_problem(
    CutGenerationInfo::ParetoLagrangianCutGeneration,
    model::Model, 
    x₀::StateInfo, 
    stateInfo::StateInfo;
    indexSets::IndexSets = indexSets,
    param::NamedTuple = param
)
    @objective(
        model, 
        Min,  
        model[:primal_objective_expression] +
        sum( 
            x₀.BinVar[:y][g] * (stateInfo.BinVar[:y][g] - model[:y_copy][g]) + 
            (stateInfo.ContStateBin !== nothing ?
                sum(
                    x₀.ContStateBin[:s][g][i] * (stateInfo.ContStateBin[:s][g][i] - model[:λ_copy][g, i]) for i in 1:param.κ[g]
                ) : x₀.ContVar[:s][g] * (stateInfo.ContVar[:s][g] - model[:s_copy][g])
            ) +
            (stateInfo.ContAugState !== nothing && haskey(stateInfo.ContAugState, :s) && haskey(stateInfo.ContAugState[:s], g) ?
                sum(
                    x₀.ContAugState[:s][g][k] * (stateInfo.ContAugState[:s][g][k] - model[:augmentVar_copy][g, k]) 
                    for k in keys(stateInfo.ContAugState[:s][g]); init = 0.0
                ) : 0.0
            ) 
            for g in indexSets.G
        )
    );
    ## ==================================================== solve the model and display the result ==================================================== ##
    optimize!(model);
    F  = JuMP.objective_value(model);
    
    negative_∇F = StateInfo(
        Dict{Any, Dict{Any, Any}}(
            :y => Dict{Any, Any}(
                g => JuMP.value(model[:y_copy][g]) - stateInfo.BinVar[:y][g] for g in indexSets.G
            )
        ), 
        nothing, 
        Dict{Any, Dict{Any, Any}}(
            :s => Dict{Any, Any}(
                g => JuMP.value(model[:s_copy][g]) - stateInfo.ContVar[:s][g] for g in indexSets.G
            )
        ), 
        nothing, 
        nothing, 
        nothing, 
        nothing, 
        nothing, 
        stateInfo.ContAugState !== nothing ? Dict{Any, Dict{Any, Dict{Any, Any}}}(
            :s => Dict{Any, Dict{Any, Any}}(
                g => Dict{Any, Any}(
                        k => JuMP.value(model[:augmentVar_copy][g, k]) - stateInfo.ContAugState[:s][g][k]
                        for k in keys(stateInfo.ContAugState[:s][g])
                    ) 
                for g in indexSets.G
            )
        ) : nothing,
        nothing,
        stateInfo.ContStateBin !== nothing ? Dict{Any, Dict{Any, Dict{Any, Any}}}(
            :s => Dict{Any, Dict{Any, Any}}(
                g => Dict{Any, Any}(
                        i => JuMP.value(model[:λ_copy][g, i]) - stateInfo.ContStateBin[:s][g][i]
                        for i in 1:param.κ[g]
                    ) 
                for g in indexSets.G
            )
        ) : nothing
    );
    currentInfo = CurrentInfo(  
        x₀, 
        - F - sum(
            (stateInfo.ContStateBin !== nothing ?
                sum(
                    x₀.ContStateBin[:s][g][i] * (CutGenerationInfo.core_point.ContStateBin[:s][g][i] - stateInfo.ContStateBin[:s][g][i]) for i in 1:param.κ[g]
                ) : x₀.ContVar[:s][g] * (CutGenerationInfo.core_point.ContVar[:s][g] - stateInfo.ContVar[:s][g])
            ) +
            x₀.BinVar[:y][g] * (CutGenerationInfo.core_point.BinVar[:y][g] - stateInfo.BinVar[:y][g]) +
            (
                stateInfo.ContAugState !== nothing && haskey(stateInfo.ContAugState, :s) && haskey(stateInfo.ContAugState[:s], g) ?
                sum(
                    x₀.ContAugState[:s][g][k] * (CutGenerationInfo.core_point.ContAugState[:s][g][k] - stateInfo.ContAugState[:s][g][k])
                    for k in keys(stateInfo.ContAugState[:s][g]); init = 0.0
                ) : 0.0
            )
            for g in indexSets.G
        ),                                                                                                                                                              ## obj function value
        Dict(
            1 => CutGenerationInfo.primal_bound - F - CutGenerationInfo.δ
        ),                                                                                                                                                              ## constraint value
        stateInfo.ContAugState !== nothing ? 
        Dict{Symbol, Dict{Int64, Any}}(
            :s => Dict(g => JuMP.value(model[:s_copy][g]) - CutGenerationInfo.core_point.ContVar[:s][g] for g in indexSets.G),
            :y => Dict(g => JuMP.value(model[:y_copy][g]) - CutGenerationInfo.core_point.BinVar[:y][g] for g in indexSets.G),
            :sur => Dict(
                g => Dict(
                    k => JuMP.value(model[:augmentVar_copy][g, k]) - CutGenerationInfo.core_point.ContAugState[:s][g][k] 
                            for k in keys(stateInfo.ContAugState[:s][g])
                ) for g in indexSets.G
            ),
            :λ => Dict(
                g => (stateInfo.ContStateBin !== nothing ?
                    Dict(
                        i => JuMP.value(model[:λ_copy][g, i]) - CutGenerationInfo.core_point.ContStateBin[:s][g][i] for i in 1:param.κ[g]
                    ) : nothing
                ) for g in indexSets.G
            )
        ) :
        Dict{Symbol, Dict{Int64, Any}}(
            :s => Dict(g => JuMP.value(model[:s_copy][g]) - CutGenerationInfo.core_point.ContVar[:s][g] for g in indexSets.G),
            :y => Dict(g => JuMP.value(model[:y_copy][g]) - CutGenerationInfo.core_point.BinVar[:y][g] for g in indexSets.G),
            :λ => Dict(
                g => (stateInfo.ContStateBin !== nothing ?
                    Dict(
                        i => JuMP.value(model[:λ_copy][g, i]) - CutGenerationInfo.core_point.ContStateBin[:s][g][i] for i in 1:param.κ[g]
                    ) : nothing
                ) for g in indexSets.G
            )
        ),                                                                                                                                                              ## obj gradient
        Dict(1 => negative_∇F )                                                                                                                                         ## constraint gradient
    );
    return (currentInfo = currentInfo, currentInfo_f = F)
end    

"""
    function solve_inner_minimization_problem(
        CutGenerationInfo::SquareMinimizationCutGeneration,
        model::Model, 
        x₀::StateInfo, 
        stateInfo::StateInfo
    )

# Arguments

    1. `CutGenerationInfo::SquareMinimizationCutGeneration` : the information of the cut that will be generated information
    2. `model::Model` : the backward model
    3. `x₀::StateInfo` : the dual information
    4. `stateInfo::StateInfo` : the last stage decision
  
# Returns
    1. `currentInfo::CurrentInfo` : the current information
  
"""
function solve_inner_minimization_problem(
    CutGenerationInfo::SquareMinimizationCutGeneration,
    model::Model, 
    x₀::StateInfo, 
    stateInfo::StateInfo;
    indexSets::IndexSets = indexSets
)
    @objective(
        model, 
        Min,  
        model[:primal_objective_expression] +
        sum(
            x₀.BinVar[:y][g] * (stateInfo.BinVar[:y][g] - model[:y_copy][g]) + 
            (stateInfo.ContStateBin !== nothing ?
                sum(
                    x₀.ContStateBin[:s][g][i] * (stateInfo.ContStateBin[:s][g][i] - model[:λ_copy][g, i]) for i in 1:param.κ[g]
                ) : x₀.ContVar[:s][g] * (stateInfo.ContVar[:s][g] - model[:s_copy][g])
            ) +
            (stateInfo.ContAugState !== nothing && haskey(stateInfo.ContAugState, :s) && haskey(stateInfo.ContAugState[:s], g) ?
                sum(
                    x₀.ContAugState[:s][g][k] * (stateInfo.ContAugState[:s][g][k] - model[:augmentVar_copy][g, k]) 
                    for k in keys(stateInfo.ContAugState[:s][g]); init = 0.0
                ) : 0.0
            ) 
            for g in indexSets.G
        )
    );
    ## ==================================================== solve the model and display the result ==================================================== ##
    optimize!(model);
    F  = JuMP.objective_value(model);
    
    negative_∇F = StateInfo(
        Dict{Any, Dict{Any, Any}}(
            :y => Dict{Any, Any}(
                g => JuMP.value(model[:y_copy][g]) - stateInfo.BinVar[:y][g] for g in indexSets.G
            )
        ), 
        nothing, 
        Dict{Any, Dict{Any, Any}}(
            :s => Dict{Any, Any}(
                g => JuMP.value(model[:s_copy][g]) - stateInfo.ContVar[:s][g] for g in indexSets.G
            )
        ), 
        nothing, 
        nothing, 
        nothing, 
        nothing, 
        nothing, 
        stateInfo.ContAugState !== nothing ? Dict{Any, Dict{Any, Dict{Any, Any}}}(
            :s => Dict{Any, Dict{Any, Any}}(
                g => Dict{Any, Any}(
                        k => JuMP.value(model[:augmentVar_copy][g, k]) - stateInfo.ContAugState[:s][g][k]
                        for k in keys(stateInfo.ContAugState[:s][g])
                    ) 
                for g in indexSets.G
            )
        ) : nothing,
        nothing,
        stateInfo.ContStateBin !== nothing ? Dict{Any, Dict{Any, Dict{Any, Any}}}(
            :s => Dict{Any, Dict{Any, Any}}(
                g => Dict{Any, Any}(
                        i => JuMP.value(model[:λ_copy][g, i]) - stateInfo.ContStateBin[:s][g][i]
                        for i in 1:param.κ[g]
                    ) 
                for g in indexSets.G
            )
        ) : nothing
    );
    currentInfo = CurrentInfo(  
        x₀, 
        1/2 * (stateInfo.ContStateBin !== nothing ?
                sum(
                    x₀.ContStateBin[:s][g][i] * x₀.ContStateBin[:s][g][i] for i in 1:param.κ[g] for g in indexSets.G
                ) : sum(x₀.ContVar[:s][g] * x₀.ContVar[:s][g] for g in indexSets.G)
        ) +
        1/2 * sum(x₀.BinVar[:y][g] * x₀.BinVar[:y][g] for g in indexSets.G) + 
        1/2 * sum(sum(x₀.ContAugState[:s][g][k] * x₀.ContAugState[:s][g][k] for k in keys(stateInfo.ContAugState[:s][g]); init = 0.0) for g in indexSets.G
        ),                                                                                                                                                              ## obj function value
        Dict(
            1 => CutGenerationInfo.primal_bound - F - CutGenerationInfo.δ
        ),                                                                                                                                                              ## constraint value
        stateInfo.ContAugState !== nothing ?
        Dict{Symbol, Dict{Int64, Any}}(
            :s => Dict(g => x₀.ContVar[:s][g] for g in indexSets.G),
            :y => Dict(g => x₀.BinVar[:y][g] for g in indexSets.G), 
            :sur => Dict(g => Dict(k => x₀.ContAugState[:s][g][k] for k in keys(stateInfo.ContAugState[:s][g])) for g in indexSets.G),
            :λ => Dict(
                g => (
                    stateInfo.ContStateBin !== nothing ?
                        Dict(i => x₀.ContStateBin[:s][g][i] for i in 1:param.κ[g]) : nothing
                ) for g in indexSets.G
            )
        ) : 
        Dict{Symbol, Dict{Int64, Any}}(
            :s => Dict(g => x₀.ContVar[:s][g] for g in indexSets.G),
            :y => Dict(g => x₀.BinVar[:y][g] for g in indexSets.G),
            :λ => Dict(
                g => (
                    stateInfo.ContStateBin !== nothing ?
                        Dict(i => x₀.ContStateBin[:s][g][i] for i in 1:param.κ[g]) : nothing
                ) for g in indexSets.G
            )
        ),                                                                                                                                                              ## obj gradient
        Dict(1 => negative_∇F )                                                                                                                                         ## constraint gradient
    );
    return (currentInfo = currentInfo, currentInfo_f = F)
end    

"""
    function solve_inner_minimization_problem(
        CutGenerationInfo::LagrangianCutGeneration,
        model::Model, 
        x₀::StateInfo, 
        stateInfo::StateInfo
    )

# Arguments

    1. `CutGenerationInfo::LagrangianCutGeneration` : the information of the cut that will be generated information
    2. `model::Model` : the backward model
    3. `x₀::StateInfo` : the dual information
    4. `stateInfo::StateInfo` : the last stage decision
  
# Returns
    1. `currentInfo::CurrentInfo` : the current information
  
"""
function solve_inner_minimization_problem(
    CutGenerationInfo::LagrangianCutGeneration,
    model::Model, 
    x₀::StateInfo, 
    stateInfo::StateInfo;
    indexSets::IndexSets = indexSets
)
    @objective(
        model, 
        Min,  
        model[:primal_objective_expression] +
        sum(
            x₀.BinVar[:y][g] * (stateInfo.BinVar[:y][g] - model[:y_copy][g]) + 
            (stateInfo.ContStateBin !== nothing ?
                sum(
                    x₀.ContStateBin[:s][g][i] * (stateInfo.ContStateBin[:s][g][i] - model[:λ_copy][g, i]) for i in 1:param.κ[g]
                ) : x₀.ContVar[:s][g] * (stateInfo.ContVar[:s][g] - model[:s_copy][g])
            ) +
            (stateInfo.ContAugState !== nothing && haskey(stateInfo.ContAugState, :s) && haskey(stateInfo.ContAugState[:s], g) ?
                sum(
                    x₀.ContAugState[:s][g][k] * (stateInfo.ContAugState[:s][g][k] - model[:augmentVar_copy][g, k]) 
                    for k in keys(stateInfo.ContAugState[:s][g]); init = 0.0
                ) : 0.0
            ) 
            for g in indexSets.G
        )
    );
    ## ==================================================== solve the model and display the result ==================================================== ##
    optimize!(model);
    F  = JuMP.objective_value(model);
    
    currentInfo = CurrentInfo(  
        x₀, 
        - F,                                                                                                                                                            ## obj function value
        Dict(
            1 => 0.0
        ),                                                                                                                                                              ## constraint value
        stateInfo.ContAugState !== nothing ?                                                                                                                            
        Dict{Symbol, Dict{Int64, Any}}(
            :s   => Dict(g => JuMP.value(model[:s_copy][g]) - stateInfo.ContVar[:s][g] for g in indexSets.G),
            :y   => Dict(g => JuMP.value(model[:y_copy][g]) - stateInfo.BinVar[:y][g]  for g in indexSets.G), 
            :sur => Dict(g => Dict(
                k => JuMP.value(model[:augmentVar_copy][g, k]) - stateInfo.ContAugState[:s][g][k] for k in keys(stateInfo.ContAugState[:s][g])
                ) for g in indexSets.G),
            :λ   => Dict(
                g => (stateInfo.ContStateBin !== nothing ?
                    Dict(i => JuMP.value(model[:λ_copy][g, i]) - stateInfo.ContStateBin[:s][g][i] for i in 1:param.κ[g]) : nothing
                ) for g in indexSets.G
            )
        ) : 
        Dict{Symbol, Dict{Int64, Any}}(
            :s   => Dict(g => JuMP.value(model[:s_copy][g]) - stateInfo.ContVar[:s][g] for g in indexSets.G),
            :y   => Dict(g => JuMP.value(model[:y_copy][g]) - stateInfo.BinVar[:y][g]  for g in indexSets.G),
            :λ   => Dict(
                g => (stateInfo.ContStateBin !== nothing ?
                    Dict(i => JuMP.value(model[:λ_copy][g, i]) - stateInfo.ContStateBin[:s][g][i] for i in 1:param.κ[g]) : nothing
                ) for g in indexSets.G
            )
        ),                                                                                                                                                              ## obj gradient
        Dict(1 => StateInfo(
            Dict{Any, Dict{Any, Any}}(:y => Dict{Any, Any}(
                g => 0.0 for g in indexSets.G)
            ), 
            nothing, 
            Dict{Any, Dict{Any, Any}}(:s => Dict{Any, Any}(
                g => 0.0 for g in indexSets.G)
            ), 
            nothing, 
            nothing, 
            nothing, 
            nothing, 
            nothing, 
            stateInfo.ContAugState !== nothing ? Dict{Any, Dict{Any, Dict{Any, Any}}}(
                :s => Dict{Any, Dict{Any, Any}}(
                    g => Dict{Any, Any}(
                            k => 0.0
                            for k in keys(stateInfo.ContAugState[:s][g])
                        ) 
                    for g in indexSets.G
                )
            ) : nothing,
            nothing,
            stateInfo.ContStateBin !== nothing ? Dict{Any, Dict{Any, Dict{Any, Any}}}(
            :s => Dict{Any, Dict{Any, Any}}(
                g => Dict{Any, Any}(
                        i => 0.0
                        for i in 1:param.κ[g]
                    ) 
                for g in indexSets.G
            )
        ) : nothing
            )
        )                                                                                                                                                               ## constraint gradient
    );
    return (currentInfo = currentInfo, currentInfo_f = F)
end    