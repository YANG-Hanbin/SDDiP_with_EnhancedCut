"""
ModelModification!(; model::Model = model)

# Arguments

    1. `model::Model` : a forward pass model of stage t
    2. `randomVariables::RandomVariables` : random variables
    3. `paramDemand::ParamDemand` : demand parameters
    4. `stateInfo::StateInfo` : the last stage decisions
  
# Modification
    1. Remove the other scenario's demand balance constraints
    2. Add the current scenario's demand balance constraints
    3. Update its last stage decision with
"""
function ModelModification!( 
    model::Model, 
    randomVariables::RandomVariables,
    paramDemand::ParamDemand,
    stateInfo::StateInfo;
    indexSets::IndexSets = indexSets
)::Nothing
    if :ContVarNonAnticipative ∉ keys(model.obj_dict)
        @constraint(
            model, 
            ContVarNonAnticipative[g in indexSets.G], 
            model[:s_copy][g] == stateInfo.ContVar[:s][g]
        );
    end

    if :BinVarNonAnticipative ∉ keys(model.obj_dict)
        @constraint(
            model, 
            BinVarNonAnticipative[g in indexSets.G], 
            model[:y_copy][g] == stateInfo.BinVar[:y][g]
        );
    end

    # power balance constraints
    for i in indexSets.B
        delete(model, model[:PowerBalance][i])
    end
    unregister(model, :PowerBalance)
    @constraint(model, PowerBalance[i in indexSets.B], 
                            sum(model[:s][g]      for g in indexSets.Gᵢ[i]) -
                            sum(model[:P][(i, j)] for j in indexSets.out_L[i]) + 
                            sum(model[:P][(j, i)] for j in indexSets.in_L[i]) .==
                            sum(paramDemand.demand[d] * randomVariables.deviation[d] * model[:x][d] for d in indexSets.Dᵢ[i])
    )

    @objective(
        model, 
        Min, 
        model[:primal_objective_expression]
    );
    return
end

"""
function sample_scenarios(; 
    numScenarios::Int64 = 10, 
    scenarioTree::ScenarioTree = scenarioTree
)

# Arguments

    1. `numScenarios`: The number of scenarios will be sampled
    2. `scenarioTree`: A scenario tree

# Returns
    1. `Ξ`: A subset of scenarios.
"""
function sample_scenarios(; 
    numScenarios::Int64 = 10, 
    scenarioTree::ScenarioTree = scenarioTree
)
    Ξ = Dict{Int64, Dict{Int64, RandomVariables}}()
    for ω in 1:numScenarios
        ξ = Dict{Int64, RandomVariables}()
        ξ[1] = scenarioTree.tree[1].nodes[1]
        n = wsample(
            collect(keys(scenarioTree.tree[1].prob)), 
            collect(values(scenarioTree.tree[1].prob)), 
            1
        )[1]
        for t in 2:length(keys(scenarioTree.tree))
            ξ[t] = scenarioTree.tree[t].nodes[n]
            n = wsample(
                collect(keys(scenarioTree.tree[t].prob)), 
                collect(values(scenarioTree.tree[t].prob)), 
            1)[1]
        end
        Ξ[ω] = ξ
    end
    return Ξ
end