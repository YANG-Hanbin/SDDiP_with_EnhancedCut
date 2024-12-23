export backwardPass
"""
RemoveContVarNonAnticipative!(model::Model)

# Arguments

    1. `model::Model` : a nodal problem
  
# Modification
    1. Remove the Non-anticipativity constraints
"""
function RemoveContVarNonAnticipative!(
    model::Model = model;
    indexSets::IndexSets = indexSets,
    param::NamedTuple = param
)::Nothing

    if :λ_copy ∉ keys(model.obj_dict)
        for g in indexSets.G
            delete(model, model[:ContVarNonAnticipative][g]);
            delete(model, model[:BinVarNonAnticipative][g]);
        end
        unregister(model, :ContVarNonAnticipative);
        unregister(model, :BinVarNonAnticipative);
    else
        for g in indexSets.G
            for i in 1:param.κ[g]
                delete(model, model[:BinarizationNonAnticipative][g, i]);
            end
        end
        unregister(model, :BinarizationNonAnticipative);
    end
    
    return
end

"""
    setup_initial_point(stateInfo::StateInfo)
# Arguments
    stateInfo::StateInfo : the parent's node decisions
    function for setting up the initial dual variables
"""
function setup_initial_point(
    stateInfo::StateInfo;
    indexSets::IndexSets = indexSets,
    param::NamedTuple = param
)::StateInfo
    BinVar = Dict{Any, Dict{Any, Any}}(:y => Dict{Any, Any}(
        g => 0.0 for g in indexSets.G)
    );
    ContVar = Dict{Any, Dict{Any, Any}}(:s => Dict{Any, Any}(
        g => 0.0 for g in indexSets.G)
    );
    if stateInfo.ContAugState == nothing 
        ContAugState = nothing
    else
        ContAugState = Dict{Any, Dict{Any, Dict{Any, Any}}}(
            :s => Dict{Any, Dict{Any, Any}}(
                g => Dict{Any, Any}(
                    k => 0.0 for k in keys(stateInfo.ContAugState[:s][g])
                ) for g in indexSets.G
            )
        );
    end

    if stateInfo.ContStateBin == nothing 
        ContStateBin = nothing
    else
        ContStateBin = Dict{Any, Dict{Any, Dict{Any, Any}}}(
            :s => Dict{Any, Dict{Any, Any}}(
                g => Dict{Any, Any}(
                    i => 0.0 for i in 1:param.κ[g]
                ) for g in indexSets.G
            )
        );
    end

    return StateInfo(
        BinVar, 
        nothing, 
        ContVar, 
        nothing, 
        nothing, 
        nothing, 
        nothing, 
        nothing, 
        ContAugState,
        nothing,
        ContStateBin
    );
end

"""
    backwardPass(backwardNodeInfo)

    function for backward pass in parallel computing
"""
function backwardPass(
    backwardNodeInfo::Tuple; 
    ModelList::Dict{Int64, SDDPModel} = ModelList,
    indexSets::IndexSets = indexSets, 
    paramDemand::ParamDemand = paramDemand, 
    paramOPF::ParamOPF = paramOPF,
    scenarioTree::ScenarioTree = scenarioTree, 
    stateInfoCollection::Dict{Any, Any} = stateInfoCollection,
    param::NamedTuple = param, param_cut::NamedTuple = param_cut, param_levelsetmethod::NamedTuple = param_levelsetmethod
)

    (i, t, n, ω, cutSelection, core_point_strategy) = backwardNodeInfo; 
    ModelModification!( 
        ModelList[t].model, 
        scenarioTree.tree[t].nodes[n],
        paramDemand,
        stateInfoCollection[i, t-1, ω];
        indexSets = indexSets
    );
    RemoveContVarNonAnticipative!(
        ModelList[t].model;
        indexSets = indexSets,
        param = param
    );

    if cutSelection == :PLC
        CutGenerationInfo = ParetoLagrangianCutGeneration{Float64}(
            core_point_strategy, 
            setup_core_point(
                stateInfoCollection[i, t-1, ω];
                indexSets = indexSets,
                paramOPF = paramOPF, 
                param_cut = param_cut,
                param = param
            ), 
            param_cut.δ, 
            stateInfoCollection[i, t, ω].StateValue
        );
    elseif cutSelection == :LC
        CutGenerationInfo = LagrangianCutGeneration{Float64}(
            stateInfoCollection[i, t, ω].StateValue
        );
    elseif cutSelection == :SMC
        CutGenerationInfo = SquareMinimizationCutGeneration{Float64}(
            param_cut.δ, 
            stateInfoCollection[i, t, ω].StateValue
        );
    else
        @warn "Invalid cutSelection value: $cutSelection. Defaulting to :SMC."
        CutGenerationInfo = SquareMinimizationCutGeneration{Float64}(
            param_cut.δ, 
            stateInfoCollection[i, t, ω].StateValue
        )
    end

    levelsetmethodOracleParam = SetupLevelSetMethodOracleParam(
        stateInfoCollection[i, t-1, ω];
        indexSets = indexSets,
        param_levelsetmethod = param_levelsetmethod
    );
    
    ((λ₀, λ₁), LMiter) = LevelSetMethod_optimization!(
        ModelList[t].model, 
        levelsetmethodOracleParam, 
        stateInfoCollection[i, t-1, ω],
        CutGenerationInfo;
        indexSets = indexSets, 
        paramDemand = paramDemand, 
        paramOPF = paramOPF, 
        param = param, param_levelsetmethod = param_levelsetmethod
    );
                                                    
    # f_star_value = λ₀ + 
    #     sum(
    #         λ₁.ContVar[:s][g] * stateInfoCollection[i, t-1, ω].ContVar[:s][g] + 
    #         λ₁.BinVar[:y][g] * stateInfoCollection[i, t-1, ω].BinVar[:y][g] + 
    #         sum(λ₁.ContAugState[:s][g][k] * stateInfoCollection[i, t-1, ω].ContAugState[:s][g][k] for k in keys(stateInfoCollection[i, t-1, ω].ContAugState[:s][g]); init = 0.0) 
    #         for g in indexSets.G
    # );

    return ((λ₀, λ₁), LMiter)  
end