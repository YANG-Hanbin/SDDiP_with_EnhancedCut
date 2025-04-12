#############################################################################################
###################################  function: backward pass ################################
#############################################################################################

"""
    function backwardModel!( 
        stageData::StageData, 
        param::NamedTuple;                                   
        θ_bound::Real = 0.0,                  
        binaryInfo::BinaryInfo = binaryInfo,                     
        tightness::Bool = false
    )::BackwardModelInfo
    
    This function is used to build the backward model for SDDP
"""
function backwardModel!( 
    stageData::StageData, 
    param::NamedTuple;                                   
    θ_bound::Real = 0.0,                  
    binaryInfo::BinaryInfo = binaryInfo,                     
    tightness::Bool = false
)::BackwardModelInfo


    model = Model(optimizer_with_attributes(
        () -> Gurobi.Optimizer(GRB_ENV), 
        "Threads" => 0)); 
    MOI.set(model, MOI.Silent(), !param.verbose);
    set_optimizer_attribute(model, "MIPGap", param.MIPGap);
    set_optimizer_attribute(model, "TimeLimit", param.TimeLimit);
    set_optimizer_attribute(model, "MIPFocus", param.MIPFocus);           
    set_optimizer_attribute(model, "FeasibilityTol", param.FeasibilityTol);

    # the number of generators will be built in this stage
    @variable(model, x[i = 1:binaryInfo.d] ≥ 0, Int)                   
    @variable(model, y[i = 1:binaryInfo.d] ≥ 0)
    # stage variable, A * Lt is total number of generators built after this stage
    @variable(model, St[i = 1:binaryInfo.d] ≥ 0, Int)                      
    @variable(model, slack ≥ 0 )
    @variable(model, θ ≥ θ_bound)
    if tightness
        @variable(model, Sc[i = 1:binaryInfo.d] ≥ 0, Int)   
        @constraint(model, Sc .≤ stageData.ū)                  
    else
        @variable(model, Sc[i = 1:binaryInfo.d] ≥ 0)    
        @constraint(model, Sc .≤ stageData.ū)                  
    end

    ## no more than max num of generators
    @constraint(model, St .≤ stageData.ū )                       
    ## no more than capacity
    @constraint(
        model, 
        stageData.h * stageData.N * (St + stageData.s₀ ) .≥ y 
    );                

    @constraint(model, demandConstraint, sum(y) + slack .≥ 0)
    ## to ensure pass a binary variable for next stage
    @constraint(model, Sc + x .== St )  
    
    # a cap constraint for two generators
    @constraint(
        model,
        [i = 4:5],
        y[i] ≤ sum(y)/5
    )

    return BackwardModelInfo(model, x, St, Sc, y, θ, slack)
end


"""
    This function is used to modify the constraint in backward pass
    to satisfy the demand
"""
function backward_Constraint_Modification!(
    backwardInfo::BackwardModelInfo, 
    demand::Vector{Float64}                     
)::Nothing

    delete(backwardInfo.model, backwardInfo.model[:demandConstraint])
    unregister(backwardInfo.model, :demandConstraint)

    # satisfy demand
    @constraint(
        backwardInfo.model, 
        demandConstraint, sum(backwardInfo.y) + backwardInfo.slack .≥ demand 
    );
    return 
end

"""
    This function is used to modify the constraint in backward pass
    to satisfy the demand
"""
function copyVariable_Constraint!(
    backwardInfo::BackwardModelInfo, 
    stageData::StageData, 
    pdemand::Vector{Float64}                 
)::Nothing

    delete(backwardInfo.model, backwardInfo.model[:hierarchicalConstraint])
    unregister(backwardInfo.model, :hierarchicalConstraint)

    # satisfy demand
    @constraint(
        backwardInfo.model, 
        hierarchicalConstraint, 
        stageData.h * stageData.N * (Sc + stageData.s₀) .≥ pdemand
    ); 
    return
end

"""
    backwardPass(backwardNodeInfo)

    function for backward pass in parallel computing
"""
function backwardPass(
    backwardNodeInfo::Tuple, 
    solCollection::Dict{Any, Any}; 
    stageDataList::Dict{Int64, StageData} = stageDataList,
    backwardInfoList::Dict{Int64, BackwardModelInfo} = backwardInfoList,
    Ω::Dict{Int64, Dict{Int64, RandomVariables}} = Ω,
    param::NamedTuple = param,
    binaryInfo::BinaryInfo = binaryInfo,
)

    (t, j, k) = backwardNodeInfo; 
    backwardInfo = backwardInfoList[t]
    backward_Constraint_Modification!(
        backwardInfo,                    
        Ω[t][j].d
    );
    (levelSetMethodParam, x₀) = setupLevelsetPara(
        forwardInfoList[t], 
        stageDataList[t], Ω[t][j].d, 
        solCollection[t-1,k].stageSolution;
        cutSelection = param.cutSelection,  # "ShrinkageLC", "LC", "ELC"                       
        binaryInfo = binaryInfo, 
        Output_Gap = param.Output_Gap, 
        max_iter = param.MaxIter,                                       
        λ = .3, ℓ1 = 0.0, ℓ2 = 1.
    );
                        
    (λ₀, λ₁) = LevelSetMethod_optimization!(
        backwardInfo, 
        x₀; 
        levelSetMethodParam = levelSetMethodParam, 
        stageData = stageDataList[t],   
        ϵ = param.ε, 
        binaryInfo = binaryInfo
    ); 
                                                    
    return [λ₀, λ₁]
end
