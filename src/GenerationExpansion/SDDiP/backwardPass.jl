#############################################################################################
##########################    auxiliary functions for backward   ############################
#############################################################################################
function add_generator_constraint(
    stageData::StageData, 
    modelInfo::BackwardModelInfo;
    binaryInfo::BinaryInfo = binaryInfo
)::Nothing
    ## no more than max num of generators
    @constraint(
        modelInfo.model, 
        binaryInfo.A * modelInfo.Lc + modelInfo.x .≤ stageData.ū 
    );
    # satisfy demand
    @constraint(
        modelInfo.model, 
        sum(modelInfo.y) + modelInfo.slack .≥ modelInfo.demand 
    );
    # no more than capacity  
    @constraint(
        modelInfo.model, 
        stageData.h * stageData.N * (binaryInfo.A * modelInfo.Lc + modelInfo.x + stageData.s₀ ) .≥ modelInfo.y 
    );  

    return 
end

function add_generator_cut(
    cutCoefficient::CutCoefficient, 
    modelInfo::BackwardModelInfo
)::Nothing
    iter = length(keys(cutCoefficient.v));  ## iter num
    k = length(keys(cutCoefficient.v[1]));  ## scenario num

    @constraint(
        modelInfo.model, 
        cut[i in 1:iter-1, m in 1:k], 
        modelInfo.θ ≥ cutCoefficient.v[i][m] + cutCoefficient.π[i][m]' * modelInfo.Lt 
    );
    return                               
end

#############################################################################################
###################################  function: backward pass ################################
#############################################################################################

"""
backwardModel!(stageData::StageData; 
    θ_bound::Real = 0.0, 
    binaryInfo::BinaryInfo = binaryInfo, 
    timelimit::Int64 = 3, 
    mipGap::Real = 1e-4, 
    tightness::Bool = false
)
"""
function backwardModel!(             
    stageData::StageData; 
    θ_bound::Real = 0.0,
    binaryInfo::BinaryInfo = binaryInfo, 
    timelimit::Int64 = 3, 
    mipGap::Real = 1e-4, 
    tightness::Bool = false
)::BackwardModelInfo


    model = Model(optimizer_with_attributes(
        ()->Gurobi.Optimizer(GRB_ENV), 
        "OutputFlag" => 0, 
        "Threads" => 0)
    ); 
    MOI.set(model, MOI.Silent(), true);
    set_optimizer_attribute(model, "MIPGap", mipGap);
    set_optimizer_attribute(model, "TimeLimit", timelimit);

    # the number of generators will be built in this stage
    @variable(model, x[i = 1:binaryInfo.d] ≥ 0, Int);                   
    # auxiliary variable (copy variable)
    if tightness
        @variable(model, Lc[i = 1:binaryInfo.n], Bin);          
    else
        @variable(model, 0 ≤ Lc[i = 1:binaryInfo.n]≤ 1);                      
    end
    # stage variable, A * Lt is total number of generators built after this stage
    @variable(model, Lt[i = 1:binaryInfo.n], Bin);                      
    @variable(model, y[i = 1:binaryInfo.d] ≥ 0);
    @variable(model, slack ≥ 0);
    @variable(model, θ ≥ θ_bound);
    ## no more than max num of generators
    @constraint(
        model, 
        binaryInfo.A * Lc + x .≤ stageData.ū
    );           
    ## no more than capacity            
    @constraint(
        model, 
        stageData.h * stageData.N * (binaryInfo.A * Lc + x + stageData.s₀ ) .≥ y 
    );   

    @constraint(
        model, 
        demandConstraint, 
        sum(y) + slack .≥ 0
    );
    ## to ensure pass a binary variable for next stage
    @constraint(
        model, 
        binaryInfo.A * Lc + x .== binaryInfo.A * Lt 
    );    
    
    # a cap constraint for two generators
    @constraint(
        model,
        [i = 4:5],
        y[i] ≤ sum(y)/5
    )
    return BackwardModelInfo(model, x, Lt, Lc, y, θ, slack)
end


function backward_Constraint_Modification!(
    backwardInfo::BackwardModelInfo,                                     
    demand::Vector{Float64}                 
)::Nothing

    delete(backwardInfo.model, backwardInfo.model[:demandConstraint])
    unregister(backwardInfo.model, :demandConstraint)
    # satisfy demand
    @constraint(
        backwardInfo.model, 
        demandConstraint, 
        sum(backwardInfo.y) + backwardInfo.slack .≥ demand 
    );
    return 
end

function copyVariable_Constraint!(
    backwardInfo::BackwardModelInfo, 
    stageData::StageData, 
    pdemand::Vector{Float64}; ## last-stage demand!
    binaryInfo::BinaryInfo = binaryInfo                  
)::Nothing

    delete(backwardInfo.model, backwardInfo.model[:hierarchicalConstraint])
    unregister(backwardInfo.model, :hierarchicalConstraint)
    # satisfy demand
    @constraint(
        backwardInfo.model, 
        hierarchicalConstraint, stageData.h * stageData.N * (binaryInfo.A * backwardInfo.Lc + stageData.s₀ ) .≥ pdemand
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
        stageDataList[t], 
        Ω[t][j].d, 
        solCollection[t-1,k].stageSolution;
        cutSelection = param.cutSelection,                    
        binaryInfo = binaryInfo,                         
        Output_Gap = param.Output_Gap,                                 
        λ = .3, 
        ℓ1 = param.ℓ1,                                              
        ℓ2 = param.ℓ2,
        nxt_bound = param.nxt_bound
    );

    (λ₀, λ₁) = LevelSetMethod_optimization!(
        backwardInfo, x₀; 
        levelSetMethodParam = levelSetMethodParam, 
        stageData = stageDataList[t], 
        ϵ = param.ε, 
        binaryInfo = binaryInfo
    )

                                                    
    return [λ₀, λ₁]
end
