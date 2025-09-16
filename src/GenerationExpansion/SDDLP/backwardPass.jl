#############################################################################################
###################################  function: backward pass ################################
#############################################################################################

"""
    This is the oracle in level set method, and it will return [model, dF]
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

    @variable(model, x[g = 1:binaryInfo.d] ≥ 0, Int)                   # the number of generators will be built in this stage
    @variable(model, y[g = 1:binaryInfo.d] ≥ 0)
    @variable(model, St[g = 1:binaryInfo.d] ≥ 0, Int)                  # stage variable, A * Lt is total number of generators built after this stage
    @variable(model, slack ≥ 0 )
    @variable(model, θ ≥ θ_bound)


    # @variable(model, sur[G, 1:max_sur], Bin)                        
    sur = Dict(
        (g, i) => @variable(model, base_name = "sur[$g, $i]", binary = true)
        for g in 1:binaryInfo.d for i in 1:1
    );
    model[:sur] = sur;                                              ## sur[g, k] is the kth surrogate variable of s[g]

    # constraints for surrogate variables
    ## Choosing one leaf node
    @constraint(model, [g in 1:binaryInfo.d], sur[g, 1] == 1)

    # auxiliary variable (copy variable) for S_{t-1} or Ŝ
    if tightness
        @variable(model, Sc[i = 1:binaryInfo.d] ≥ 0, Int)   
        @constraint(model, Sc .≤ stageData.ū)                
        sur_copy = Dict(
            (g, i) => @variable(model, base_name = "sur_copy[$g, $i]", binary = true)
            for g in 1:binaryInfo.d for i in 1:1
        );
        model[:sur_copy] = sur_copy;     
    else
        @variable(model, Sc[i = 1:binaryInfo.d] ≥ 0)   
        @constraint(model, Sc .≤ stageData.ū)                
        sur_copy = Dict(
            (g, i) => @variable(model, base_name = "sur_copy[$g, $i]", lower_bound = 0, upper_bound = 1)
            for g in 1:binaryInfo.d for i in 1:1
        );
        model[:sur_copy] = sur_copy;     
    end
    @constraint(model, [g in 1:binaryInfo.d], sur_copy[g, 1] == 1)


    @constraint(model, St .≤ stageData.ū )                       ## no more than max num of generators
    @constraint(model, stageData.h * stageData.N 
                * (St + stageData.s₀ ) .≥ y )                ## no more than capacity

    @constraint(model, demandConstraint, sum(y) + slack .≥ 0);
    ## to ensure pass a binary variable for next stage
    @constraint(model, Sc + x .== St);     
    
    # a cap constraint for two generators
    @constraint(
        model,
        [i = 4:5],
        y[i] ≤ sum(y)/5
    )
    

    return BackwardModelInfo(model, x, St, Sc, y, θ, slack)
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
    pdemand::Vector{Float64}                 
)::Nothing

    delete(backwardInfo.model, backwardInfo.model[:hierarchicalConstraint])
    unregister(backwardInfo.model, :hierarchicalConstraint)

    # satisfy demand
    @constraint(
        backwardInfo.model, 
        hierarchicalConstraint, 
        stageData.h * stageData.N * (Sc + stageData.s₀ ) .≥ pdemand
    );    
    return
end

function get_Benders_coefficient!(
    backwardInfo::BackwardModelInfo, 
    stateInfo::NamedTuple;
    binaryInfo::BinaryInfo = binaryInfo
)::Dict
    @constraint(backwardInfo.model, StateNonAnticipativity, backwardInfo.Sc .== stateInfo.stageSolution);  
    @constraint(backwardInfo.model, AugNonAnticipativity[g in 1:binaryInfo.d, i in keys(stateInfo.stageSur[g])], backwardInfo.model[:sur_copy][g, i] .== stateInfo.stageSur[g][i]);  
    lp_model = relax_integrality(backwardInfo.model);
    optimize!(backwardInfo.model);
    # obtain the Benders cut coefficients
    x₀ = Dict( 
        :St => dual.(backwardInfo.model[:StateNonAnticipativity]),                 
        :sur => Dict(
            g => Dict(
                i => dual.(backwardInfo.model[:AugNonAnticipativity][g, i]) for i in keys(stateInfo.stageSur[g])
            ) for g in 1:binaryInfo.d
        )
    );
    delete(backwardInfo.model, backwardInfo.model[:StateNonAnticipativity])
    for g in 1:binaryInfo.d, i in keys(stateInfo.stageSur[g])
        delete(backwardInfo.model, backwardInfo.model[:AugNonAnticipativity][g, i])
    end
    unregister(backwardInfo.model, :StateNonAnticipativity)
    unregister(backwardInfo.model, :AugNonAnticipativity)

    lp_model();

    return x₀
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
        backwardInfo,
        stageDataList[t], 
        Ω[t][j].d, 
        solCollection[t-1,k], 
        param;
        binaryInfo = binaryInfo, 
        Output_Gap = param.Output_Gap,
        λ = .5
    );

    (λ₀, λ₁) = LevelSetMethod_optimization!(
        backwardInfo, 
        x₀; 
        levelSetMethodParam = levelSetMethodParam, 
        stageData = stageDataList[t],     
        binaryInfo = binaryInfo,
        param = param
    ); 
                                               
    return [λ₀, λ₁]
end

