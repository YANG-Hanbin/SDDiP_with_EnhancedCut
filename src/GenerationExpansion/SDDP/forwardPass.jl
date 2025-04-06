#############################################################################################
###################################  function: forward pass #################################
#############################################################################################

"""     forward_step_optimize(StageCoefficient[t], b, x_ancestor, cut_collection, t)
        1. Solve the forward problem and return the useful info

        cutCoefficient: is the cut info given stage t
        stageData: is the param info given stage t
"""
function forwardModel!(
    stageData::StageData, 
    param::NamedTuple;       
    θ_bound::Real = 0.0,                                 
    binaryInfo::BinaryInfo = binaryInfo
)::ForwardModelInfo
                            
    ## construct forward problem (3.1)
    model = Model(
        optimizer_with_attributes(
            () -> Gurobi.Optimizer(GRB_ENV), 
            "Threads" => 0
        )
    ); 
    MOI.set(model, MOI.Silent(), !param.verbose);
    set_optimizer_attribute(model, "MIPGap", param.MIPGap);
    set_optimizer_attribute(model, "TimeLimit", param.TimeLimit);
    # set_optimizer_attribute(model, "MIPFocus", param.MIPFocus);           
    # set_optimizer_attribute(model, "FeasibilityTol", param.FeasibilityTol);
    
    @variable(model, x[i = 1:binaryInfo.d] ≥ 0, Int)        ## for current state, x is the number of generators will be built in this stage
    @variable(model, y[i = 1:binaryInfo.d] ≥ 0)             ## amount of electricity
    @variable(model, St[i = 1:binaryInfo.d] ≥ 0, Int)       ## stage variable, A * Lt is total number of generators built after this stage
    @variable(model, slack ≥ 0 )
    @variable(model, θ ≥ θ_bound)

    @objective(
        model, 
        Min, 
        stageData.c1'* x + stageData.c2' * y + stageData.penalty * slack + θ 
    );

    ## no more than max num of generators
    @constraint(
        model, 
        limitationConstraint,    
        St .≤ stageData.ū
    ); 

    # satisfy demand
    @constraint(
        model, 
        demandConstraint,        
        sum(y) + slack .≥ 0 
    );

    # no more than capacity
    @constraint(
        model, 
        capacityConstraint,      
        stageData.h * stageData.N * (St + stageData.s₀ ) .≥ y 
    );  

    @constraint(
        model, 
        totalGenerators,         
        0.0 .+ x .== St
    );

    return ForwardModelInfo(model, x, St, y, θ, slack)
end


"""
    function forward_modify_constraints!(
        forwardInfo::ForwardModelInfo, 
        stageData::StageData, 
        demand::Vector{Float64}, 
        Ŝ::Vector{Float64};
        binaryInfo::BinaryInfo = binaryInfo 
    )::Nothing
    1. modify the forward problem (3.1) for each stage t
    2. update the demand constraint and total generator constraint
    3. update the demand and total generator constraint
"""
function forward_modify_constraints!(
    forwardInfo::ForwardModelInfo, 
    stageData::StageData, 
    demand::Vector{Float64}, 
    Ŝ::Vector{Float64};
    binaryInfo::BinaryInfo = binaryInfo 
)::Nothing

    delete(
        forwardInfo.model, 
        forwardInfo.model[:demandConstraint]
    );
    
    delete(
        forwardInfo.model, 
        forwardInfo.model[:totalGenerators]
    );

    unregister(forwardInfo.model, :demandConstraint)
    unregister(forwardInfo.model, :totalGenerators)

    # satisfy demand
    @constraint(
        forwardInfo.model, 
        demandConstraint,        
        sum(forwardInfo.y) + forwardInfo.slack .≥ demand 
    );

    @constraint(
        forwardInfo.model, 
        totalGenerators,         
        Ŝ + forwardInfo.x .== forwardInfo.St
    );

    return 
end

"""
forwardPass(k): function for forward pass in parallel computing

# Arguments

  1. `k`: A sampled scenario path

# Returns
  1. `solCollection`: forward pass solution collection
"""
function forwardPass(
    k::Int64, 
    Scenarios::Dict{Int64, Vector{Int64}};
    param::NamedTuple = param,
    Ω::Dict{Int64, Dict{Int64, RandomVariables}} = Ω,
    binaryInfo::BinaryInfo = binaryInfo,
    stageDataList::Dict{Int, StageData} = stageDataList,
    forwardInfoList::Dict{Int, ForwardModelInfo} = forwardInfoList
)::Dict
    solCollection = Dict();  # to store every iteration results
    Ŝ = [0.0 for i in 1:binaryInfo.d];  
    for t in 1:param.T
        forwardInfo = forwardInfoList[t];
        ## realization of k-th scenario at stage t
        ω = Scenarios[k][t];
        ## the following function is used to (1). change the problem coefficients for different node within the same stage t.
        forward_modify_constraints!(
            forwardInfo,     
            stageDataList[t], 
            Ω[t][ω].d, 
            Ŝ, 
            binaryInfo = binaryInfo                                
        );
        optimize!(forwardInfo.model);

        solCollection[t, k] = ( 
            stageSolution = round.(JuMP.value.(forwardInfo.St), digits = 3), 
            stageValue = JuMP.objective_value(forwardInfo.model) - JuMP.value(forwardInfo.θ),
            OPT = JuMP.objective_value(forwardInfo.model)    
        );

        Ŝ  = solCollection[t, k].stageSolution;
    end
    return solCollection
end