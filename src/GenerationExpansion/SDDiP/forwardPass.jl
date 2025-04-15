#############################################################################################
###################################  function: forward pass #################################
#############################################################################################
"""     
    forwardModel!(
        stageData::StageData;                            
        θ_bound::Real = 0.0, 
        binaryInfo::BinaryInfo = binaryInfo, 
        timelimit::Int64 = 10, 
        mipGap::Real = 1e-4
    )::ForwardModelInfo  
        
    1. Solve the forward problem and return the useful info
    2. cutCoefficient: is the cut info given stage t
    3. stageData: is the param info given stage t
"""
function forwardModel!(
    stageData::StageData;                            
    θ_bound::Real = 0.0, 
    binaryInfo::BinaryInfo = binaryInfo, 
    timelimit::Int64 = 10, 
    mipGap::Real = 1e-4
)::ForwardModelInfo                   
    ## construct forward problem (3.1)
    model = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => 0, 
            "Threads" => 0
        )
    ); 
    MOI.set(model, MOI.Silent(), true);
    set_optimizer_attribute(model, "MIPGap", mipGap);
    set_optimizer_attribute(model, "TimeLimit", timelimit);

    @variable(model, x[i = 1:binaryInfo.d] ≥ 0, Int)    ## for current state, x is the number of generators will be built in this stage
    @variable(model, y[i = 1:binaryInfo.d] ≥ 0)         ## amount of electricity
    @variable(model, Lt[i = 1:binaryInfo.n], Bin)       ## stage variable, A * Lt is total number of generators built after this stage
    @variable(model, slack ≥ 0 )
    @variable(model, θ ≥ θ_bound)


    @expression(
        model, 
        primal_objective_expression, 
        stageData.c1'* x + stageData.c2' * y + stageData.penalty * slack + θ
    );

    @objective(
        model, 
        Min, 
        primal_objective_expression
    );

    ## no more than max num of generators
    @constraint(
        model, 
        limitationConstraint, 
        0.0 .+ x .≤ stageData.ū
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
        stageData.h * stageData.N * (0.0 .+ x + stageData.s₀ ) .≥ y 
    ); 

    @constraint(
        model, 
        totalGenerators,         
        0.0 .+ x .== binaryInfo.A * Lt
    );

    # a cap constraint for two generators
    @constraint(
        model,
        [i = 4:5],
        y[i] ≤ sum(y)/5
    );
    return ForwardModelInfo(model, x, Lt, y, θ, slack)
end


"""
    function forward_modify_constraints!(
        forwardInfo::ForwardModelInfo, 
        stageData::StageData, 
        demand::Vector{Float64}, 
        L̂::Vector{Float64};
        binaryInfo::BinaryInfo = binaryInfo                     ## realization of the random time
    )::Nothing

    1. modify the constraints of forward model
    2. update the demand and capacity constraints
    3. update the total number of generators built after this stage
    4. update the slack variable
    5. update the objective function
"""
function forward_modify_constraints!(
    forwardInfo::ForwardModelInfo, 
    stageData::StageData, 
    demand::Vector{Float64}, 
    L̂::Vector{Float64};
    binaryInfo::BinaryInfo = binaryInfo                     ## realization of the random time
)::Nothing

    delete(forwardInfo.model, forwardInfo.model[:limitationConstraint])
    delete(forwardInfo.model, forwardInfo.model[:demandConstraint])
    delete(forwardInfo.model, forwardInfo.model[:capacityConstraint])
    delete(forwardInfo.model, forwardInfo.model[:totalGenerators])

    unregister(forwardInfo.model, :limitationConstraint)
    unregister(forwardInfo.model, :demandConstraint)
    unregister(forwardInfo.model, :capacityConstraint)
    unregister(forwardInfo.model, :totalGenerators)

    ## no more than max num of generators
    @constraint(
        forwardInfo.model, 
        limitationConstraint,    
        binaryInfo.A * L̂ + forwardInfo.x .≤ stageData.ū
    ); 

    # satisfy demand
    @constraint(
        forwardInfo.model, 
        demandConstraint,        
        sum(forwardInfo.y) + forwardInfo.slack .≥ demand 
    );

    # no more than capacity
    @constraint(
        forwardInfo.model, 
        capacityConstraint,      
        stageData.h * stageData.N * (binaryInfo.A * L̂ + forwardInfo.x + stageData.s₀ ) .≥ forwardInfo.y 
    );  

    @constraint(
        forwardInfo.model, 
        totalGenerators,         
        binaryInfo.A * L̂ + forwardInfo.x .== binaryInfo.A * forwardInfo.Lt
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
    L̂= [0.0 for i in 1:binaryInfo.n];  ## for the first-stage subproblem, we create a zero vector as 'x_ancestor'
    for t in 1:param.T
        forwardInfo = forwardInfoList[t];
        ## realization of k-th scenario at stage t
        ω = Scenarios[k][t];
        ## the following function is used to (1). change the problem coefficients for different node within the same stage t.
        forward_modify_constraints!(
            forwardInfo,           
            stageDataList[t],                  
            Ω[t][ω].d,                          
            L̂,      
            binaryInfo = binaryInfo                
        );
        optimize!(forwardInfo.model);

        solCollection[t, k] = ( 
            stageSolution = round.(JuMP.value.(forwardInfo.Lt), digits = 3), 
            stageValue = JuMP.objective_value(forwardInfo.model) - JuMP.value(forwardInfo.θ),                       
            OPT = JuMP.objective_value(forwardInfo.model)
        );

        L̂ = round.(JuMP.value.(forwardInfo.Lt), digits = 3);
    end
    return solCollection
end
