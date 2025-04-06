#############################################################################################
###################################  function: forward pass #################################
#############################################################################################

"""     
    function forwardModel!(
        stageData::StageData;
        θ_bound::Real = 0.0, 
        binaryInfo::BinaryInfo = binaryInfo, 
        timelimit::Int64 = 10, 
        mipGap::Real = 1e-4
    )::ForwardModelInfo
        1. construct the forward model
        2. return the forward model
        3. return the variables in the forward model
        4. return the objective function value of the forward model
"""
function forwardModel!(
    stageData::StageData;
    θ_bound::Real = 0.0, 
    binaryInfo::BinaryInfo = binaryInfo, 
    timelimit::Int64 = 10, 
    mipGap::Real = 1e-4
)::ForwardModelInfo
                            
    ## construct forward problem (3.1)
    model = Model(optimizer_with_attributes(
        ()->Gurobi.Optimizer(GRB_ENV), 
        "OutputFlag" => 0, 
        "Threads" => 0)
    ); 
    MOI.set(model, MOI.Silent(), true);
    set_optimizer_attribute(model, "MIPGap", mipGap);
    set_optimizer_attribute(model, "TimeLimit", timelimit);
    
    @variable(model, x[g = 1:binaryInfo.d] ≥ 0, Int)        ## for current state, x is the number of generators will be built in this stage
    @variable(model, y[g = 1:binaryInfo.d] ≥ 0)             ## amount of electricity
    @variable(model, St[g = 1:binaryInfo.d] ≥ 0, Int)       ## stage variable, A * Lt is total number of generators built after this stage
    @variable(model, slack ≥ 0 )
    @variable(model, θ ≥ θ_bound)

    sur = Dict(
        (g, i) => @variable(model, base_name = "sur[$g, $i]", binary = true)
        for g in 1:binaryInfo.d for i in 1:1
    );
    model[:sur] = sur;
    # constraints for surrogate variables
    ## Choosing one leaf node
    @constraint(
        model, 
        [g in 1:binaryInfo.d], 
        sur[g, 1] == 1
    );


    @objective(
        model, 
        Min, 
        stageData.c1'* x + stageData.c2' * y + stageData.penalty * slack + θ 
    );

    ## no more than max num of generators
    @constraint(model, limitationConstraint,    St .≤ stageData.ū ) 

    # satisfy demand
    @constraint(model, demandConstraint,        sum(y) + slack .≥ 0 )

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
        binaryInfo::BinaryInfo = binaryInfo                     ## realization of the random time
    )::Nothing
    1. modify the constraints in the forward model
    2. delete the constraints in the forward model
    3. add the new constraints in the forward model
    4. return nothing
"""
function forward_modify_constraints!(
    forwardInfo::ForwardModelInfo, 
    stageData::StageData, 
    demand::Vector{Float64}, 
    Ŝ::Vector{Float64};
    binaryInfo::BinaryInfo = binaryInfo                     ## realization of the random time
)::Nothing

    # delete(forwardInfo.model, forwardInfo.model[:limitationConstraint])
    delete(forwardInfo.model, forwardInfo.model[:demandConstraint])
    # delete(forwardInfo.model, forwardInfo.model[:capacityConstraint])
    delete(forwardInfo.model, forwardInfo.model[:totalGenerators])

    # unregister(forwardInfo.model, :limitationConstraint)
    unregister(forwardInfo.model, :demandConstraint)
    # unregister(forwardInfo.model, :capacityConstraint)
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
)
    solCollection = Dict();  # to store every iteration results
    Ŝ = [0.0 for i in 1:binaryInfo.d];  ## for the first-stage subproblem, we create a zero vector as 'x_ancestor'
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
            stageSur = Dict{Int64, Dict{Int64, Float64}}(
                g => Dict(i => round.(JuMP.value(forwardInfo.model[:sur][g, i]), digits = 3) for i in StateVarList[t].leaf[g]) 
                for g in 1:binaryInfo.d
                ),
            stageValue = JuMP.objective_value(forwardInfo.model) - JuMP.value(forwardInfo.θ),
            OPT = JuMP.objective_value(forwardInfo.model)
        );
        Ŝ  = solCollection[t, k].stageSolution;
    end
    return solCollection
end