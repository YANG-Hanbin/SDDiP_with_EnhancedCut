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
    Q = Model(optimizer_with_attributes(
        () -> Gurobi.Optimizer(GRB_ENV), 
        "Threads" => 0)); 
    MOI.set(Q, MOI.Silent(), !param.verbose);
    set_optimizer_attribute(Q, "MIPGap", param.MIPGap);
    set_optimizer_attribute(Q, "TimeLimit", param.TimeLimit);
    # set_optimizer_attribute(Q, "MIPFocus", param.MIPFocus);           
    # set_optimizer_attribute(Q, "FeasibilityTol", param.FeasibilityTol);
    
    @variable(Q, x[i = 1:binaryInfo.d] ≥ 0, Int)        ## for current state, x is the number of generators will be built in this stage
    @variable(Q, y[i = 1:binaryInfo.d] ≥ 0)             ## amount of electricity
    @variable(Q, St[i = 1:binaryInfo.d] ≥ 0, Int)       ## stage variable, A * Lt is total number of generators built after this stage
    @variable(Q, slack ≥ 0 )
    @variable(Q, θ ≥ θ_bound)

    @objective(Q, Min, stageData.c1'* x + stageData.c2' * y + stageData.penalty * slack + θ )

    ## no more than max num of generators
    @constraint(Q, limitationConstraint,    St .≤ stageData.ū ) 

    # satisfy demand
    @constraint(Q, demandConstraint,        sum(y) + slack .≥ 0 )

    # no more than capacity
    @constraint(Q, capacityConstraint,      stageData.h * stageData.N 
                                                            * (St + stageData.s₀ ) .≥ y )  

    @constraint(Q, totalGenerators,         0.0 .+ x .== St)
    

    forwardModelInfo = ForwardModelInfo(Q, x, St, y, θ, slack)

    return forwardModelInfo
    # return [round.(JuMP.value.(Lt)), JuMP.value.(y), JuMP.value(θ), JuMP.objective_value(Q) - JuMP.value(θ)]  ## returen [Lt, y, θ, f]
end


"""
    forward_modify_constraints!(forwardInfo::ForwardModelInfo, 
                                        stageData::StageData, 
                                        demand::Vector{Float64}, 
                                        L̂::Vector{Float64};
                                        binaryInfo::BinaryInfo = binaryInfo                     ## realization of the random time
                                        )

TBW
"""
function forward_modify_constraints!(forwardInfo::ForwardModelInfo, 
                                        stageData::StageData, 
                                        demand::Vector{Float64}, 
                                        Ŝ::Vector{Float64};
                                        binaryInfo::BinaryInfo = binaryInfo                     ## realization of the random time
                                        )

    # delete(forwardInfo.model, forwardInfo.model[:limitationConstraint])
    delete(forwardInfo.model, forwardInfo.model[:demandConstraint])
    # delete(forwardInfo.model, forwardInfo.model[:capacityConstraint])
    delete(forwardInfo.model, forwardInfo.model[:totalGenerators])

    # unregister(forwardInfo.model, :limitationConstraint)
    unregister(forwardInfo.model, :demandConstraint)
    # unregister(forwardInfo.model, :capacityConstraint)
    unregister(forwardInfo.model, :totalGenerators)

    # ## no more than max num of generators
    # @constraint(forwardInfo.model, limitationConstraint,    Ŝ + forwardInfo.x .≤ stageData.ū ) 

    # satisfy demand
    @constraint(forwardInfo.model, demandConstraint,        sum(forwardInfo.y) + forwardInfo.slack .≥ demand )

    # # no more than capacity
    # @constraint(forwardInfo.model, capacityConstraint,      stageData.h * stageData.N 
    #                                                         * (Ŝ + forwardInfo.x + stageData.s₀ ) .≥ forwardInfo.y )  

    @constraint(forwardInfo.model, totalGenerators,         Ŝ + forwardInfo.x .== forwardInfo.St)

end
