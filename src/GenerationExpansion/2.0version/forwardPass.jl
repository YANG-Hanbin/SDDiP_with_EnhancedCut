#############################################################################################
###################################  function: forward pass #################################
#############################################################################################

"""     forward_step_optimize(StageCoefficient[t], b, x_ancestor, cut_collection, t)
        1. Solve the forward problem and return the useful info

        cutCoefficient: is the cut info given stage t
        stageData: is the param info given stage t
"""


function forwardModel!(stageData::StageData;
                                θ_bound::Real = 0.0, 
                                binaryInfo::BinaryInfo = binaryInfo )
                            
    ## construct forward problem (3.1)
    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => 0, 
                                          "Threads" => 0) 
                                          )
    @variable(Q, x[i = 1:binaryInfo.d] ≥ 0, Int)    ## for current state, x is the number of generators will be built in this stage
    @variable(Q, y[i = 1:binaryInfo.d] ≥ 0)         ## amount of electricity
    @variable(Q, Lt[i = 1:binaryInfo.n], Bin)       ## stage variable, A * Lt is total number of generators built after this stage
    @variable(Q, slack ≥ 0 )
    @variable(Q, θ ≥ θ_bound)

    @objective(Q, Min, stageData.c1'* x + stageData.c2' * y + stageData.penalty * slack + θ )

    ## no more than max num of generators
    @constraint(Q, limitationConstraint,    0.0 .+ x .≤ stageData.ū ) 

    # satisfy demand
    @constraint(Q, demandConstraint,        sum(y) + slack .≥ 0 )

    # no more than capacity
    @constraint(Q, capacityConstraint,      stageData.h * stageData.N 
                                                            * (0.0 .+ x + stageData.s₀ ) .≥ y )  

    @constraint(Q, totalGenerators,         0.0 .+ x .== binaryInfo.A * Lt)
    

    forwardModelInfo = ForwardModelInfo(Q, x, Lt, y, θ, slack)

    return forwardModelInfo
    # return [round.(JuMP.value.(Lt)), JuMP.value.(y), JuMP.value(θ), JuMP.objective_value(Q) - JuMP.value(θ)]  ## returen [Lt, y, θ, f]
end


function forward_modify_constraints!(forwardInfo::ForwardModelInfo, 
                                        stageData::StageData, 
                                        demand::Vector{Float64}, 
                                        sum_generator::Vector{Float64};
                                        binaryInfo::BinaryInfo = binaryInfo                     ## realization of the random time
                                        )

    delete(forwardInfo.model, forwardInfo.model[:limitationConstraint])
    delete(forwardInfo.model, forwardInfo.model[:demandConstraint])
    delete(forwardInfo.model, forwardInfo.model[:capacityConstraint])
    delete(forwardInfo.model, forwardInfo.model[:totalGenerators])

    unregister(forwardInfo.model, :limitationConstraint)
    unregister(forwardInfo.model, :demandConstraint)
    unregister(forwardInfo.model, :capacityConstraint)
    unregister(forwardInfo.model, :totalGenerators)

    ## no more than max num of generators
    @constraint(forwardInfo.model, limitationConstraint,    binaryInfo.A * sum_generator + forwardInfo.x .≤ stageData.ū ) 

    # satisfy demand
    @constraint(forwardInfo.model, demandConstraint,        sum(forwardInfo.y) + forwardInfo.slack .≥ demand )

    # no more than capacity
    @constraint(forwardInfo.model, capacityConstraint,      stageData.h * stageData.N 
                                                            * (binaryInfo.A * sum_generator + forwardInfo.x + stageData.s₀ ) .≥ forwardInfo.y )  

    @constraint(forwardInfo.model, totalGenerators,         binaryInfo.A * sum_generator + forwardInfo.x .== binaryInfo.A * forwardInfo.Lt)

end



