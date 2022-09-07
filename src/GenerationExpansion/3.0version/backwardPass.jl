#############################################################################################
##########################    auxiliary functions for backward   ############################
#############################################################################################

function add_generator_constraint(stageData::StageData, modelInfo::BackwardModelInfo;
                                                        binaryInfo::BinaryInfo = binaryInfo)

    @constraint(modelInfo.model, binaryInfo.A * modelInfo.Lc + modelInfo.x .≤ stageData.ū )  ## no more than max num of generators
    @constraint(modelInfo.model, sum(modelInfo.y) + modelInfo.slack .≥ modelInfo.demand )  # satisfy demand
    @constraint(modelInfo.model, stageData.h * stageData.N 
                            * (binaryInfo.A * modelInfo.Lc + modelInfo.x + stageData.s₀ ) .≥ modelInfo.y )  # no more than capacity

end




function add_generator_cut(cutCoefficient::CutCoefficient, modelInfo::BackwardModelInfo)


    iter = length(keys(cutCoefficient.v))  ## iter num
    k = length(keys(cutCoefficient.v[1]))  ## scenario num

    @constraint(modelInfo.model, cut[i in 1:iter-1, m in 1:k], modelInfo.θ ≥ cutCoefficient.v[i][m] + 
                                                cutCoefficient.π[i][m]' * modelInfo.Lt )
                                                
end






#############################################################################################
###################################  function: backward pass ################################
#############################################################################################

"""
    This is the oracle in level set method, and it will return [F, dF]
"""
function backwardModel!(             stageData::StageData; 
                                            θ_bound::Real = 0.0,
                                            binaryInfo::BinaryInfo = binaryInfo, 
                                            timelimit::Int64 = 1)

    (A, n, d) = (binaryInfo.A, binaryInfo.n, binaryInfo.d)

    F = Model(
        optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                    "OutputFlag" => 0,
                                    "Threads" => 0,
                                    "MIPGap" => 1e-2, 
                                    "TimeLimit" => timelimit )
            )

    @variable(F, x[i = 1:binaryInfo.d] ≥ 0, Int)                   # the number of generators will be built in this stage
    # @variable(F, 0 ≤ Lc[i = 1:binaryInfo.n]≤ 1)                    # auxiliary variable (copy variable)
    @variable(F, Lc[i = 1:binaryInfo.n], Bin)                    # auxiliary variable (copy variable)
    @variable(F, Lt[i = 1:binaryInfo.n], Bin)                      # stage variable, A * Lt is total number of generators built after this stage
    @variable(F, y[i = 1:binaryInfo.d] ≥ 0)
    @variable(F, slack ≥ 0 )
    @variable(F, θ ≥ θ_bound)

    @constraint(F, binaryInfo.A * Lc + x .≤ stageData.ū )                       ## no more than max num of generators
    @constraint(F, stageData.h * stageData.N 
                * (binaryInfo.A * Lc + x + stageData.s₀ ) .≥ y )                ## no more than capacity

    @constraint(F, demandConstraint, sum(y) + slack .≥ 0)

    @constraint(F, hierarchicalConstriant, stageData.h * stageData.N 
                            * (binaryInfo.A * Lc + stageData.s₀ ) .≥ y)         ## constraints from last-stage: to restrict the region of local copy
    @constraint(F, binaryInfo.A * Lc + x .== binaryInfo.A * Lt )                ## to ensure pass a binary variable for next stage
    

    return BackwardModelInfo(F, x, Lt, Lc, y, θ, slack)
end


function backward_Constraint_Modification!(backwardInfo::BackwardModelInfo, 
                                        demand::Vector{Float64}                 
                                        )

    delete(backwardInfo.model, backwardInfo.model[:demandConstraint])
    unregister(backwardInfo.model, :demandConstraint)

    # satisfy demand
    @constraint(backwardInfo.model, demandConstraint, sum(backwardInfo.y) + backwardInfo.slack .≥ demand )
end

function copyVariable_Constraint!(backwardInfo::BackwardModelInfo, 
                                        stageData::StageData, 
                                        pdemand::Vector{Float64}; ## last-stage demand!
                                        binaryInfo::BinaryInfo = binaryInfo                  
                                        )

    delete(backwardInfo.model, backwardInfo.model[:hierarchicalConstriant])
    unregister(backwardInfo.model, :hierarchicalConstriant)

    # satisfy demand

    @constraint(backwardInfo.model, hierarchicalConstriant, stageData.h * stageData.N 
                                                                * (binaryInfo.A * backwardInfo.Lc + stageData.s₀ ) .≥ pdemand)     
end


