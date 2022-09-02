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
                                            binaryInfo::BinaryInfo = binaryInfo )

    (A, n, d) = (binaryInfo.A, binaryInfo.n, binaryInfo.d)

    F = Model(
        optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                    "OutputFlag" => 0,
                                    "Threads" => 0 )
            )

    @variable(F, x[i = 1:binaryInfo.d] ≥ 0, Int)                   # the number of generators will be built in this stage
    @variable(F, 0 ≤ Lc[i = 1:binaryInfo.n]≤ 1)                    # auxiliary variable (copy variable)
    @variable(F, Lt[i = 1:binaryInfo.n], Bin)                      # stage variable, A * Lt is total number of generators built after this stage
    @variable(F, y[i = 1:binaryInfo.d] ≥ 0)
    @variable(F, slack ≥ 0 )
    @variable(F, θ ≥ θ_bound)

    @constraint(F, binaryInfo.A * Lc + x .≤ stageData.ū )           ## no more than max num of generators
    # @constraint(F, sum(y) + slack .≥ demand )                     # satisfy demand
    @constraint(F, stageData.h * stageData.N 
                * (binaryInfo.A * Lc + x + stageData.s₀ ) .≥ y )    # no more than capacity

    @constraint(F, demandConstraint, sum(y) + slack .≥ 0)
    @constraint(F, binarizationConstraint, x .== binaryInfo.A * Lt )               ## to ensure pass a binary variable for next stage

    return BackwardModelInfo(F, x, Lt, Lc, y, θ, slack)
end


function backward_Constraint_Modification!(backwardInfo::BackwardModelInfo, 
                                        stageData::StageData, 
                                        demand::Vector{Float64}, 
                                        L̂::Vector{Float64};
                                        binaryInfo::BinaryInfo = binaryInfo                  
                                        )

    delete(backwardInfo.model, backwardInfo.model[:demandConstraint])
    unregister(backwardInfo.model, :demandConstraint)

    delete(backwardInfo.model, backwardInfo.model[:binarizationConstraint])
    unregister(backwardInfo.model, :binarizationConstraint)

    # satisfy demand
    @constraint(backwardInfo.model, demandConstraint, sum(backwardInfo.y) + backwardInfo.slack .≥ demand )
    @constraint(backwardInfo.model, binarizationConstraint, binaryInfo.A * L̂ + backwardInfo.x .== binaryInfo.A * backwardInfo.Lt )               ## to ensure pass a binary variable for next stage

end


