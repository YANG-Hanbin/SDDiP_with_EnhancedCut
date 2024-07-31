#############################################################################################
###################################  function: backward pass ################################
#############################################################################################

"""
    This is the oracle in level set method, and it will return [F, dF]
"""
function backwardModel!(             stageData::StageData; 
                                            θ_bound::Real = 0.0,
                                            binaryInfo::BinaryInfo = binaryInfo, 
                                            timelimit::Int64 = 3, mipGap::Real = 1e-4, tightness::Bool = false)

    F = Model(
        optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                    "OutputFlag" => 0,
                                    "Threads" => 0,
                                    "MIPGap" => mipGap, 
                                    "TimeLimit" => timelimit )
            )

    @variable(F, x[g = 1:binaryInfo.d] ≥ 0, Int)                   # the number of generators will be built in this stage
    @variable(F, y[g = 1:binaryInfo.d] ≥ 0)
    @variable(F, St[g = 1:binaryInfo.d] ≥ 0, Int)                      # stage variable, A * Lt is total number of generators built after this stage
    @variable(F, slack ≥ 0 )
    @variable(F, θ ≥ θ_bound)


    # @variable(model, sur[G, 1:max_sur], Bin)                        
    sur = Dict(
        (g, i) => @variable(F, base_name = "sur[$g, $i]", binary = true)
        for g in 1:binaryInfo.d for i in 1:1
    );
    F[:sur] = sur;                                              ## sur[g, k] is the kth surrogate variable of s[g]

    # constraints for surrogate variables
    ## Choosing one leaf node
    @constraint(F, [g in 1:binaryInfo.d], sur[g, 1] == 1)

    # auxiliary variable (copy variable) for S_{t-1} or Ŝ
    if tightness
        @variable(F, Sc[i = 1:binaryInfo.d] ≥ 0, Int)   
        @constraint(F, Sc .≤ stageData.ū)                
        sur_copy = Dict(
            (g, i) => @variable(F, base_name = "sur_copy[$g, $i]", binary = true)
            for g in 1:binaryInfo.d for i in 1:1
        );
        F[:sur_copy] = sur_copy;     
    else
        @variable(F, Sc[i = 1:binaryInfo.d] ≥ 0)   
        @constraint(F, Sc .≤ stageData.ū)                
        sur_copy = Dict(
            (g, i) => @variable(F, base_name = "sur_copy[$g, $i]", lower_bound = 0, upper_bound = 1)
            for g in 1:binaryInfo.d for i in 1:1
        );
        F[:sur_copy] = sur_copy;     
    end
    @constraint(F, [g in 1:binaryInfo.d], sur_copy[g, 1] == 1)


    @constraint(F, St .≤ stageData.ū )                       ## no more than max num of generators
    @constraint(F, stageData.h * stageData.N 
                * (St + stageData.s₀ ) .≥ y )                ## no more than capacity

    @constraint(F, demandConstraint, sum(y) + slack .≥ 0)

    # @constraint(F, hierarchicalConstriant, stageData.h * stageData.N 
    #                                         * (Sc + stageData.s₀ ) .≥ y)         ## constraints from last-stage: to restrict the region of local copy
    @constraint(F, Sc + x .== St )                ## to ensure pass a binary variable for next stage
    

    return BackwardModelInfo(F, x, St, Sc, y, θ, slack)
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
                                        pdemand::Vector{Float64}                 
                                        )

    delete(backwardInfo.model, backwardInfo.model[:hierarchicalConstriant])
    unregister(backwardInfo.model, :hierarchicalConstriant)

    # satisfy demand
    @constraint(backwardInfo.model, hierarchicalConstriant, stageData.h * stageData.N 
                                                                * (Sc + stageData.s₀ ) .≥ pdemand)     
end


