#############################################################################################
##########################    auxiliary functions for forward    ############################
#############################################################################################

function add_generator_constraint(stageData::StageData, modelInfo::ForwardModelInfo;
                                  binaryInfo::BinaryInfo = binaryInfo)    

  @constraint(modelInfo.model, binaryInfo.A * modelInfo.sum_generator + modelInfo.x .≤ stageData.ū )  ## no more than max num of generators
  
  @constraint(modelInfo.model, sum(modelInfo.y) + modelInfo.slack .≥ modelInfo.demand )  # satisfy demand
  
  @constraint(modelInfo.model, stageData.h * stageData.N 
                               * (binaryInfo.A * modelInfo.sum_generator + modelInfo.x + stageData.s₀ ) .≥ modelInfo.y )  # no more than capacity
  
end



function add_generator_cut(cutCoefficient::CutCoefficient, modelInfo::ForwardModelInfo)

  iter = length(keys(cutCoefficient.v))  ## iter num
  k = length(keys(cutCoefficient.v[1]))  ## scenario num

  @constraint(modelInfo.model, cut[i in 1:iter-1, m in 1:k], modelInfo.θ ≥ cutCoefficient.v[i][m] + 
                                                                    cutCoefficient.π[i][m]' * modelInfo.Lt )

end






#############################################################################################
###################################  function: forward pass #################################
#############################################################################################

"""     forward_step_optimize(StageCoefficient[t], b, x_ancestor, cut_collection, t)
        1. Solve the forward problem and return the useful info

        cutCoefficient: is the cut info given stage t
        stageData: is the param info given stage t
"""


function forward_step_optimize!(stageData::StageData, 
                                demand::Vector{Float64}, 
                                sum_generator::Vector{Float64}, 
                                cutCoefficient::CutCoefficient;
                                θ_bound::Real = 0.0, binaryInfo::BinaryInfo = binaryInfo )

    (A, n, d) = (binaryInfo.A, binaryInfo.n, binaryInfo.d)                            
    ## construct forward problem (3.1)
    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => 0, 
                                          "Threads" => 0) 
                                          )
    @variable(Q, x[i = 1:d] ≥ 0, Int)   ## for current state, x is the number of generators will be built in this stage
    @variable(Q, y[i = 1:d] ≥ 0)        ## amount of electricity
    @variable(Q, Lt[i = 1:n], Bin)       ## stage variable, A * Lt is total number of generators built after this stage
    @variable(Q, slack ≥ 0 )
    @variable(Q, θ ≥ θ_bound)
    modelInfo = ForwardModelInfo(Q, x, Lt, y, θ, demand, slack, sum_generator)

    add_generator_constraint(stageData, modelInfo, binaryInfo = binaryInfo)
    add_generator_cut(cutCoefficient, modelInfo)
    @constraint(Q, A * sum_generator + x.== A * Lt)

    @objective(Q, Min, stageData.c1'* x + stageData.c2' * y + stageData.penalty * slack + θ )
    optimize!(Q)

    return [round.(JuMP.value.(Lt)), JuMP.value.(y), JuMP.value(θ), JuMP.objective_value(Q) - JuMP.value(θ)]  ## returen [Lt, y, θ, f]
end

