#############################################################################################
###################################  function: forward pass #################################
#############################################################################################

"""     forward_step_optimize(StageCoefficient[t], b, x_ancestor, cut_collection, t)
        1. Solve the forward problem and return the useful info

        cutCoefficient: is the cut info given stage t
        stageData: is the param info given stage t
"""


function forward_stage1_Model!(; stageData::StageData = stageData, indexSets::IndexSets = indexSets,
                                prob::Prob = prob,
                                timelimit::Int64 = 10, mipGap::Real = 1e-4 )
                            
    ## construct forward master problem
    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => 0, 
                                          "Threads" => 0,
                                          "MIPGap" => mipGap, 
                                          "TimeLimit" => timelimit) 
                                          )
    @variable(model, x[indexSets.D], Bin) # first-stage binary variables denote if a sensor is installed on arc (i,j) ∈ D
    @variable(model, θ[indexSets.Ω] ≥ 0)  # approximation of value functions for each scenario

    @constraint(model, [ω in indexSets.Ω], θ[ω] ≤ 1)
    @constraint(model, sum(stageData.c[(i,j)] * x[(i,j)] for (i,j) in indexSets.D) ≤ stageData.b )    

    @objective(model, Min, sum(prob.p[ω] * θ[ω] for ω in indexSets.Ω))

    forward_stage1_Info = Forward_stage1_Info(model, x, θ)

    return forward_stage1_Info
    # return [round.(JuMP.value.(Lt)), JuMP.value.(y), JuMP.value(θ), JuMP.objective_value(Q) - JuMP.value(θ)]  ## returen [Lt, y, θ, f]
end



"""     forward_step_optimize(StageCoefficient[t], b, x_ancestor, cut_collection, t)
        1. Solve the forward problem and return the useful info

        cutCoefficient: is the cut info given stage t
        stageData: is the param info given stage t
"""


function forward_stage2_Model!(scenarioData::ScenarioData, x̂::Any;
                                prob::Prob = prob, indexSets::IndexSets = indexSets,
                                timelimit::Int64 = 10, mipGap::Real = 1e-4 )
                            
    ## construct forward second-stage problem 
    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => 0, 
                                          "Threads" => 0,
                                          "MIPGap" => mipGap, 
                                          "TimeLimit" => timelimit) 
                                          )

    @variable(model, π[indexSets.N] ≥ 0)  # gives the probability that the evader can travel from i to t undetected
    @variable(model, 0 ≤ x[indexSets.D] ≤ 1)  # local copy variable

    @constraint(model, destinationConstraint, π[scenarioData.t] == 1)
    @constraint(model, [(i,j) in indexSets.D ], π[i] - prob.q[(i,j)] * π[j] ≥ 0)
    @constraint(model, [(i,j) in indexSets.Dᶜ], π[i] - prob.r[(i,j)] * π[j] ≥ 0)
    @constraint(model, SensorConstraint[(i,j) in indexSets.D ], π[i] - prob.r[(i,j)] * π[j] ≥ 
                                                    (prob.q[(i,j)] - prob.r[(i,j)]) * scenarioData.ϕ[j] * x[(i,j)]
                    )
    @constraint(model, nonantipativity[i in 1:length(keys(indexSets.D))], x[indexSets.D[keys(indexSets.D)[i]]] .== x̂[i])
    @objective(model, Min, π[scenarioData.s])

    forward_stage2_Info = Forward_stage2_Info(model, x, π)

    return forward_stage2_Info
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
function forward_modify_constraints!(forward_stage2_Info::Forward_stage2_Info, 
                                        scenarioData::ScenarioData, 
                                        x̂::Vector{Float64}; 
                                        prob::Prob = prob, indexSets::IndexSets = indexSets               
                                        )

    delete(forward_stage2_Info.model, forward_stage2_Info.model[:destinationConstraint])
    for i in 1:length(keys(indexSets.D))
        delete(forward_stage2_Info.model, forward_stage2_Info.model[:SensorConstraint][indexSets.D[keys(indexSets.D)[i]]])
        delete(forward_stage2_Info.model, forward_stage2_Info.model[:nonantipativity][i])
    end
    

    unregister(forward_stage2_Info.model, :destinationConstraint)
    unregister(forward_stage2_Info.model, :SensorConstraint)
    unregister(forward_stage2_Info.model, :nonantipativity)

    ## no more than max num of generators
    @constraint(forward_stage2_Info.model, destinationConstraint,  forward_stage2_Info.π[scenarioData.t] == 1) 
    @constraint(forward_stage2_Info.model, SensorConstraint[(i,j) in indexSets.D ], forward_stage2_Info.π[i] - prob.r[(i,j)] * forward_stage2_Info.π[j] ≥ 
                                                                (prob.q[(i,j)] - prob.r[(i,j)]) * scenarioData.ϕ[j] * forward_stage2_Info.model[:x][(i,j)]
                                        )

    @constraint(forward_stage2_Info.model, nonantipativity[i in 1:length(keys(indexSets.D))], forward_stage2_Info.model[:x][indexSets.D[keys(indexSets.D)[i]]] .== x̂[i])

    @objective(forward_stage2_Info.model, Min, forward_stage2_Info.π[scenarioData.s])
end
