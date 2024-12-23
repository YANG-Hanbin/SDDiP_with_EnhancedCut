#############################################################################################
###################################  function: backward pass ################################
#############################################################################################

"""
    This is the oracle in level set method, and it will return [F, dF]
"""

function backwardModel!(scenarioData::ScenarioData, x̂::Any;
                        prob::Prob = prob, 
                        indexSets::IndexSets = indexSets,
                        stageData::StageData = stageData,
                        timelimit::Int64 = 10, mipGap::Real = 1e-4,  
                        tightness::Bool = true
                        )

    ## construct forward second-stage problem 
    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => 0, 
                                          "Threads" => 0,
                                          "MIPGap" => mipGap, 
                                          "TimeLimit" => timelimit) 
                                          )

    @variable(model, π[indexSets.N] ≥ 0)  # gives the probability that the evader can travel from i to t undetected
    if tightness
        @variable(model, x[indexSets.D], Bin)  # local copy variable
        @constraint(model, sum(stageData.c[(i,j)] * x[(i,j)] for (i,j) in indexSets.D) ≤ stageData.b)
    else
        @variable(model, 0 ≤ x[indexSets.D] ≤ 1)  # local copy variable
    end
    @constraint(model, destinationConstraint, π[scenarioData.t] == 1)
    @constraint(model, [(i,j) in indexSets.D ], π[i] - prob.q[(i,j)] * π[j] ≥ 0)
    @constraint(model, [(i,j) in indexSets.Dᶜ], π[i] - prob.r[(i,j)] * π[j] ≥ 0)
    @constraint(model, SensorConstraint[(i,j) in indexSets.D ], π[i] - prob.r[(i,j)] * π[j] ≥ 
                                                    (prob.q[(i,j)] - prob.r[(i,j)]) * scenarioData.ϕ[j] * x[(i,j)]
                    )
    @constraint(model, nonantipativity[i in 1:length(keys(indexSets.D))], x[indexSets.D[keys(indexSets.D)[i]]] .== x̂[i])


    backward_stage2_Info = Backward_stage2_Info(model, x, π)
    return backward_stage2_Info
end


function backward_modify_constraints!(backward_stage2_Info::Backward_stage2_Info, 
                                        scenarioData::ScenarioData, 
                                        x̂::Vector{Float64}; 
                                        indexSets::IndexSets = indexSets,
                                        prob::Prob = prob                
                                        )

    delete(backward_stage2_Info.model, backward_stage2_Info.model[:destinationConstraint])
    for i in 1:length(keys(indexSets.D))
        delete(backward_stage2_Info.model, backward_stage2_Info.model[:SensorConstraint][indexSets.D[keys(indexSets.D)[i]]])
        delete(backward_stage2_Info.model, backward_stage2_Info.model[:nonantipativity][i])
    end

    unregister(backward_stage2_Info.model, :destinationConstraint)
    unregister(backward_stage2_Info.model, :SensorConstraint)
    unregister(backward_stage2_Info.model, :nonantipativity)

    ## no more than max num of generators
    @constraint(backward_stage2_Info.model, destinationConstraint,  backward_stage2_Info.π[scenarioData.t] == 1) 
    @constraint(backward_stage2_Info.model, SensorConstraint[(i,j) in indexSets.D ], backward_stage2_Info.π[i] - prob.r[(i,j)] * backward_stage2_Info.π[j] ≥ 
                                                                (prob.q[(i,j)] - prob.r[(i,j)]) * scenarioData.ϕ[j] * backward_stage2_Info.model[:x][(i,j)]
                                        )

    @constraint(backward_stage2_Info.model, nonantipativity[i in 1:length(keys(indexSets.D))], backward_stage2_Info.model[:x][indexSets.D[keys(indexSets.D)[i]]] .== x̂[i])


    @objective(backward_stage2_Info.model, Min, backward_stage2_Info.π[scenarioData.s])
end


