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
        @constraint(model, sum(stageData.c[(i,j)] * x[(i,j)] for (i,j) in indexSet.D) ≤ stageData.b)
    else
        @variable(model, 0 ≤ x[indexSets.D] ≤ 1)  # local copy variable
    end
    @constraint(model, destinationConstraint, π[scenarioData.t] == 1)
    @constraint(model, [(i,j) in indexSets.D ], π[i] - prob.q[(i,j)] * π[j] ≥ 0)
    @constraint(model, [(i,j) in indexSets.Dᶜ], π[i] - prob.r[(i,j)] * π[j] ≥ 0)
    @constraint(model, SensorConstraint, [(i,j) in indexSets.D ], π[i] - prob.r[(i,j)] * π[j] ≥ 
                                                    (prob.q[(i,j)] - prob.r[(i,j)]) * scenarioData.ϕ[j] * x[(i,j)]
                    )
    @constraint(model, nonantipativity, x .== x̂)

    backward_stage2_Info = Backward_stage2_Info(model, x, π)
    return backward_stage2_Info
end


function backward_modify_constraints!(backward_stage2_Info::Backward_stage2_Info, 
                                        scenarioData::ScenarioData, 
                                        x̂::Any; 
                                        prob::Prob = prob                 
                                        )

    delete(backward_stage2_Info.model, backward_stage2_Info.model[:destinationConstraint])
    delete(backward_stage2_Info.model, backward_stage2_Info.model[:SensorConstraint])
    delete(backward_stage2_Info.model, backward_stage2_Info.model[:nonantipativity])

    unregister(backward_stage2_Info.model, :destinationConstraint)
    unregister(backward_stage2_Info.model, :SensorConstraint)
    unregister(backward_stage2_Info.model, :nonantipativity)

    ## no more than max num of generators
    @constraint(backward_stage2_Info.model, destinationConstraint,  backward_stage2_Info.π[scenarioData.t] == 1) 
    @constraint(backward_stage2_Info.model, SensorConstraint,  [(i,j) in indexSets.D ], backward_stage2_Info.π[i] - prob.r[(i,j)] * π[j] ≥ 
                                                                (prob.q[(i,j)] - prob.r[(i,j)]) * scenarioData.ϕ[j] * backward_stage2_Info.x[(i,j)]
                                        )
    @constraint(model, nonantipativity, backward_stage2_Info.x .== x̂)

    @objective(model, Min, backward_stage2_Info.π[scenarioData.s])
end

