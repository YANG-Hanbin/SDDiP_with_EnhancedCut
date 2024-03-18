"""
forwardModel!(; indexSets::IndexSets = indexSets, 
                            paramDemand::ParamDemand = paramDemand, 
                                paramOPF::ParamOPF = paramOPF, 
                                    stageRealization::StageRealization = stageRealization,
                                    θ_bound::Real = 0.0, outputFlag::Int64 = 0imelimit::Real = 100, mipGap::Float64 = 1e-2
                            )

# Arguments

    1. `indexSets::IndexSets` : index sets for the power network
    2. `paramDemand::ParamDemand` : demand parameters
    3. `paramOPF::ParamOPF` : OPF parameters
    4. `stageRealization::StageRealization` : realization of the stage
  
# Returns
    1. `model` : a forward pass model of stage t
  
"""
function forwardModel!(; indexSets::IndexSets = indexSets, 
                            paramDemand::ParamDemand = paramDemand, 
                                paramOPF::ParamOPF = paramOPF, 
                                    stageRealization::StageRealization = stageRealization,
                                    θ_bound::Real = 0.0, outputFlag::Int64 = 0, timelimit::Real = 3, mipGap::Float64 = 1e-3
                            )
    (D, G, L, B) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T) 
    (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L) 
    N = keys(stageRealization.prob)

    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                            "OutputFlag" => outputFlag, 
                                            "Threads" => 0, 
                                            "MIPGap" => mipGap, 
                                            "TimeLimit" => timelimit)
                                ) 
    @variable(model, θ_angle[B])      ## phase angle of the bus i
    @variable(model, P[L])            ## real power flow on line l; elements in L is Tuple (i, j)
    @variable(model, s[G])            ## real power generation at generator g
    @variable(model, 0 ≤ x[D] ≤ 1)    ## load shedding

    @variable(model, y[G], Bin)                 ## binary variable for generator commitment status
    @variable(model, v[G], Bin)                 ## binary variable for generator startup decision
    @variable(model, w[G], Bin)                 ## binary variable for generator shutdowm decision

    @variable(model, θ[N] ≥ θ_bound)            ## auxiliary variable for approximation of the value function

    # power flow constraints
    for l in L
        i = l[1]
        j = l[2]
        @constraint(model, P[l] ≤ - paramOPF.b[l] * (θ_angle[i] - θ_angle[j]))
        @constraint(model, P[l] ≥ - paramOPF.b[l] * (θ_angle[i] - θ_angle[j]))
    end
    
    # power flow limitation
    @constraint(model, [l in L], P[l] ≥ - paramOPF.W[l])
    @constraint(model, [l in L], P[l] ≤   paramOPF.W[l])
    # genertor limitation
    @constraint(model, [g in G], s[g] ≥ paramOPF.smin[g] * y[g])
    @constraint(model, [g in G], s[g] ≤ paramOPF.smax[g] * y[g])

    # power balance constriant
    @constraint(model, PowerBalance[i in B], sum(s[g] for g in Gᵢ[i]) -
                                                sum(P[(i, j)] for j in out_L[i]) + 
                                                    sum(P[(j, i)] for j in in_L[i]) 
                                                        .== sum(paramDemand.demand[d] * x[d] for d in Dᵢ[i]) )
    
    # on/off status with startup and shutdown decision
    @constraint(model, ShutUpDown[g in indexSets.G], v[g] - w[g] == model[:y][g])
    @constraint(model, Ramping1[g in indexSets.G], s[g] <= paramOPF.smin[g] * v[g])
    @constraint(model, Ramping2[g in indexSets.G], s[g] >= - paramOPF.M[g] * y[g] - paramOPF.smin[g] * w[g])

    # objective function
    @objective(model, Min, sum(paramOPF.slope[g] * s[g] +
                                paramOPF.intercept[g] * y[g] +
                                    paramOPF.C_start[g] * v[g] + 
                                        paramOPF.C_down[g] * w[g] for g in G) + 
                                            sum(paramDemand.w[d] * (1 - x[d]) for d in D) + sum(θ)
                )
    # @objective(model, Min, sum(paramOPF.slope[g] * s[g] for g in G) + 
    #                             sum(paramOPF.intercept[g] * y[g] for g in G) + 
    #                                 sum(paramOPF.C_start[g] * v[g] for g in G) + 
    #                                     sum(paramOPF.C_down[g] * w[g] for g in G) + 
    #                                         sum(paramDemand.w[d] * (1 - x[d]) for d in D) + sum(θ[n] for n in N)
    #             )
    return model
end

"""
forwardModification!(; model::Model = model)

# Arguments

    1. `model::Model` : a forward pass model of stage t
  
# Modification
    1. Remove the other scenario's demand balance constraints
    2. Add the current scenario's demand balance constriants
    3. Update its last stage decision with
"""

function forwardModification!(; model::Model = model, 
                            randomVariables::RandomVariables = randomVariables,
                                    paramOPF::ParamOPF = paramOPF, 
                                        stageDecision::Dict{Symbol, Dict{Int64, Float64}} = stageDecision, 
                                            indexSets::IndexSets = indexSets
                                        )

    for g in indexSets.G
        delete(model, model[:ShutUpDown][g])
        delete(model, model[:Ramping1][g])
        delete(model, model[:Ramping2][g])
    end
    unregister(model, :ShutUpDown)
    unregister(model, :Ramping1)
    unregister(model, :Ramping2)

    # on/off status with startup and shutdown decision
    @constraint(model, ShutUpDown[g in indexSets.G], model[:v][g] - model[:w][g] == model[:y][g] - stageDecision[:y][g])
    @constraint(model, Ramping1[g in indexSets.G], model[:s][g] - stageDecision[:s][g] <= paramOPF.M[g] * stageDecision[:y][g] + paramOPF.smin[g] * model[:v][g])
    @constraint(model, Ramping2[g in indexSets.G], model[:s][g] - stageDecision[:s][g] >= - paramOPF.M[g] * model[:y][g] - paramOPF.smin[g] * model[:w][g])

    # power balance constriant
    for i in indexSets.B
        delete(model, model[:PowerBalance][i])
    end
    unregister(model, :PowerBalance)
    @constraint(model, PowerBalance[i in indexSets.B], sum(model[:s][g] for g in indexSets.Gᵢ[i]) -
                                                            sum(model[:P][(i, j)] for j in indexSets.out_L[i]) + 
                                                                sum(model[:P][(j, i)] for j in indexSets.in_L[i]) 
                                                                    .== sum(paramDemand.demand[d] * randomVariables.deviation[d] * model[:x][d] for d in indexSets.Dᵢ[i]) )
end