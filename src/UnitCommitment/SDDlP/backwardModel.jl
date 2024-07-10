"""
backwardModel!(; indexSets::IndexSets = indexSets, 
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
function backwardModel!(; tightness::Bool = true, indexSets::IndexSets = indexSets, 
                            paramDemand::ParamDemand = paramDemand, 
                                paramOPF::ParamOPF = paramOPF, 
                                    stageRealization::StageRealization = stageRealization,
                                        max_sur::Int64 = 100,
                                            θ_bound::Real = 0.0, outputFlag::Int64 = 0, timelimit::Real = 5, mipGap::Float64 = 1e-4 
                            )
    (D, G, L, B) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T) 
    (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L) 
    N = keys(stageRealization.prob)

    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                            "OutputFlag" => outputFlag, 
                                            "Threads" => 0, 
                                            "MIPGap" => mipGap, 
                                            "TimeLimit" => timelimit);
                                ) 
    @variable(model, θ_angle[B])                                    ## phase angle of the bus i
    @variable(model, P[L])                                          ## real power flow on line l; elements in L is Tuple (i, j)
    @variable(model, 0 ≤ s[g in G] ≤ paramOPF.smax[g])              ## real power generation at generator g
    @variable(model, 0 ≤ x[D] ≤ 1)                                  ## load shedding

    @variable(model, y[G], Bin)                                     ## binary variable for generator commitment status
    @variable(model, v[G], Bin)                                     ## binary variable for generator startup decision
    @variable(model, w[G], Bin)                                     ## binary variable for generator shutdowm decision

    @variable(model, θ[N] ≥ θ_bound)                                ## auxiliary variable for approximation of the value function
    @variable(model, sur[G, 1:max_sur], Bin)                        ## sur[g, k] is the kth surrogate variable of s[g]

    # constraints for surrogate variables
    ## Choosing one leaf node
    @constraint(model, [g in G], sur[g, 1] == 1)

    # copy variables: :s, :y, :sur
    if tightness
        @variable(model, 0 ≤ s_copy[g in G] ≤ paramOPF.smax[g])
        @variable(model, y_copy[G], Bin)        
        @variable(model, sur_copy[G, 1:max_sur], Bin)                        ## sur[g, k] is the kth surrogate variable of s[g]
    else
        @variable(model, 0 ≤ s_copy[g in G] ≤ paramOPF.smax[g])
        @variable(model, 0 ≤ y_copy[G] ≤ 1)  
        @variable(model, 0 ≤ sur_copy[G, 1:max_sur] ≤ 1)                     ## sur[g, k] is the kth surrogate variable of s[g]
    end
    @constraint(model, [g in G], sur_copy[g, 1] == 1)

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
    @constraint(model, ShutUpDown[g in G], v[g] - w[g] == y[g] - y_copy[g])
    @constraint(model, Ramping1[g in G], s[g] - s_copy[g] <= paramOPF.M[g] * y_copy[g] + paramOPF.smin[g] * v[g])
    @constraint(model, Ramping2[g in G], s[g] - s_copy[g] >= - paramOPF.M[g] * y[g] - paramOPF.smin[g] * w[g])

    return model
end

"""
backwardModification!(; model::Model = model)

# Arguments

    1. `model::Model` : a backward pass model of stage t
  
# Modification
    1. Remove the other scenario's demand balance constraints
    2. Add the current scenario's demand balance constriants
"""
function backwardModification!(; model::Model = model, 
                            randomVariables::RandomVariables = randomVariables,
                                paramDemand::ParamDemand = paramDemand, 
                                    paramOPF::ParamOPF = paramOPF, 
                                        # stageDecision::Dict{Symbol, Dict{Int64, Any}} = stageDecision, 
                                            indexSets::IndexSets = indexSets
                                        )

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
"""
getValueFunctionLB(; model::Model = model, 
                            paramDemand::ParamDemand = paramDemand, 
                                paramOPF::ParamOPF = paramOPF, 
                                    stageDecision::StageDecision = stageDecision, 
                                        indexSets::IndexSets = indexSets)

# Arguments

    1. `model::Model` : a backward pass model of stage t
  
# Modification
    1. Add Nonanticipativity constraints and objective function
    2. Obtain the LB
    3. Remove the Nonanticipativity constraints
# Return 
    LB
"""
function getValueFunctionLB(; model::Model = model, 
                                        paramDemand::ParamDemand = paramDemand, 
                                            paramOPF::ParamOPF = paramOPF, 
                                                stageDecision::Dict{Symbol, Dict{Int64, Any}} = stageDecision, 
                                                    indexSets::IndexSets = indexSets
                                                )
    @constraint(model, Nonanticipativity[g in indexSets.G, k in keys(stageDecision[:sur][g])], model[:sur][g, k] == stageDecision[:sur][g][k]);
    @objective(model, Min,  sum(paramOPF.slope[g] * model[:s][g] + paramOPF.intercept[g] * model[:y][g] +
                                                                        paramOPF.C_start[g] * model[:v][g] + 
                                                                            paramOPF.C_down[g] * model[:w][g] for g in indexSets.G) + 
                                                                                sum(paramDemand.w[d] * (1 - model[:x][d]) for d in indexSets.D) + sum(model[:θ]));
    optimize!(model); f_star_value = JuMP.objective_value(model);
    for g in indexSets.G
        for k in keys(stageDecision[:sur][g])
            delete(model, model[:Nonanticipativity][g,k]);
        end
    end
    unregister(model, :Nonanticipativity);
    return f_star_value 
end