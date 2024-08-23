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
function backwardModel!(; tightness::Bool = tightness, indexSets::IndexSets = indexSets, 
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
    @variable(model, θ_angle[B])                                    ## phase angle of the bus i
    @variable(model, P[L])                                          ## real power flow on line l; elements in L is Tuple (i, j)
    @variable(model, 0 ≤ s[g in G] ≤ paramOPF.smax[g])              ## real power generation at generator g
    @variable(model, 0 ≤ x[D] ≤ 1)                                  ## load shedding

    @variable(model, y[G], Bin)                 ## binary variable for generator commitment status
    @variable(model, v[G], Bin)                 ## binary variable for generator startup decision
    @variable(model, w[G], Bin)                 ## binary variable for generator shutdowm decision

    @variable(model, θ[N] ≥ θ_bound)            ## auxiliary variable for approximation of the value function

    # copy variables: :s, :y
    if tightness
        # @variable(model, s_copy[G])
        @variable(model, 0 ≤ s_copy[g in G] ≤ paramOPF.smax[g])
        @variable(model, y_copy[G], Bin)        
    else
        @variable(model, s_copy[G])
        @variable(model, 0 ≤ y_copy[G] ≤ 1)       
    end

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

    1. `model::Model` : a forward pass model of stage t
  
# Modification
    1. Remove the other scenario's demand balance constraints
    2. Add the current scenario's demand balance constriants
"""
function backwardModification!(; model::Model = model, 
                            randomVariables::RandomVariables = randomVariables,
                                    paramOPF::ParamOPF = paramOPF, 
                                        paramDemand::ParamDemand = paramDemand,
                                        # stageDecision::StageDecision = stageDecision, 
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
backwardPass(backwardNodeInfo)

function for backward pass in parallel computing
"""
function backwardPass(backwardNodeInfo::Tuple; 
                            indexSets::IndexSets = indexSets, 
                            paramDemand::ParamDemand = paramDemand, 
                            paramOPF::ParamOPF = paramOPF, max_iter::Int64 = max_iter, Output_Gap::Bool = Output_Gap, tightness::Bool = tightness, δ::Float64 = δ,
                            backwardInfoList::Dict{Int64, Model} = backwardInfoList, forwardInfoList::Dict{Int64, Model} = forwardInfoList, scenarioTree::ScenarioTree = scenarioTree, solCollection::Dict{Any, Any} = solCollection
                            )
    (i, t, n, ω, cutSelection) = backwardNodeInfo; 

    # forwardModification!(model = forwardInfoList[t], randomVariables = scenarioTree.tree[t].nodes[n], stageDecision = solCollection[i, t-1, ω].stageSolution, paramOPF = paramOPF, indexSets = indexSets, paramDemand = paramDemand);
    # optimize!(forwardInfoList[t]); f_star_value = JuMP.objective_value(forwardInfoList[t]);

    backwardModification!(model = backwardInfoList[t], randomVariables = scenarioTree.tree[t].nodes[n], paramOPF = paramOPF, indexSets = indexSets, paramDemand = paramDemand);

    (x_interior, levelSetMethodParam, x₀) = setupLevelSetMethod(stageDecision = solCollection[i, t-1, ω].stageSolution, 
                                                        f_star_value = solCollection[i, t, ω].OPT, 
                                                            cutSelection = "LC", max_iter = max_iter,
                                                                Output_Gap = Output_Gap, ℓ = .0, λ = .3);
    # model = backwardInfoList[t]; stageDecision = solCollection[i, t-1, ω].stageSolution;
    ((λ₀, λ₁), LMiter) = LevelSetMethod_optimization!(levelSetMethodParam = levelSetMethodParam, model = backwardInfoList[t],
                                            cutSelection = "LC",
                                                stageDecision = solCollection[i, t-1, ω].stageSolution, tightness = tightness,
                                                        x_interior = x_interior, x₀ = x₀, indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF, δ = δ);

    
    if cutSelection != "LC"  # && gap ≥ 5e-2
        f_star_value = λ₀ + sum(λ₁[:s][g] * solCollection[i, t-1, ω].stageSolution[:s][g] + 
                                λ₁[:y][g] * solCollection[i, t-1, ω].stageSolution[:y][g] for g in indexSets.G );
        (x_interior, levelSetMethodParam, x₀) = setupLevelSetMethod(stageDecision = solCollection[i, t-1, ω].stageSolution, 
                                                                    f_star_value = f_star_value, 
                                                                    cutSelection = cutSelection, max_iter = max_iter,
                                                                    Output_Gap = Output_Gap, ℓ = .0, λ = .1 );
        # model = backwardInfoList[t]; stageDecision = solCollection[i, t-1, ω].stageSolution;
        ((λ₀, λ₁), LMiter) = LevelSetMethod_optimization!(levelSetMethodParam = levelSetMethodParam, model = backwardInfoList[t],
                                            cutSelection = cutSelection,
                                                stageDecision = solCollection[i, t-1, ω].stageSolution, tightness = tightness,
                                                        x_interior = x_interior, x₀ = x₀, indexSets = indexSets, paramDemand = paramDemand, paramOPF = paramOPF, δ = δ);
    end

    return ((λ₀, λ₁), LMiter)  
end