"""
    function forwardModel!( 
        paramDemand::ParamDemand, 
        paramOPF::ParamOPF, 
        stageRealization::StageRealization;
        indexSets::IndexSets = indexSets, 
        para::NamedTuple = para
    )
# Arguments

    1. `indexSets::IndexSets` : index sets for the power network
    2. `paramDemand::ParamDemand` : demand parameters
    3. `paramOPF::ParamOPF` : OPF parameters
    4. `stageRealization::StageRealization` : realization of the stage
  
# Returns
    1. `SDDPLModel`
"""
function forwardModel!( 
    paramDemand::ParamDemand, 
    paramOPF::ParamOPF, 
    stageRealization::StageRealization;
    indexSets::IndexSets = indexSets, 
    param::NamedTuple = param
)::SDDPLModel
    ## build the forward model
    model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), 
                                                    "Threads" => 0)); 
    MOI.set(model, MOI.Silent(), !param.verbose);
    set_optimizer_attribute(model, "MIPGap", param.MIPGap);
    set_optimizer_attribute(model, "TimeLimit", param.TimeLimit);
                
    ## define variables
    @variable(model, θ_angle[indexSets.B])                                                ## phase angle of the bus i
    @variable(model, P[indexSets.L])                                                      ## real power flow on line l; elements in L is Tuple (i, j)
    @variable(model, 0 ≤ s[g in indexSets.G] ≤ paramOPF.smax[g])                          ## real power generation at generator g
    @variable(model, 0 ≤ x[indexSets.D] ≤ 1)                                              ## load shedding

    @variable(model, y[indexSets.G], Bin)                                                 ## binary variable for generator commitment status
    @variable(model, v[indexSets.G], Bin)                                                 ## binary variable for generator startup decision
    @variable(model, w[indexSets.G], Bin)                                                 ## binary variable for generator shutdowm decision

    @variable(model, h[indexSets.G] ≥ 0);                                                 ## production cost at generator g

    @variable(model, θ[keys(stageRealization.prob)] ≥ param.θ̲)                             ## auxiliary variable for approximation of the value function

    ## augmented variables
    # define the augmented variables for cont. variables
    augmentVar = Dict(
        (g, k) => @variable(model, base_name = "augmentVar[$g, $k]", binary = true)
        for g in indexSets.G for k in 1:1
    );
    model[:augmentVar] = augmentVar;
    
    ## define copy variables
    if param.tightness
        @variable(model, 0 ≤ s_copy[g in indexSets.G] ≤ paramOPF.smax[g])
        @variable(model, y_copy[indexSets.G], Bin)        
        augmentVar_copy = Dict(
            (g, k) => @variable(model, base_name = "augmentVar_copy[$g, $k]", binary = true)
            for g in indexSets.G for k in 1:1
        );
        model[:augmentVar_copy] = augmentVar_copy;     
    else
        @variable(model, 0 ≤ s_copy[g in indexSets.G] ≤ paramOPF.smax[g])
        @variable(model, 0 ≤ y_copy[indexSets.G] ≤ 1)  
        augmentVar_copy = Dict(
            (g, k) => @variable(model, base_name = "augmentVar_copy[$g, $k]", lower_bound = 0, upper_bound = 1)
            for g in indexSets.G for k in 1:1
        );
        model[:augmentVar_copy] = augmentVar_copy;     
    end

    # constraints for augmented variables: Choosing one leaf node
    @constraint(model, [g in indexSets.G, k in [1]], augmentVar[g, k] == 1)


    ## problem constraints:
    # power flow constraints
    for l in indexSets.L
        i = l[1]
        j = l[2]
        @constraint(model, P[l] ≤ - paramOPF.b[l] * (θ_angle[i] - θ_angle[j]))
        @constraint(model, P[l] ≥ - paramOPF.b[l] * (θ_angle[i] - θ_angle[j]))
    end
    
    # power flow limitation
    @constraint(model, [l in indexSets.L], P[l] ≥ - paramOPF.W[l])
    @constraint(model, [l in indexSets.L], P[l] ≤   paramOPF.W[l])
    # generator limitation
    @constraint(model, [g in indexSets.G], s[g] ≥ paramOPF.smin[g] * y[g])
    @constraint(model, [g in indexSets.G], s[g] ≤ paramOPF.smax[g] * y[g])

    # power balance constraints
    @constraint(model, PowerBalance[i in indexSets.B], 
                    sum(s[g] for g in indexSets.Gᵢ[i]) -
                    sum(P[(i, j)] for j in indexSets.out_L[i]) + 
                    sum(P[(j, i)] for j in indexSets.in_L[i]) 
                    .== sum(paramDemand.demand[d] * x[d] for d in indexSets.Dᵢ[i]) 
    )
    
    # on/off status with startup and shutdown decision
    @constraint(
        model, 
        ShutUpDown[g in indexSets.G], 
        v[g] - w[g] == y[g] - y_copy[g]
    );
    @constraint(
        model, 
        Ramping1[g in indexSets.G], 
        s[g] - s_copy[g] <= paramOPF.M[g] * y_copy[g] + paramOPF.smin[g] * v[g]
    );
    @constraint(
        model, 
        Ramping2[g in indexSets.G], 
        s[g] - s_copy[g] >= - paramOPF.M[g] * y[g] - paramOPF.smin[g] * w[g]
    );

    # production cost
    @constraint(
        model, 
        production[g in indexSets.G, o in keys(paramOPF.slope[g])], 
        h[g] ≥ paramOPF.slope[g][o] * s[g] + paramOPF.intercept[g][o] * y[g]
    );

    # objective function
    @expression(
        model, 
        primal_objective_expression, 
        sum(h[g] + paramOPF.C_start[g] * v[g] + 
        paramOPF.C_down[g] * w[g] for g in indexSets.G) + 
        sum(paramDemand.w[d] * (1 - x[d]) for d in indexSets.D) + 
        sum(θ)
    );

    # @expression(
    #     model, 
    #     dual_objective_expression, 
    #     sum(h[g] + paramOPF.C_start[g] * v[g] + 
    #     paramOPF.C_down[g] * w[g] for g in indexSets.G) + 
    #     sum(paramDemand.w[d] * (1 - x[d]) for d in indexSets.D) + 
    #     sum(θ)
    # );
    @objective(
        model, 
        Min, 
        primal_objective_expression
    );
    
    return SDDPLModel(
                model, 
                Dict{Any, Dict{Any, VariableRef}}(:y => Dict{Any, VariableRef}(g => y[g] for g in indexSets.G)), 
                nothing, 
                Dict{Any, Dict{Any, VariableRef}}(:s => Dict{Any, VariableRef}(g => s[g] for g in indexSets.G)), 
                nothing, 
                Dict(:s => Dict{Any, Dict{Any, Dict{Symbol, Any}}}(
                        g => Dict(
                            k => Dict(
                                :lb => 0.0, 
                                :ub => paramOPF.smax[g], 
                                :parent => nothing, 
                                :sibling => nothing, 
                                :var => augmentVar[g,1]) for k in 1:1
                                ) for g in indexSets.G
                    )
                )
    )
end



"""
forwardPass(ξ): function for forward pass in parallel computing

# Arguments

  1. `ξ`: A sampled scenario path

# Returns
  1. `scenario_solution_collection`: cut coefficients

"""
function forwardPass(
    ξ::Dict{Int64, RandomVariables};
    ModelList::Dict{Int, SDDPLModel} = ModelList,
    paramDemand::ParamDemand = paramDemand, 
    paramOPF::ParamOPF = paramOPF, 
    indexSets::IndexSets = indexSets, 
    initialStateInfo::StateInfo = initialStateInfo
)
    stateInfoList = Dict();
    stateInfoList[0] = deepcopy(initialStateInfo);
    for t in 1:indexSets.T
        ModelModification!( 
            ModelList[t].model, 
            ξ[t],
            paramDemand,
            stateInfoList[t-1];
            indexSets = indexSets
        )
        optimize!(ModelList[t].model);

        # record the solution
        BinVar = Dict{Any, Dict{Any, Any}}(:y => Dict{Any, Any}(
            g => JuMP.value(ModelList[t].BinVar[:y][g]) for g in indexSets.G)
        );
        ContVar = Dict{Any, Dict{Any, Any}}(:s => Dict{Any, Any}(
            g => JuMP.value(ModelList[t].ContVar[:s][g]) for g in indexSets.G)
        );
        ContVarLeaf = Dict{Any, Dict{Any, Dict{Any, Dict{Symbol, Any}}}}(
            :s => Dict{Any, Dict{Any, Dict{Symbol, Any}}}(
                g => Dict(
                    k => Dict(
                        :var => JuMP.value(ModelList[t].ContVarLeaf[:s][g][k][:var])
                    ) for k in keys(ModelList[t].ContVarLeaf[:s][g]) if JuMP.value(ModelList[t].ContVarLeaf[:s][g][k][:var]) > .5
                ) for g in indexSets.G
            )
        );
        
        stageValue = JuMP.objective_value(ModelList[t].model) - sum(JuMP.value.(ModelList[t].model[:θ])); 
        stateValue = JuMP.objective_value(ModelList[t].model);
        ContAugState = Dict{Any, Dict{Any, Dict{Any, Any}}}(:s => Dict{Any, Dict{Any, Any}}(g => Dict{Any, Any}() for g in indexSets.G));
        stateInfoList[t] = StateInfo(
            BinVar, 
            nothing, 
            ContVar, 
            nothing, 
            ContVarLeaf, 
            stageValue, 
            stateValue, 
            nothing, 
            ContAugState
        );
    end  
    return stateInfoList  
end
