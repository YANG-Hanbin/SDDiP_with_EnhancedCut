"""
ModelModification!(; model::Model = model)

# Arguments

    1. `model::Model` : a forward pass model of stage t
    2. `randomVariables::RandomVariables` : random variables
    3. `paramDemand::ParamDemand` : demand parameters
    4. `stateInfo::StateInfo` : the last stage decisions
  
# Modification
    1. Remove the other scenario's demand balance constraints
    2. Add the current scenario's demand balance constraints
    3. Update its last stage decision with
"""
function ModelModification!( 
    model::Model, 
    randomVariables::RandomVariables,
    paramDemand::ParamDemand,
    stateInfo::StateInfo;
    indexSets::IndexSets = indexSets
)::Nothing
    if :ContVarNonAnticipative ∉ keys(model.obj_dict)
        @constraint(
            model, 
            ContVarNonAnticipative[g in indexSets.G], 
            model[:s_copy][g] == stateInfo.ContVar[:s][g]
        );
    end

    if :BinVarNonAnticipative ∉ keys(model.obj_dict)
        @constraint(
            model, 
            BinVarNonAnticipative[g in indexSets.G], 
            model[:y_copy][g] == stateInfo.BinVar[:y][g]
        );
    end

    # power balance constraints
    for i in indexSets.B
        delete(model, model[:PowerBalance][i])
    end
    unregister(model, :PowerBalance)
    @constraint(model, PowerBalance[i in indexSets.B], 
                            sum(model[:s][g]      for g in indexSets.Gᵢ[i]) -
                            sum(model[:P][(i, j)] for j in indexSets.out_L[i]) + 
                            sum(model[:P][(j, i)] for j in indexSets.in_L[i]) .==
                            sum(paramDemand.demand[d] * randomVariables.deviation[d] * model[:x][d] for d in indexSets.Dᵢ[i])
    )

    @objective(
        model, 
        Min, 
        model[:primal_objective_expression]
    );
    return
end

"""
function sample_scenarios(; 
    numScenarios::Int64 = 10, 
    scenarioTree::ScenarioTree = scenarioTree
)

# Arguments

    1. `numScenarios`: The number of scenarios will be sampled
    2. `scenarioTree`: A scenario tree

# Returns
    1. `Ξ`: A subset of scenarios.
"""
function sample_scenarios(; 
    numScenarios::Int64 = 10, 
    scenarioTree::ScenarioTree = scenarioTree
)
    Ξ = Dict{Int64, Dict{Int64, RandomVariables}}()
    for ω in 1:numScenarios
        ξ = Dict{Int64, RandomVariables}()
        ξ[1] = scenarioTree.tree[1].nodes[1]
        n = wsample(
            collect(keys(scenarioTree.tree[1].prob)), 
            collect(values(scenarioTree.tree[1].prob)), 
            1
        )[1]
        for t in 2:length(keys(scenarioTree.tree))
            ξ[t] = scenarioTree.tree[t].nodes[n]
            n = wsample(
                collect(keys(scenarioTree.tree[t].prob)), 
                collect(values(scenarioTree.tree[t].prob)), 
            1)[1]
        end
        Ξ[ω] = ξ
    end
    return Ξ
end

"""

function print_iteration_info(
    i::Int, 
    LB::Float64, 
    UB::Float64,
    gap::Float64, 
    iter_time::Float64, 
    LM_iter::Int, 
    total_time::Float64
)

# Arguments

    1. `i`: The current iteration number
    2. `LB`: The lower bound
    3. `UB`: The upper bound
    4. `gap`: The gap between the lower and upper bounds
    5. `iter_time`: The time spent on the current iteration
    6. `LM_iter`: The number of Lagrangian multipliers updated in the current iteration
    7. `total_time`: The total time spent on the algorithm

# Prints

"""
function print_iteration_info(
    i::Int64, 
    LB::Float64, 
    UB::Float64,
    gap::Float64, 
    iter_time::Float64, 
    LM_iter::Int, 
    total_Time::Float64
)::Nothing
    @printf("%4d | %12.2f     | %12.2f     | %9.2f%%     | %9.2f s     | %6d     | %10.2f s     \n", 
                i, LB, UB, gap, iter_time, LM_iter, total_Time); 
    return 
end

function print_iteration_info_bar()::Nothing
    println("------------------------------------------ Iteration Info ------------------------------------------------")
    println("Iter |        LB        |        UB        |       Gap      |      i-time     |    #LM     |     T-Time")
    println("----------------------------------------------------------------------------------------------------------")
    return 
end

function save_info(
    param::NamedTuple, 
    sddpResults::Dict
)::Nothing
    case = param.case; cutSelection = param.cutSelection; num = param.num; T = param.T; algorithm = param.algorithm;
    save(
        "/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-$cutSelection.jld2", 
        "sddpResults", 
        sddpResults
    );
    return 
end

function param_setup(;
    terminate_threshold::Float64 = 1e-3,
    MaxIter::Int64 = 3000,
    tightness::Bool = true,
    numScenarios::Int64 = 3,
    LiftIterThreshold::Int64 = 10,
    branch_threshold::Float64 = 1e-3,
    cutSelection::Symbol = :PLC, 
    algorithm::Symbol = :SDDPL,
    T::Int64 = 12,
    num::Int64 = 10,
    case::String = "case30pwl"
)::NamedTuple
    param = (
        verbose             = false,
        MIPGap              = 1e-4,
        TimeLimit           = 3600,
        terminate_threshold = terminate_threshold,
        MaxIter             = MaxIter,
        θ̲                   = 0.0,
        OPT                 = 0.0,
        tightness           = tightness,
        numScenarios        = numScenarios,
        LiftIterThreshold   = LiftIterThreshold,
        branch_threshold    = branch_threshold,
        ## "interval_mid", "exact_point"
        med_method          = "interval_mid",   
        ## :PLC, :SMC, :LC
        cutSelection        = cutSelection,             
        ## :SDDPL, :SDDP, :SDDiP
        algorithm           = algorithm,   
        T                   = T, 
        num                 = num, 
        case                = case, # "case_RTS_GMLC", "case30", "case30pwl", "case24_ieee_rts"     
    )
    return param
end