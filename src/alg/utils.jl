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
    terminate_time::Any = 3600,
    TimeLimit::Any = 10,
    terminate_threshold::Float64 = 1e-3,
    ε::Float64 = 0.125,
    MaxIter::Int64 = 3000,
    tightness::Bool = true,
    numScenarios::Int64 = 3,
    LiftIterThreshold::Int64 = 10,
    branch_threshold::Float64 = 1e-3,
    cutSelection::Symbol = :PLC, 
    algorithm::Symbol = :SDDPL,
    T::Int64 = 12,
    num::Int64 = 10,
    med_method::Symbol = :ExactPoint,
    case::String = "case30pwl"
)::NamedTuple

    return (
        verbose             = false,
        MIPGap              = 1e-4,
        TimeLimit           = TimeLimit,
        terminate_time      = terminate_time,
        terminate_threshold = terminate_threshold,
        MaxIter             = MaxIter,
        θ̲                   = 0.0,
        OPT                 = 0.0,
        ε                   = ε,
        κ                   = Dict{Int64, Int64}(),
        tightness           = tightness,
        numScenarios        = numScenarios,
        LiftIterThreshold   = LiftIterThreshold,
        branch_threshold    = branch_threshold,
        ## "interval_mid", "exact_point"
        med_method          = med_method,   
        ## :PLC, :SMC, :LC
        cutSelection        = cutSelection,             
        ## :SDDPL, :SDDP, :SDDiP
        algorithm           = algorithm,   
        T                   = T, 
        num                 = num, 
        case                = case, # "case_RTS_GMLC", "case30", "case30pwl", "case24_ieee_rts"     
    )
end


function param_levelsetmethod_setup(;
    μ::Float64 = 0.9,
    λ::Float64 = 0.5,
    threshold::Float64 = 1e-4,
    nxt_bound::Float64 = 1e10,
    MaxIter::Int64 = 200,
    verbose::Bool = false
)::NamedTuple
    return (
        μ             = μ,
        λ             = λ,
        threshold     = threshold,
        nxt_bound     = nxt_bound,
        MaxIter       = MaxIter,
        verbose       = verbose,
    )
end


function param_cut_setup(;
    core_point_strategy::String = "Eps", # "Mid", "Eps"
    δ::Float64 = 1e-3,
    ℓ::Float64 = .0,
)::NamedTuple
    return (
        core_point_strategy = core_point_strategy, # "Mid", "Eps"
        δ                   = δ,
        ℓ                   = ℓ,
    )
end