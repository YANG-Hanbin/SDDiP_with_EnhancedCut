## binarize stage variable, x = A * L, where L ∈ {0, 1}ⁿ
function intergerBinarization(ū::Vector{Float64})
    row_num = size(ū)[1];

    var_num = floor.(Int, log.(2,ū)) .+ 1; 

    col_num = sum(Int, var_num);

    A = zeros(Int64, row_num, col_num)
 
    for i in 1:row_num
        l = sum(var_num[l] for l in 1:i)
        for j in (l+1-var_num[i]):(l+var_num[i]-var_num[i])
            A[i,j] = 2^(j-(l+1-var_num[i]))
        end
    end
    return BinaryInfo(A, col_num, row_num)
end






################### nonanticipativity for multistage problem #######################

function recursion_scenario_tree(
    pathList::Vector{Int64},                                 
    P::Float64,                                 
    scenario_sequence::Dict{Int64, Dict{Int64, Any}},                                 
    t::Int64;                                   
    Ω::Dict{Int64,Dict{Int64,RandomVariables}} = Ω, 
    prob::Dict{Int64,Vector{Float64}} = prob, 
    T::Int64 = 2
)

    if t ≤ T
        for ω_key in keys(Ω[t])

            pathList_copy = copy(pathList)
            P_copy = copy(P)

            push!(pathList_copy, ω_key)
            P_copy = P_copy * prob[t][ω_key]

            recursion_scenario_tree(pathList_copy, P_copy, scenario_sequence, t+1, Ω = Ω, prob = prob, T = T)
        end
    else
        if haskey(scenario_sequence, 1)
            scenario_sequence[maximum(keys(scenario_sequence))+1] = Dict(1 => pathList, 2 => P)
        else
            scenario_sequence[1] = Dict(1 => pathList, 2 => P)
        end
        return scenario_sequence
    end

end

# scenario_sequence = Dict{Int64, Dict{Int64, Any}}()  ## the first index is for scenario index, the second one is for stage
# pathList = Vector{Int64}()
# push!(pathList, 1)
# P = 1.0

# recursion_scenario_tree(pathList, P, scenario_sequence, 2, T = T)
# scenario_tree = scenario_sequence



## sampling function 
function DrawSamples(scenario_sequence::Dict{Int64, Dict{Int64, Any}})
    # draw f, A, B, C, b from Ωₜ according to distribution P
    P = Vector{Float64}()
    for key in keys(scenario_sequence)
        push!(P, scenario_sequence[key][2])
    end
    items = [i for i in keys(scenario_sequence)]
    weights = Weights(P)
    j = sample(items, weights)
    return j
end


## form scenarios
function SampleScenarios(
    scenario_sequence::Dict{Int64, Dict{Int64, Any}}; 
    T::Int64 = 5, 
    M::Int64 = 30
)
    ## a dict to store realization for each stage t in scenario k
    scenarios = Dict{Int64, Int64}()
    for k in 1:M
          scenarios[k] = DrawSamples(scenario_sequence)
    end
    return scenarios
end

function SampleScenarios(
    Ω::Dict{Int64, Dict{Int64, RandomVariables}}, 
    probList::Dict{Int64, Vector{Float64}}; 
    M::Int64 = 30
)
    ## a dict to store realization for each stage t in scenario k
    T = length(keys(Ω));
    scenarios = Dict{Int64, Vector{Int64}}()
    for k in 1:M
        scenario_path = [1];
        for t in 2:T
            items = [i for i in keys(Ω[t])]
            weights = Weights(probList[t])
            j = sample(items, weights)
            push!(scenario_path, j)
        end
        scenarios[k] = scenario_path
    end
    return scenarios
end

## rounding data
function round!(a::Float64)               ## a = 1.3333e10
    b = floor(log10(a))                   ## b = 10
    c = round(a/10^b,digits = 2)          ## c = 1.33
    d = c * 10^b                          ## d = 1.33e10
    return [b, c, d]
end


## setup coefficients
function dataGeneration(; 
    T::Int64 = 2, 
    r::Float64 = 0.08, ## the annualized interest rate
    N::Matrix{Float64} = N, ## Generator rating
    ū::Vector{Float64} = ū, ## maximum number of each type of generators
    c::Vector{Float64} = c, # c_g from table 4, cost/MW to build a generator of type g
    mg::Vector{Int64} = mg,
    fuel_price::Vector{Float64} = fuel_price,
    heat_rate::Vector{Int64} = heat_rate,
    eff::Vector{Float64} = eff,
    om_cost::Vector{Float64} = om_cost, 
    s₀::Vector{Int64} = s₀,
    penalty::Float64 = 1e10, 
    initial_demand::Float64 = 1e7, 
    seed::Int64 = 1234, num_Ω::Int64 = 10
)::NamedTuple

    binaryInfo = intergerBinarization(ū)
    (A, n, d) = (binaryInfo.A, binaryInfo.n, binaryInfo.d)

    #  compute c1
    c1 = [[c[i]*mg[i]/(1+r)^j for i in 1:6 ]  for j in 1:T ] 
    #  compute c2
    c2 = [[fuel_price[i]*heat_rate[i]*1e-3*eff[i] for i in 1:6]*(1.02)^j + om_cost*(1.03)^j for j in 1:T]


    stageDataList = Dict{Int64,StageData}()
    for t in 1:T 
        stageDataList[t] = StageData(c1[t], c2[t], ū, 8760., N, s₀, penalty)
    end

    ##########################################################################################
    ############################  To generate random variable  ###############################
    ##########################################################################################
    N_rv = Vector{Int64}()  # the number of realization of each stage
    N_rv = [num_Ω for t in 1:T] 
    # N_rv = round.(rand(T) * 10) 

    Random.seed!(seed)

    Ω = Dict{Int64,Dict{Int64,RandomVariables}}()   # each stage t, its node are included in Ω[t]

    for t in 1:T 
        Ω[t] = Dict{Int64,RandomVariables}()
        for i in 1:N_rv[t]
            if t == 1
                Ω[t][i]= RandomVariables([initial_demand])
            else
                Ω[t][i]= RandomVariables( rand(Uniform(1, 1.5))*Ω[t-1][i].d )
            end
        end
    end


    probList = Dict{Int64,Vector{Float64}}()  # P(node in t-1 --> node in t ) = prob[t]
    for t in 1:T 
        # randomVector = round.(rand(N_rv[t]),digits = 2)
        # probList[t] = round.(randomVector/sum(randomVector),digits = 2)
        probList[t] = [1/N_rv[t] for i in 1:N_rv[t]]
    end

    return (probList = probList, 
            stageDataList = stageDataList, 
            Ω = Ω, 
            binaryInfo = binaryInfo)
end



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
    println("Iter |        LB        |        UB        |       Gap      |      i-time     |    #D.     |     T-Time")
    println("----------------------------------------------------------------------------------------------------------")
    return 
end

function save_info(
    param::NamedTuple, 
    sddpResults::Dict;
    logger_save::Bool = true
)::Nothing
    if logger_save == true
        cutSelection = param.cutSelection; num = param.num; T = param.T; tightness = param.tightness; algorithm = param.algorithm;
        save(
            "/Users/aaron/SDDiP_with_EnhancedCut/src/GenerationExpansion/logger/Periods$T-Real$num/$algorithm-$cutSelection-$tightness.jld2", 
            "sddpResults", 
            sddpResults
        );
    end
    return 
end


function param_setup(;
    terminate_time::Any             = 3600,
    terminate_threshold::Float64    = 1e-3,
    MaxIter::Int64                  = 3000,
    M::Int64                        = 5, 
    ε::Float64                      = 0.125,
    tightness::Bool                 = true,
    cutSelection::String            = "LC", 
    T::Int64                        = 12,
    num::Int64                      = 10,
    Output_Gap::Bool                = false,
    ℓ1::Float64                     = 1.0,
    ℓ2::Float64                     = 1.0,
    nxt_bound::Float64              = 1e8,
    logger_save::Bool               = true,
    algorithm::Symbol               = :SDDiP
)::NamedTuple
    return (
        terminate_time      = terminate_time,
        terminate_threshold = terminate_threshold,
        MaxIter             = MaxIter,
        M                   = M, 
        ε                   = ε,
        tightness           = tightness,
        cutSelection        = cutSelection, 
        Output_Gap          = Output_Gap,
        T                   = T, 
        num                 = num, 
        ℓ1                  = ℓ1,
        ℓ2                  = ℓ2,
        nxt_bound           = nxt_bound,
        logger_save         = logger_save,
        algorithm           = algorithm  
    )
end