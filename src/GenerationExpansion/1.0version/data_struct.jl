#############################################################################################
####################################   Data Structure   #####################################
#############################################################################################
struct CutCoefficient
    v               ::Dict{Int64,Dict{Int64, Float64}} # [i][k] where i is the iteration index, and k is the scenario index
    π               ::Dict{Int64,Dict{Int64, Vector{Float64}}}  # [[1.,2.,3.],[1.,2.,3.]]  -- push!(Π, π) to add a new element
    # Enhand_Cut      ::Bool
end



struct RandomVariables
    d   ::Vector{Float64}
end



struct StageData ## with the assumption that only b has stochasticity
    c1       ::Vector{Float64}
    c2       ::Vector{Float64}
    ū        ::Vector{Float64}
    h        ::Float64
    N        ::Matrix{Float64}
    s₀       ::Vector{Float64}
    penalty  ::Float64
end



struct ForwardModelInfo
    model           :: Model
    x               :: Vector{VariableRef} ## for current state, x is the number of generators
    Lt              :: Vector{VariableRef} ## stage variable, A * Lt is total number of generators built
    y               :: Vector{VariableRef} ## amount of electricity
    θ               :: VariableRef
    demand          :: Vector{Float64}
    slack           :: VariableRef
    sum_generator   :: Vector{Float64}
end



struct BackwardModelInfo
    model           :: Model
    x               :: Vector{VariableRef} ## for current state, x is the number of generators
    Lt              :: Vector{VariableRef} ## stage variable, A * Lt is total number of generators built
    Lc              :: Vector{VariableRef} ## local cooy variable
    y               :: Vector{VariableRef} ## amount of electricity
    θ               :: VariableRef
    demand          :: Vector{Float64}
    slack           :: VariableRef
    sum_generator   :: Vector{Float64}
end



## data structure for levelset method
mutable struct FunctionInfo
    x_his        :: Dict{Int64, Vector{Float64}}  ## record every x_j point
    G_max_his    :: Dict{Int64, Float64}          ## record max(g[k] for k in 1:m)(x_j)
    f_his        :: Dict{Int64, Float64}          ## record f(x_j)
    df           :: Vector{Float64}
    dG           :: Dict{Int64, Vector{Float64}}  ## actually is a matrix.  But we use dict to store it
    G            :: Dict{Int64, Float64}          
end




struct ModelInfo
    model :: Model
    x     :: Vector{VariableRef}
    y     :: VariableRef
    z     :: VariableRef
end




## binarize stage variable, x = A * L, where L ∈ {0, 1}ⁿ
function binarize_gen(ū::Vector{Float64})
    row_num = size(ū); row_num = row_num[1];

    var_num = floor.(Int, log.(2,ū) ) .+ 1; col_num = sum(Int, var_num);

    A = zeros(Int64, row_num, col_num)
 
    for i in 1:row_num
        l = sum(var_num[l] for l in 1:i)
        for j in (l+1-var_num[i]):(l+var_num[i]-var_num[i])
            A[i,j] = 2^(j-(l+1-var_num[i]))
        end
    end

    return Dict(1 =>A, 2 =>col_num, 3=> row_num)
end

# binaryDict = binarize_gen(ū)
# (A, n, d) = (binaryDict[1], binaryDict[2], binaryDict[3])







################### nonanticipativity for multistage problem #######################

function recursion_scenario_tree(pathList::Vector{Int64}, P::Float64, scenario_sequence::Dict{Int64, Dict{Int64, Any}}, t::Int64;   
    Ω::Dict{Int64,Dict{Int64,RandomVariables}} = Ω, prob::Dict{Int64,Vector{Float64}} = prob, T::Int64 = 2)

    if t <= T
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
function SampleScenarios(scenario_sequence::Dict{Int64, Dict{Int64, Any}}; T::Int64 = 5, M::Int64 = 30)
    ## a dict to store realization for each stage t in scenario k
    scenarios = Dict{Int64, Int64}()
    for k in 1:M
          scenarios[k] = DrawSamples(scenario_sequence)
    end
    return scenarios
end







