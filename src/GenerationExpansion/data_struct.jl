#############################################################################################
####################################   Data Structure   #####################################
#############################################################################################
struct CutCoefficient
    v   ::Dict{Int64,Dict{Int64, Float64}} # [i][k] where i is the iteration index, and k is the scenario index
    π   ::Dict{Int64,Dict{Int64, Vector{Float64}}}  # [[1.,2.,3.],[1.,2.,3.]]  -- push!(Π, π) to add a new element
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
    l               :: Vector{VariableRef}
    y               :: Vector{VariableRef}
    θ               :: VariableRef
    demand          :: Vector{Float64}
    slack           :: VariableRef
    sum_generator   :: Vector{Float64}
end


struct BackwardModelInfo
    model           :: Model
    l               :: Vector{VariableRef}
    L               :: Vector{VariableRef}
    y               :: Vector{VariableRef}
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



## sampling function 
function DrawSamples(P)
    # draw f, A, B, C, b from Ωₜ according to distribution P
    items = [i for i in 1:length(P)]
    weights = Weights(P)
    j = sample(items, weights)
    return j
end

## form scenarios
function SampleScenarios( ;T::Int64 = 5, M::Int64 = 30)
    ## a dict to store realization for each stage t in scenario k
    scenarios = Dict()
    for k in 1:M 
        for t in 1:T 
          ## draw sample
          scenarios[k, t] = DrawSamples(prob[t])
        end
    end
    return scenarios
end


