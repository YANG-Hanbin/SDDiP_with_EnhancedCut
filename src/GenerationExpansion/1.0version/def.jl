#############################################################################################
####################################   Data Structure   #####################################
#############################################################################################
struct CutCoefficient
    v               ::Dict{Int64,Dict{Int64, Float64}} # [i][k] where i is the iteration index, and k is the scenario index
    π               ::Dict{Int64,Dict{Int64, Vector{Float64}}}  # [[1.,2.,3.],[1.,2.,3.]]  -- push!(Π, π) to add a new element
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



struct GurobiModelInfo
    model           :: Model
    x               :: Array{VariableRef} ## for current state, x is the number of generators
    y               :: Array{VariableRef} ## amount of electricity
    slack           :: Matrix{VariableRef}
    num_Ω           :: Int64
end


## data structure for levelset method
mutable struct FunctionInfo
    x_his        :: Dict{Int64, Vector{Float64}}  ## record every x_j point
    G_max_his    :: Dict{Int64, Float64}          ## record max(g[k] for k in 1:m)(x_j)
    f_his        :: Dict{Int64, Float64}          ## record f(x_j)
    df           :: Vector{Float64}
    dG           :: Dict{Int64, Vector{Float64}}  ## actually is a matrix.  But we use dict to store it
    G            :: Dict{Int64, Float64}          ## the index k is k-th constraint
end




struct ModelInfo
    model :: Model
    x     :: Vector{VariableRef}
    y     :: VariableRef
    z     :: VariableRef
end




struct BinaryInfo
    A     ::Matrix{Int64}
    n     ::Int64
    d     ::Int64
end






struct LevelSetMethodParam
    μ             ::Float64   ## param for adjust α
    λ             ::Float64   ## param for adjust level
    threshold     ::Float64   ## threshold for Δ
    nxt_bound     ::Float64   ## lower bound for solving next iteration point π
    max_iter      ::Int64     
    Output        ::Int64     ## Gurobi Output parameter
    Output_Gap    ::Bool      ## if True will print Δ info
    Adj           ::Bool      ## whether adjust oracle lower bound
end


