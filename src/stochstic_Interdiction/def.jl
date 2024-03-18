#############################################################################################
####################################   Data Structure   #####################################
#############################################################################################
struct LagrangianCut
    v               ::Float64                         # [i][k] where i is the iteration index, and k is the scenario index
    π               ::Vector{Float64}                 # [[1.,2.,3.],[1.,2.,3.]]  -- push!(Π, π) to add a new element
    x               ::Union{Vector{Real}, Nothing}
end

# the scenario data
struct ScenarioData 
    s        ::  Int64                                # s is the origin of the intruder under scenario ω
    t        ::  Int64                                # t is the desination of the intruder under scenario ω
    ϕ        ::  Dict{Int64, Real}                    # ϕ[j] is the value of maximum-reliablility path from j to t
end


struct IndexSets
    N        ::  Vector{Int64}                       # set of nodes
    A        ::  Vector{Tuple{Int64, Int64}}         # set of arcs of the network
    D        ::  Vector{Tuple{Int64, Int64}}         # subset of arcs on which sensors are allowed to be placed
    Dᶜ       ::  Vector{Tuple{Int64, Int64}}         # subset of arcs on which sensors are NOt allowed to be placed
    Ω        ::  Vector{Int64}                       # set of scenarios
end




# the first-stage data
struct StageData 
    c        ::  Dict{Tuple{Int64, Int64}, Real}     # c is the cost of installing a sensor on arc (i,j)
    b        ::  Real                                # total budge
end



# the probability data
struct Prob 
    p        ::  Dict{Int64, Real}                   # p[ω] is the realization probability of scenario ω
    r        ::  Dict{Tuple{Int64, Int64}, Real}     # r[i,j] is the probability of intruder avoiding detection with a sensor
    q        ::  Dict{Tuple{Int64, Int64}, Real}     # q[i,j] is the probability of intruder avoiding detection without a sensor
end


struct Forward_stage1_Info
    model           :: Model
    x               :: Vector{VariableRef} 
    θ               :: Vector{VariableRef}
end

struct Forward_stage2_Info
    model           :: Model
    x               :: Vector{VariableRef} 
    π               :: Vector{VariableRef}
end


struct Backward_stage2_Info
    model           :: Model
    x               :: Vector{VariableRef} 
    π               :: Vector{VariableRef}
end



struct BackwardModelInfo
    model           :: Model
    x               :: Vector{VariableRef} ## for current state, x is the number of generators
    Lt              :: Vector{VariableRef} ## stage variable, A * Lt is total number of generators built
    Lc              :: Vector{VariableRef} ## local cooy variable
    y               :: Vector{VariableRef} ## amount of electricity
    θ               :: VariableRef
    # demand          :: Vector{Float64}
    slack           :: VariableRef
    # sum_generator   :: Vector{Float64}
end



struct GurobiModelInfo
    model           :: Model
    x               :: Array{VariableRef} ## for current state, x is the number of generators
    y               :: Array{VariableRef} ## amount of electricity
    slack           :: Matrix{VariableRef}
    num_Ω           :: Int64
end


## data structure for levelset method
## data structure for levelset method
mutable struct FunctionHistory
    f_his        :: Dict{Int64, Float64}          ## record f(x_j)     
    G_max_his    :: Dict{Int64, Float64}          ## record max(g[k] for k in 1:m)(x_j)   
end


mutable struct CurrentInfo
    x            :: Vector{Float64}                 ## record x point
    f            :: Float64                         ## record f(x_j)
    G            :: Dict{Int64, Float64} 
    df           :: Vector{Float64}
    dG           :: Dict{Int64, Vector{Float64}}    ## actually is a matrix.  But we use dict to store it
    L_at_x̂       :: Float64                         ## only for solving dual problem
end




struct ModelInfo
    model :: Model
    x     :: Vector{VariableRef}
    y     :: VariableRef
    z     :: VariableRef
end




struct LevelSetMethodParam
    μ             ::Float64                         ## param for adjust α
    λ             ::Union{Float64, Nothing}         ## param for adjust level
    threshold     ::Float64                         ## threshold for Δ
    nxt_bound     ::Float64                         ## lower bound for solving next iteration point π
    max_iter      ::Int64     
    Output        ::Int64                           ## Gurobi Output parameter
    Output_Gap    ::Bool                            ## if True will print Δ info

    ## the following parameters are problem-specific
    x̂             ::Vector{Float64}                 ## first stage solution
    cutSelection  ::String
    x̃             ::Union{Vector{Float64}, Nothing} ## interior point
    f_star_value  ::Union{Float64, Nothing}         ## subproblem optimal value
end


