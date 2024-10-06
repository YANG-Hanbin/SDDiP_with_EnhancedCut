# define data structure

## ====================================================================================== ##
## =================================DC-OPF: Index Sets ================================== ##
## ====================================================================================== ##
mutable struct IndexSets
    D       :: Vector{Int64}                    ## set of load demand
    G       :: Vector{Int64}                    ## set of generators
    L       :: Vector{Tuple{Int64, Int64}}      ## set of transmission lines
    B       :: Vector{Int64}                    ## set of buses
    T       :: Int64                            ## set of time periods  1:T
    Dᵢ      :: Dict{Int64,Vector{Int64}}       
    Gᵢ      :: Dict{Int64,Vector{Int64}}
    out_L   :: Dict{Int64,Vector{Int64}} 
    in_L    :: Dict{Int64,Vector{Int64}}
end


mutable struct ParamDemand
    demand        :: Dict{Int64, Float64}                    ## mean of load demand for d ∈ D
    w             :: Dict{Int64, Float64}                    ## priority level of load D
    cb            :: Dict{Int64, Float64}                    ## set of fire damage cost cᵢ at :b ∈ B
    cg            :: Dict{Int64, Float64}                    ## set of fire damage cost cᵢ at :g ∈ G
    cl            :: Dict{Tuple{Int64, Int64}, Float64}      ## set of fire damage cost cᵢ at :l ∈ L
    penalty       :: Float64                                 ## penalty paramter for constraints b and c
end


mutable struct ParamOPF  ## included in a period dict
    b           :: Dict{Tuple{Int64, Int64}, Float64}       ## :l ∈ L  total line charging susceptance
    θmax        :: Float64                                  ## angle difference
    θmin        :: Float64
    W           :: Dict{Tuple{Int64, Int64}, Float64}       ## :l ∈ L
    smax        :: Dict{Int64, Float64}                     ## :g ∈ G  Pmax
    smin        :: Dict{Int64, Float64}                     ## :g ∈ G
    M           :: Dict{Int64, Float64}                     ## :g ∈ G  multi-period ramping limit 
    slope       :: Dict{Int64, Dict{Int64, Float64}}        ## :g ∈ G, :o ∈ O  cost function slope
    intercept   :: Dict{Int64, Dict{Int64, Float64}}        ## :g ∈ G, :o ∈ O  cost function intercept
    C_start     :: Dict{Int64, Float64}                     ## :g ∈ G
    C_down      :: Dict{Int64, Float64}                     ## :g ∈ G
end


## for period t with realization ω
struct RandomVariables  
    deviation        :: Dict{Int64, Float64}                   ## :d ∈ D, realization deviation from the mean demand
end

struct StageRealization
    nodes            :: Dict{Int64, RandomVariables}           ## :n ∈ N_t, nodes of stage t
    prob             :: Dict{Int64, Float64}                   ## :n ∈ N_{t+1}, probabilities for next stage nodes
end

struct ScenarioTree
    tree      :: Dict{Int64, StageRealization}                ## :t ∈ [1, ..., T], realization percentage of the demand
end

## ====================================================================================== ##
## ================================== Cut Coefficient =================================== ##
## ====================================================================================== ##
struct CutCoefficient
    v               ::Dict{Int64,Dict{Int64, Float64}} # [i][k] where i is the iteration index, and k is the scenario index
    πb               ::Dict{Int64,Dict{Int64, Vector{Float64}}}  # [[1.,2.,3.],[1.,2.,3.]]  -- push!(Π, π) to add a new element
    πg               ::Dict{Int64,Dict{Int64, Vector{Float64}}}
    πl               ::Dict{Int64,Dict{Int64, Vector{Float64}}}
end



## ====================================================================================== ##
## ================================== Levelset Method =================================== ##
## ====================================================================================== ##

## data structure for levelset method
mutable struct FunctionHistory
    f_his        :: Dict{Int64, Float64}          ## record f(x_j)     
    G_max_his    :: Dict{Int64, Float64}          ## record max(g[k] for k in 1:m)(x_j)   
end


mutable struct CurrentInfo
    x            :: Dict                          ## record x point
    f            :: Float64                       ## record f(x_j)
    G            :: Dict
    df           :: Dict
    dG           :: Dict                          ## actually is a matrix.  But we use dict to store it
end


struct ModelInfo
    model :: Model
    xλ    :: JuMP.Containers.SparseAxisArray{VariableRef, 2, Tuple{Int64, Int64}}
    xy    :: JuMP.Containers.DenseAxisArray{VariableRef}
    y     :: VariableRef
    z     :: VariableRef
end


struct LevelSetMethodParam
    μ             ::Float64                     ## param for adjust α
    λ             ::Union{Float64, Nothing}     ## param for adjust level
    threshold     ::Float64                     ## threshold for Δ
    nxt_bound     ::Float64                     ## lower bound for solving next iteration point π
    max_iter      ::Int64     
    Output        ::Int64                       ## Gurobi Output parameter
    Output_Gap    ::Bool                        ## if True will print Δ info
    valValFun     ::Float64                     ## value of the approximate value function
    # Adj           ::Union{Bool, Nothing}      ## whether adjust oracle lower bound
end


mutable struct BackwardInfo
    model           ::Model
    x               ::JuMP.Containers.DenseAxisArray{VariableRef, 2, Tuple{Vector{Int64}, Base.OneTo{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}, JuMP.Containers._AxisLookup{Base.OneTo{Int64}}}}
    νb              ::JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}
    νg              ::JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}
    νl              ::JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Tuple{Int64, Int64}}}, Tuple{JuMP.Containers._AxisLookup{Dict{Tuple{Int64, Int64}, Int64}}}}
    zg              ::JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}
    zb              ::JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}
    zl              ::JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Tuple{Int64, Int64}}}, Tuple{JuMP.Containers._AxisLookup{Dict{Tuple{Int64, Int64}, Int64}}}} 
end


mutable struct ForwardInfo
    model           ::Model
    θ               ::JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}
    zg              ::JuMP.Containers.DenseAxisArray{VariableRef, 2, Tuple{Vector{Int64}, Base.OneTo{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}, JuMP.Containers._AxisLookup{Base.OneTo{Int64}}}}
    zb              ::JuMP.Containers.DenseAxisArray{VariableRef, 2, Tuple{Vector{Int64}, Base.OneTo{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}, JuMP.Containers._AxisLookup{Base.OneTo{Int64}}}}
    zl              ::JuMP.Containers.DenseAxisArray{VariableRef, 2, Tuple{Vector{Tuple{Int64, Int64}}, Base.OneTo{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Tuple{Int64, Int64}, Int64}}, JuMP.Containers._AxisLookup{Base.OneTo{Int64}}}}
end


mutable struct Forward2Info
    model           :: Model
    yb              :: JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}
    yg              :: JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}
    yl              :: JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Tuple{Int64, Int64}}}, Tuple{JuMP.Containers._AxisLookup{Dict{Tuple{Int64, Int64}, Int64}}}}
    νb              :: JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}
    νg              :: JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}
    νl              :: JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Tuple{Int64, Int64}}}, Tuple{JuMP.Containers._AxisLookup{Dict{Tuple{Int64, Int64}, Int64}}}}
  end