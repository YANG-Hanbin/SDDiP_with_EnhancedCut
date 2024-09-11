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


struct ParamOPF  ## included in a period dict
    b           :: Dict{Tuple{Int64, Int64}, Float64}       ## :l ∈ L  total line charging susceptance
    θmax        :: Float64                                  ## angle difference
    θmin        :: Float64
    W           :: Dict{Tuple{Int64, Int64}, Float64}       ## :l ∈ L
    smax        :: Dict{Int64, Float64}                     ## :g ∈ G  Pmax
    smin        :: Dict{Int64, Float64}                     ## :g ∈ G
    M           :: Dict{Int64, Float64}                     ## :g ∈ G  multi-period ramping limit 
    slope       :: Dict{Int64, Float64}                     ## :g ∈ G  cost function
    intercept   :: Dict{Int64, Float64}                     ## :g ∈ G
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