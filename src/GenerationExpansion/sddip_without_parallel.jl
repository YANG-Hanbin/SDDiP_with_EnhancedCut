# using Distributed
# addprocs(4)

using JuMP, Test, Statistics, StatsBase, Gurobi, Distributed, ParallelDataTransfer, Random

const GRB_ENV = Gurobi.Env()


include("data_struct.jl")
include("backward_pass.jl")
include("forward_pass.jl")
include("data_file.jl")

#############################################################################################
####################################    main function   #####################################
#############################################################################################
max_iter = 20; ϵ = 1e-2; Enhand_Cut = false

# @broadcast max_iter = 20; ϵ = 1e-2; 
# @broadcast Enhand_Cut = false;

function SDDiP_algorithm(Ω::Dict{Int64,Dict{Int64,RandomVariables}}, prob::Dict{Int64,Vector{Float64}}, StageCoefficient::Dict{Int64, StageData}; 
    ϵ::Float64 = 0.001, M::Int64 = 30, max_iter::Int64 = 20, Enhand_Cut::Bool = true, n::Int64 = 2, d::Int64 = 1, A::Matrix{Int64} = [1 2;])
    ## d: x dim
    ## M: num of scenarios when doing one iteration
    
    T = length(keys(Ω))
    
    i = 1
    LB = - Inf 
    UB = Inf
    
    cut_collection = Dict{Int64, CutCoefficient}()  # here, the index is stage t

    for t in 1:T

        cut_collection[t] = CutCoefficient(
                                            Dict(1=> Dict(1=> 0.0), 2 => Dict()), 
                                            Dict(1=> Dict(1=> zeros(Float64, n)), 2 => Dict()) 
                                          )
    end

    add_result = Dict()

    while true
        M = 100
        Sol_collection = Dict()  # to store every iteration results
        u = Vector{Float64}(undef, M)  # to compute upper bound
        
        Random.seed!(i*3)
        Scenarios = SampleScenarios(T = T, M = M)
        
        ## Forward Step
        for k in 1:M
            sum_generator= [0.0 for i in 1:n]
            for t in 1:T
                Sol_collection[t, k] = forward_step_optimize!(StageCoefficient[t], 
                                                                Ω[t][Scenarios[k, t]].d, 
                                                                sum_generator, 
                                                                cut_collection[t], 
                                                                A = A, d = d, n = n)
                sum_generator = Sol_collection[t,k][1]
            end
            u[k] = sum(Sol_collection[t, k][4] for t in 1:T)
        end

        # @passobj 1 workers() Sol_collection

        ## compute the upper bound
        μ̄ = mean(u)
        σ̂² = var(u)
        UB = μ̄ + 1.96 * sqrt(σ̂²/M)

        ##################################### Parallel Computation for backward step ###########################

        M = 4

        for t = reverse(2:T)
            cut_collection[t-1].v[i] = Dict()
            cut_collection[t-1].π[i] = Dict()
            for k in 1:M 
                c = [0, zeros(Float64,n)]
                for j in keys(Ω[t])
                    @info "$t $k $j"
                    λ_value = .5
                    ϵ_value = 1e-3
                    nxt_bound = 1e14
                    Adj = false
                    @time c = c + prob[t][j] * LevelSetMethod_optimization!(StageCoefficient[t], 
                                                                                    Ω[t][j].d, 
                                                                                    Sol_collection[t-1,k][1], 
                                                                                    cut_collection[t], 
                                                                                    max_iter = 4000, 
                                                                                    Enhand_Cut = Enhand_Cut,  
                                                                                    nxt_bound = nxt_bound,
                                                                                    μ = 0.95, λ = λ_value, ϵ = ϵ_value,
                                                                                    Output_Gap = true, Adj = Adj, threshold = 1e-3, 
                                                                                    A = A, d = d, n = m ) 
                end
                # add cut
                cut_collection[t-1].v[i][k] = c[1]
                cut_collection[t-1].π[i][k] = c[2]
            end
        end

        ########################################################################################################



        ## compute the lower bound
        _LB = forward_step_optimize!(StageCoefficient[1], Ω[1][1].d, [0.0 for i in 1:n], cut_collection[1], A = A, n = n, d = d)
        LB = _LB[3] + _LB[4]
        add_result[i] = _LB  # store the result
        i = i + 1

        @info "LB is $LB, UB is $UB"
        if UB-LB <= ϵ * UB || i > max_iter
            return [LB,UB]
        end

    end

end































