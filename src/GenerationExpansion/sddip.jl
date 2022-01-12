using JuMP, Test, Statistics, StatsBase, Gurobi, Distributed, Random  #, Flux

const GRB_ENV = Gurobi.Env()

include("/Users/aaron/SDDiP_with_EnhancedCut/src/GenerationExpansion/data_struct.jl")
include("/Users/aaron/SDDiP_with_EnhancedCut/src/GenerationExpansion/backward_pass.jl")
include("/Users/aaron/SDDiP_with_EnhancedCut/src/GenerationExpansion/forward_pass.jl")


#############################################################################################
####################################    main function   #####################################
#############################################################################################

function SDDiP_algorithm(Ω::Dict{Int64,Dict{Int64,RandomVariables}}, prob::Dict{Int64,Vector{Float64}}, StageCoefficient::Dict{Int64, StageData}; 
    ϵ::Float64 = 0.001, M::Int64 = 30, d::Int64 = 6, max_iter::Int64 = 2000, Enhand_Cut::Bool = true)
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
                                            Dict(1=> Dict(1=> zeros(Float64, d)), 2 => Dict()) 
                                          )
    end

    while true
        Sol_collection = Dict()  # to store every iteration results
        u = Vector{Float64}(undef, M)  # to compute upper bound
        
        ## Forward Step
        for k in 1:M
            sum_generator= [0.0 for i in 1:d]
            for t in 1:T
                demand = Ω[t][DrawSamples(prob[t])].d
                Sol_collection[t, k] = forward_step_optimize!(StageCoefficient[t], demand, sum_generator, cut_collection[t])
                sum_generator = Sol_collection[t,k][1]
            end
            u[k] = sum(Sol_collection[t, k][4] for t in 1:T)
        end


        ## compute the upper bound
        μ̄ = mean(u)
        σ̂² = var(u)
        UB = μ̄ + 1.96 * sqrt(σ̂²/M)

        ## Backward Step
        for t = reverse(2:T)
            cut_collection[t-1].v[i] = Dict()
            cut_collection[t-1].π[i] = Dict()
            for k in 1:M 
                c = [0, zeros(Float64,d)]
                for j in keys(Ω[t])
                    @info "$t $k $j"
                    λ_value = .5
                    ϵ_value = 1e-3
                    # compute average cut coefficient
                    demand = Ω[t][j].d
                    if t == 5
                        oracle_bound = 1e13
                        nxt_bound = 1e14
                        if demand[1] > 3.1e8
                            λ_value = .9
                            ϵ_value = .1
                        elseif 2.9e8 <= demand[1] <= 3.1e8
                            λ_value = .7
                            ϵ_value = .1
                        elseif 2.9e8 > demand[1]
                            λ_value = .3
                            ϵ_value = .1
                        end
                    end

                    @time c = c + prob[t][j] * LevelSetMethod_generator_v(StageCoefficient[t], 
                                                                                    Ω[t][j].d, 
                                                                                    Sol_collection[t-1,k][1], 
                                                                                    cut_collection[t], 
                                                                                    max_iter = 3000, 
                                                                                    Enhand_Cut = Enhand_Cut, 
                                                                                    oracle_bound = oracle_bound,
                                                                                    nxt_bound = nxt_bound,
                                                                                    μ = 0.93, λ = λ_value, ϵ = ϵ_value,
                                                                                    Output_Gap = false )  
                end
                # add cut
                cut_collection[t-1].v[i][k] = c[1]
                cut_collection[t-1].π[i][k] = c[2]
            end
        end

        ## compute the lower bound
        _LB = forward_step_optimize!(StageCoefficient[1], Ω[1][1].d, [0.0 for i in 1:d], cut_collection[1])
        LB = _LB[3] + _LB[4]
        i = i + 1

        @info "LB is $LB, UB is $UB"
        if UB-LB <= ϵ * UB || i > max_iter
            return [LB,UB]
        end

    end

end


lagrangian_cut = cut_collection
































