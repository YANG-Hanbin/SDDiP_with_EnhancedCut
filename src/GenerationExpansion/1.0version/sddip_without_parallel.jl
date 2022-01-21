using JuMP, Test, Statistics, StatsBase, Gurobi, Distributed, ParallelDataTransfer, Random

const GRB_ENV = Gurobi.Env()


include("data_struct.jl")
include("backward_pass.jl")
include("forward_pass.jl")
include("runtests_small3.jl")
include("gurobiTest.jl")

#############################################################################################
####################################    main function   #####################################
#############################################################################################
max_iter = 200; ϵ = 1e-2; Enhand_Cut = true


function SDDiP_algorithm(Ω::Dict{Int64,Dict{Int64,RandomVariables}}, prob::Dict{Int64,Vector{Float64}}, StageCoefficient::Dict{Int64, StageData}; 
    scenario_sequence::Dict{Int64, Dict{Int64, Any}} = scenario_sequence, ϵ::Float64 = 0.001, M::Int64 = 30, max_iter::Int64 = 200, 
    Enhand_Cut::Bool = true, binaryInfo::BinaryInfo = binaryInfo)
    ## d: x dim
    ## M: num of scenarios when doing one iteration
    
    T = length(keys(Ω))
    
    i = 1
    LB = - Inf 
    UB = Inf
    
    cut_collection = Dict{Int64, CutCoefficient}()  # here, the index is stage t

    for t in 1:T

        cut_collection[t] = CutCoefficient(
                                            Dict(1=> Dict(1=> 0.0), 2 => Dict()),   ## v
                                            Dict(1=> Dict(1=> zeros(Float64, n)), 2 => Dict())  ## π
                                          )
    end

    add_result = Dict()

    while true
        M = 100
        Sol_collection = Dict()  # to store every iteration results
        u = Vector{Float64}(undef, M)  # to compute upper bound
        
        Random.seed!(i*3)
        Scenarios = SampleScenarios(scenario_sequence, T = T, M = M)
        
        ## Forward Step
        for k in 1:M
            sum_generator= [0.0 for i in 1:n]
            for t in 1:T
                if k == 1  ## create data struct for next cut coefficients
                    cut_collection[t].v[i] = Dict{Int64, Float64}(1=> 0.0)
                    cut_collection[t].π[i] = Dict{Int64, Vector{Float64}}(1=> zeros(Float64, n))
                end
                j = scenario_sequence[Scenarios[k]][1][t]  ## realization of k-th scenario at stage t

                Sol_collection[t, k] = forward_step_optimize!(StageCoefficient[t], Ω[t][j].d, sum_generator, cut_collection[t], 
                                                                binaryInfo = binaryInfo)
                sum_generator = Sol_collection[t,k][1]
            end
            u[k] = sum(Sol_collection[t, k][4] for t in 1:T)
        end

        ## compute the upper bound
        μ̄ = mean(u)
        σ̂² = var(u)
        UB = μ̄ + 1.96 * sqrt(σ̂²/M)

        ##################################### Parallel Computation for backward step ###########################

        M = 4
        

        for t = reverse(2:T)
            for k in 1:M 
                c = [0, zeros(Float64,n)]
                for j in keys(Ω[t])
                    # @info "$t $k $j"
                    ϵ_value = 1e-5
                    λ_value = .5; threshold = 1e2; Output = 0; Output_Gap = true; Adj = false;
                    levelSetMethodParam = LevelSetMethodParam(0.95, λ_value, threshold, 1e14, 3e3, Output, Output_Gap, Adj)
                    c = c + prob[t][j] * LevelSetMethod_optimization!(StageCoefficient[t], Ω[t][j].d, Sol_collection[t-1,k][1], cut_collection[t], 
                                                                        levelSetMethodParam = levelSetMethodParam, ϵ = ϵ_value, 
                                                                        Enhand_Cut = Enhand_Cut, 
                                                                        binaryInfo = binaryInfo) 
                end
                # add cut
                cut_collection[t-1].v[i][k] = c[1]
                cut_collection[t-1].π[i][k] = c[2]
            end
        end

        ########################################################################################################


        ## compute the lower bound
        _LB = forward_step_optimize!(StageCoefficient[1], Ω[1][1].d, [0.0 for i in 1:n], cut_collection[1], binaryInfo = binaryInfo)
        LB = _LB[3] + _LB[4]
        add_result[i] = _LB  # store the result
        i = i + 1

        @info "iter num is $(i-1), LB is $LB, UB is $UB, first stage decision is $(A*_LB[1])"
        if UB-LB <= ϵ * UB || i > max_iter
            return [LB,UB]
        end

    end

end



using JLD2, FileIO
result_enhanced = copy(add_result)
cut_enhanced = copy(cut_collection)

@save "runtests_small2_enhanced.jld2" result_enhanced cut_enhanced
# @load "runtests_small2_enhanced.jld2" result_enhanced cut_enhanced




result_LC = copy(add_result)
cut_LC = copy(cut_collection)

@save "runtests_small2_LC.jld2" result_LC cut_LC
# @load "runtests_small2_LC.jld2" result_LC cut_LC























