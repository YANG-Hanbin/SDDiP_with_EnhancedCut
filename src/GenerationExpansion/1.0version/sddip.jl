using Distributed
addprocs(5)

@everywhere using JuMP, Test, Statistics, StatsBase, Gurobi, Distributed, ParallelDataTransfer, Random

@everywhere const GRB_ENV = Gurobi.Env()

@everywhere begin
    include("data_struct.jl")
    include("backward_pass.jl")
    include("forward_pass.jl")
    include("data_file.jl")
end
#############################################################################################
####################################    main function   #####################################
#############################################################################################
# n = 21; d = 6; max_iter = 20; ϵ = 1e-2; Enhand_Cut = true

# @broadcast n = 21; d = 6; max_iter = 20; ϵ = 1e-2; Enhand_Cut = true;
# @broadcast Enhand_Cut = true;

function SDDiP_algorithm(Ω::Dict{Int64,Dict{Int64,RandomVariables}}, prob::Dict{Int64,Vector{Float64}}, StageCoefficient::Dict{Int64, StageData}; 
    scenario_sequence::Dict{Int64, Dict{Int64, Any}} = scenario_sequence, ϵ::Float64 = 0.001, M::Int64 = 30, max_iter::Int64 = 20, 
    n::Int64 = 2, d::Int64 = 1, A::Matrix{Int64} = [1 2;], Enhand_Cut::Bool = true)
    ## d: x dim
    ## M: num of scenarios when doing one iteration
    
    @broadcast T = length(keys(Ω))
    
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


    ## an auxiliary function for backward iteration
    @everywhere begin 
        function inner_func_backward(j, k, t)
            @info "$t $k $j"
            λ_value = .1
            ϵ_value = 5e-4
            nxt_bound = 1e14
            # compute average cut coefficient
            demand = Ω[t][j].d
            if t == 5
                if demand[1] > 3.0e8
                    λ_value = .3
                    ϵ_value = 1e-2
                elseif 2.8e8 <= demand[1] <= 3.0e8
                    λ_value = .2
                    ϵ_value = 1e-2
                elseif 2.8e8 > demand[1]
                    λ_value = .01
                    ϵ_value = 1e-2
                end
            end
            @time c =  prob[t][j] * LevelSetMethod_optimization!(StageCoefficient[t], 
                                                                Ω[t][j].d, 
                                                                Sol_collection[t-1,k][1], 
                                                                cut_collection[t], 
                                                                max_iter = 3000, 
                                                                threshold = 1e4,
                                                                Enhand_Cut = Enhand_Cut,  
                                                                nxt_bound = nxt_bound,
                                                                μ = 0.95, λ = λ_value, ϵ = ϵ_value,
                                                                Output_Gap = false, Adj = Adj, threshold = 1e-3, 
                                                                A = A, d = d, n = n ) 
        end
    end


    while true
        M = 30
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

        @passobj 1 workers() Sol_collection

        ## compute the upper bound
        μ̄ = mean(u)
        σ̂² = var(u)
        UB = μ̄ + 1.96 * sqrt(σ̂²/M)

        ##################################### Parallel Computation for backward step ###########################

        M = 4

        for t = reverse(2: T)
            cut_collection[t-1].v[i] = Dict()
            cut_collection[t-1].π[i] = Dict()
            @passobj 1 workers() cut_collection
            for k in 1:M 
                p = pmap(inner_func_backward, 1:10, [k for i in keys(Ω[t])], [t for i in keys(Ω[t])])
                c = sum(p)
                cut_collection[t-1].v[i][k] = c[1]
                cut_collection[t-1].π[i][k] = c[2]
            end
            @passobj 1 workers() cut_collection
        end


        ########################################################################################################



        ## compute the lower bound
        _LB = forward_step_optimize!(StageCoefficient[1], Ω[1][1].d, [0.0 for i in 1:n], cut_collection[1], A = A, d = d, n = n )
        LB = _LB[3] + _LB[4]
        add_result[i] = _LB  # store the result
        i = i + 1

        @info "LB is $LB, UB is $UB"
        if UB-LB <= ϵ * UB || i > max_iter
            return [LB,UB]
        end

    end

end































