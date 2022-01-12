function LevelSetMethod_generator_v(StageProblemData::StageData, demand::Vector{Float64}, sum_generator::Vector{Float64}, cut_coefficient::CutCoefficient; 
    μ::Float64 = 0.95, λ::Float64 = 0.3, Output::Int64 = 0, threshold::Float64 = 1e3, ϵ::Float64 = 1e-2, interior_value ::Float64 = 0.5,
    d::Int64 = 6, n::Int64 = 21, max_iter::Int64 = 3e3, Enhand_Cut::Bool = true, nxt_bound::Float64 = 1e13, oracle_bound::Float64 = 1e14,
    A::Matrix{Int64} = [1 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                        0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0 0 0 0 1 2 4 8 16 32 0 0 0;
                        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 4;]
                        )

    λ_value = .5
    ϵ_value = 1e-3
    # compute average cut coefficient
    demand = Ω[t][j].d
    if t == 5
        oracle_bound = 1e13
        nxt_bound = 1e14
        if demand[1] > 3.1e8
            λ = .9
            ϵ = 5e-2  
        elseif 2.9e8 <= demand[1] <= 3.1e8
            λ = .7
            ϵ = 1e-1
        elseif 2.9e8 > demand[1]
            λ = .3
            ϵ = 1e-1
        end
    end

    ######################################################################################################################
    ###############################   auxiliary function for function information   ######################################
    ######################################################################################################################
    l_interior= [interior_value for i in 1:n]

    f_star = forward_step_optimize!(StageCoefficient[t], Ω[t][j].d, Sol_collection[t-1,k][1], cut_collection[t])
    f_star_value = f_star[3] + f_star[4]

    function compute_f_G(π::Vector{Float64}; Enhand_Cut::Bool = true )
        F_solution = backward_step_F(StageCoefficient[t], Ω[t][j].d, Sol_collection[t-1,k][1], π, cut_collection[t])
        if Enhand_Cut
            function_value_info  = Dict(1 => -F_solution[1] - π' * (A * l_interior .- Sol_collection[t-1,k][1]),
                                        2 => -F_solution[2] - (A * l_interior .- Sol_collection[t-1,k][1]),
                                        3 => Dict(1 => (1- ϵ) * f_star_value - F_solution[1]),
                                        4 => Dict(1 => - F_solution[2]),
                                        )
        else
            function_value_info  = Dict(1 => - F_solution[1] - π' *  Sol_collection[t-1,k][1],
                                        2 => - F_solution[2] - Sol_collection[t-1,k][1],
                                        3 => Dict(1 => 0.0 ),
                                        4 => Dict(1 => - F_solution[2] * 0),
                                        )
        end
        return function_value_info
        ## Com_f = function_value_info[1], Com_grad_f = function_value_info[2], 
        ## Com_G = function_value_info[3], Com_grad_G = function_value_info[4], Com_max_g = function_value_info[3]
    end

    ######################################################################################################################
    ##############################################   level set method   ##################################################
    ######################################################################################################################
    
    x₀ = ones(d)

    iter = 1
    α = 1/2
    
    ## trajectory
    function_value_info = compute_f_G(x₀, Enhand_Cut = Enhand_Cut)
    function_info = FunctionInfo(   Dict(1 => x₀), 
                                    function_value_info[3], 
                                    Dict(1 => function_value_info[1]), 
                                    function_value_info[2], 
                                    function_value_info[4], 
                                    function_value_info[3]
                                    )

    ## model for oracle
    model_oracle = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 1)
            )

    # if Enhand_Cut
    #     @variable(model_oracle, z >= - oracle_bound)
    # else 
    #     para_oracle_bound =  α * function_info.f_his[1] + (1-α) * function_info.G_max_his[1] 
    #     @variable(model_oracle, z >= - 10^(ceil(log10(-para_oracle_bound))))
    # end
    @variable(model_oracle, z)
    @variable(model_oracle, x[i = 1:d])
    @variable(model_oracle, y <= 0)

    para_oracle_bound =  abs(function_info.f_his[1])
    z_rhs = 5 * 10^(ceil(log10(para_oracle_bound)))
    @constraint(model_oracle, oracle_bound, z >= - z_rhs)

    @objective(model_oracle, Min, z)
    oracle_info = ModelInfo(model_oracle, x, y, z)


    while true
        if true
            param_z_rhs = abs(function_info.f_his[iter])
            if z_rhs <  param_z_rhs
                z_rhs = 1.2 * z_rhs
            end

            if z_rhs > 2 * param_z_rhs
                z_rhs = 0.8 * z_rhs
            end 
            set_normalized_rhs(oracle_bound, - z_rhs)  
        end
        add_constraint(function_info, oracle_info, iter)
        optimize!(model_oracle)

        st = termination_status(model_oracle)
        if st != MOI.OPTIMAL
            @info "oracle is infeasible"
            # break
        end

        f_star = JuMP.objective_value(model_oracle)

        ## formulate alpha model

        result = Δ_model_formulation(function_info, f_star, iter, Output = Output)
        Δ, a_min, a_max = result[1], result[2], result[3]
        @info "Gap is $Δ, iter num is $iter"
        ## update α
        if μ/2 <= (α-a_min)/(a_max-a_min) .<= 1-μ/2
            α = α
        else
            α = (a_min+a_max)/2
        end

        ## update level
        w = α * f_star
        W = minimum( α * function_info.f_his[j] + (1-α) * function_info.G_max_his[j] for j in 1:iter) 
        # @info "w is $w, W is $W"
        level = w + λ * (W - w)
        
        ######################################################################################################################
        #########################################     next iteration point   #################################################
        ######################################################################################################################

        ## obtain the next iteration point
        model_nxt = Model(
            optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 1)
            )

        @variable(model_nxt, x1[i = 1:d])
        @variable(model_nxt, z1 >= - nxt_bound)
        @variable(model_nxt, y1)

        @constraint(model_nxt, level_constraint, α * z1 + (1 - α) * y1 <= level)
        @constraint(model_nxt, z1 .>= function_info.f_his[iter] + function_info.df' * (x1 - function_info.x_his[iter]) )
        @constraint(model_nxt, [k in keys(function_info.G)], y1 .>= function_info.G[k] + function_info.dG[k]' * (x1 - function_info.x_his[iter]) )
        @objective(model_nxt, Min, (x1 - function_info.x_his[iter])' * (x1 - function_info.x_his[iter]))
        optimize!(model_nxt)
        st = termination_status(model_nxt)
        if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
            x_nxt = JuMP.value.(x1)
        elseif st == MOI.NUMERICAL_ERROR ## need to figure out why this case happened and fix it
            if Enhand_Cut
                return [ - function_value_info[1], function_info.x_his[iter]] 
            else
                return [ - function_value_info[1] - function_info.x_his[iter]' * Sol_collection[t-1,k][1], 
                                                function_info.x_his[iter]]
            end
        else
            set_normalized_rhs( level_constraint, w + 0.95 * (W - w))
            optimize!(model_nxt)
            x_nxt = JuMP.value.(x1)
            # break   
        end

        ## stop rule
        if Δ < threshold || iter > max_iter
            if Enhand_Cut
                return [ - function_value_info[1], function_info.x_his[iter]] 
            else
                return [ - function_value_info[1] - function_info.x_his[iter]' * Sol_collection[t-1,k][1], 
                                                function_info.x_his[iter]]
            end
        end
        ######################################################################################################################
        #####################################################    end   #######################################################
        ######################################################################################################################

        ## save the trajectory
        function_value_info = compute_f_G(x_nxt)
        iter = iter + 1
        function_info.x_his[iter]     = x_nxt
        function_info.G_max_his[iter] = function_value_info[3][1]
        function_info.f_his[iter]     = function_value_info[1]
        function_info.df              = function_value_info[2]
        function_info.dG              = function_value_info[4]
        function_info.G               = function_value_info[3]
        # @info "Gap is $Δ, iter num is $iter"

    end

end