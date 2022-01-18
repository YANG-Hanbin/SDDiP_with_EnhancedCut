LevelSetMethod_optimization!(StageCoefficient[t], 
                                                                                    Ω[t][j].d, 
                                                                                    Sol_collection[t-1,k][1], 
                                                                                    cut_collection[t], 
                                                                                    max_iter = 4000, 
                                                                                    Enhand_Cut = Enhand_Cut,  
                                                                                    nxt_bound = nxt_bound,
                                                                                    μ = 0.95, λ = λ_value, ϵ = ϵ_value,
                                                                                    Output_Gap = true, Adj = Adj ) 


function LevelSetMethod_optimization!(StageProblemData::StageData, demand::Vector{Float64}, sum_generator::Vector{Float64}, cut_coefficient::CutCoefficient; 
                                        μ::Float64 = 0.95, λ::Float64 = 0.3, ϵ::Float64 = 1e-2, interior_value::Float64 = 0.5,
                                        d::Int64 = 1, n::Int64 = 2, threshold::Float64 = 1e3, nxt_bound::Float64 = 1e13, Adj::Bool = true, 
                                        max_iter::Int64 = 3e3, Enhand_Cut::Bool = true, Output_Gap::Bool = false, Output::Int64 = 0)
    ######################################################################################################################
    μ= 0.95;ϵ=1e-5;λ = .3;interior_value =.5;threshold=1e1;Output = 0;Output_Gap = true;Enhand_Cut = true;max_iter = 3e3;Adj = false;
    ###############################   auxiliary function for function information   ######################################
    ######################################################################################################################
    l_interior= [interior_value for i in 1:n]

    f_star = forward_step_optimize!(StageCoefficient[t], Ω[t][j].d, Sol_collection[t-1,k][1], cut_collection[t])
    f_star_value = f_star[3] + f_star[4]


    ## collect the information from the objecive f, and constraints G
    function compute_f_G(π::Vector{Float64}; Enhand_Cut::Bool = true, f_star_value::Float64 = f_star_value, 
                            StageProblemData::StageData = StageCoefficient[t], demand::Vector{Float64} = Ω[t][j].d, 
                            sum_generator::Vector{Float64} = Sol_collection[t-1,k][1], cut_coefficient::CutCoefficient = cut_collection[t]   )

        F_solution = backward_step_F(StageCoefficient[t], Ω[t][j].d, Sol_collection[t-1,k][1], π, cut_collection[t], Enhand_Cut = Enhand_Cut)

        if Enhand_Cut
            function_value_info  = Dict(1 => - F_solution[1] - π' * (l_interior .- Sol_collection[t-1,k][1]),
                2 => - F_solution[2] - (l_interior .- Sol_collection[t-1,k][1]),
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

    x₀ = ones(n)

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



    @variable(model_oracle, z)
    @variable(model_oracle, x[i = 1:n])
    @variable(model_oracle, y <= 0)

    # para_oracle_bound =  abs(α * function_info.f_his[1] + (1-α) * function_info.G_max_his[1] )
    # @variable(model_oracle, z >= - 10^(ceil(log10(-para_oracle_bound))))
    para_oracle_bound = abs(function_info.f_his[1])
    z_rhs = 5.3 * 10^(ceil(log10(para_oracle_bound)))
    @constraint(model_oracle, oracle_bound, z >= - z_rhs)

    @objective(model_oracle, Min, z)
    oracle_info = ModelInfo(model_oracle, x, y, z)



    while true
        add_constraint(function_info, oracle_info, iter)
        @info "$model_oracle"
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

        ## update α
        if μ/2 <= (α-a_min)/(a_max-a_min) .<= 1-μ/2
            α = α
        else
            α = (a_min+a_max)/2
        end

        ## update level
        w = α * f_star
        W = minimum( α * function_info.f_his[j] + (1-α) * function_info.G_max_his[j] for j in 1:iter) 
        level = w + λ * (W - w)

        if Output_Gap == true
            @info "Gap is $Δ, iter num is $iter, func_val is $( - function_value_info[1]), alpha is $α, w is $w, W is $W"
            @info "Constraint is $(function_info.G_max_his[iter])"
        end

        ######################################################################################################################
        #########################################     next iteration point   #################################################
        ######################################################################################################################

        ## obtain the next iteration point
            model_nxt = Model(
            optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 1)
            )

            @variable(model_nxt, x1[i = 1:n])
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
                    return [ - function_info.f_his[iter], function_info.x_his[iter]] 
                else
                    return [ - function_info.f_his[iter] - function_info.x_his[iter]' * Sol_collection[t-1,k][1], 
                            function_info.x_his[iter]]
                end
            else
                set_normalized_rhs( level_constraint, w + 1 * (W - w))
                optimize!(model_nxt)
                x_nxt = JuMP.value.(x1)
                # break   
            end

            ## stop rule
            if Δ < threshold || iter > max_iter 
                if Enhand_Cut
                    return [ - function_info.f_his[iter], function_info.x_his[iter]] 
                else
                    return [ - function_info.f_his[iter] - function_info.x_his[iter]' * Sol_collection[t-1,k][1], 
                            function_info.x_his[iter]]
                end
            end
            ######################################################################################################################
            #####################################################    end   #######################################################
            ######################################################################################################################

            ## save the trajectory
            function_value_info = compute_f_G(x_nxt, Enhand_Cut = Enhand_Cut)
            iter = iter + 1
            function_info.x_his[iter]     = x_nxt
            function_info.G_max_his[iter] = function_value_info[3][1]
            function_info.f_his[iter]     = function_value_info[1]
            function_info.df              = function_value_info[2]
            function_info.dG              = function_value_info[4]
            function_info.G               = function_value_info[3]

    end

end