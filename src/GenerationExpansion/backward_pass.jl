#############################################################################################
###########################  auxiliary functions for level set method #######################
#############################################################################################

"""
This function is to constraint the model for solving gap and alpha
"""

function Δ_model_formulation(function_info::FunctionInfo, f_star::Float64, iter::Int64; Output::Int64 = 0)
    
    model_alpha = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV),
            "OutputFlag" => Output, 
            "Threads" => 1)
            )

    @variable(model_alpha, z)
    @variable(model_alpha, 0 <= α <= 1)
    @constraint(model_alpha, con[j = 1:iter], z <=  α * ( function_info.f_his[j] - f_star) + (1 - α) * function_info.G_max_his[j] )
    
    # we first compute gap Δ
    @objective(model_alpha, Max, z)
    optimize!(model_alpha)
    st = termination_status(model_alpha)
    Δ = JuMP.objective_value(model_alpha)

    
    ## then we modify above model to compute alpha
    # alpha_min
    @constraint(model_alpha, z .>= 0)
    @objective(model_alpha, Min, α)
    optimize!(model_alpha)
    a_min = JuMP.value(α)

    # alpha_max
    @objective(model_alpha, Max, α)
    optimize!(model_alpha)
    a_max = JuMP.value(α)

    return Dict(1 => Δ, 2 => a_min, 3 => a_max)

end


"""
    This function is to add constraints for the model f_star and nxt pt.
"""
function add_constraint(function_info::FunctionInfo, model_info::ModelInfo, iter::Int64)
    m = length(function_info.G)

    xⱼ = function_info.x_his[iter]
    # add constraints
    @constraint(model_info.model, model_info.z .>= function_info.f_his[iter] + function_info.df' * (model_info.x - xⱼ) )
    @constraint(model_info.model, [k = 1:m], model_info.y .>= function_info.G[k] + function_info.dG[k]' * (model_info.x - xⱼ) )
end



#############################################################################################
##########################    auxiliary functions for backward   ############################
#############################################################################################

function add_generator_constraint(StageProblemData::StageData, model_info::BackwardModelInfo;
    A::Matrix{Int64} = [1 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                        0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0 0 0 0 1 2 4 8 16 32 0 0 0;
                        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 4;]
                    )
    @constraint(model_info.model, A * (model_info.l + model_info.L) .<= StageProblemData.ū )  ## no more than max num of generators
    @constraint(model_info.model, sum(model_info.y) + model_info.slack .>= model_info.demand )  # satisfy demand
    @constraint(model_info.model, StageProblemData.h * StageProblemData.N 
                            * (A * (model_info.l + model_info.L) + StageProblemData.s₀ ) .>= model_info.y )  # no more than capacity

end




function add_generator_cut(cut_coefficient::CutCoefficient, model_info::BackwardModelInfo; 
    Enhand_Cut::Bool = true,
    d::Int64 = 8, n::Int64 = 21, interior_value ::Float64 = 0.5,
    A::Matrix{Int64} = [1 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                        0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0 0 0 0 1 2 4 8 16 32 0 0 0;
                        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 4;]
                        )

    l_interior= [interior_value for i in 1:n]

    iter = length(keys(cut_coefficient.v))  ## iter num
    k = length(keys(cut_coefficient.v[1]))  ## scenario num
    if Enhand_Cut 
        @constraint(model_info.model, cut[i in 1:iter-1, m in 1:k], model_info.θ >= cut_coefficient.v[i][m] + 
                                                cut_coefficient.π[i][m]' * A * (model_info.l .- l_interior))
    else
        @constraint(model_info.model, cut[i in 1:iter-1, m in 1:k], model_info.θ >= cut_coefficient.v[i][m] + 
                                                cut_coefficient.π[i][m]' * A * model_info.l )
    end

end






#############################################################################################
###################################  function: backward pass ################################
#############################################################################################


function backward_step_F(StageProblemData::StageData, demand::Vector{Float64}, sum_generator::Vector{Float64}, π::Vector{Float64}, cut_coefficient::CutCoefficient; 
                        θ_bound::Real = 0.0, interior_value ::Float64 = 0.5, Enhand_Cut::Bool = true,
                        d::Int64 = 6, n::Int64 = 21, 
                        A::Matrix{Int64} = [1 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                            0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                            0 0 0 0 0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0;
                                            0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
                                            0 0 0 0 0 0 0 0 0 0 0 0 1 2 4 8 16 32 0 0 0;
                                            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 4;]
                        )
    F = Model(
        optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                    "OutputFlag" => 0,
                                    "Threads"=>1)
            )

    @variable(F, l[i = 1:n], Bin)      # current stage variable l_t
    @variable(F, 0 <= L[i = 1:n]<= 1)  # auxiliary variable (copy variable)
    @variable(F, y[i = 1:d] >= 0)
    @variable(F, slack >= 0 )
    @variable(F, θ >= θ_bound)

    model_info = BackwardModelInfo(F, l, L, y, θ, demand, slack, sum_generator)
    add_generator_constraint(StageProblemData, model_info)
    add_generator_cut(cut_coefficient, model_info, Enhand_Cut = Enhand_Cut)

    if Enhand_Cut 
        @objective(F, Min, StageProblemData.c1' * A * l + StageProblemData.c2' * y +  θ + StageProblemData.penalty * slack + 
                                                            π' * (sum_generator .- A * L) )
        optimize!(F)
        result = [ JuMP.objective_value(F), sum_generator .- A * JuMP.value.(L) ]
    else
        @objective(F, Min, StageProblemData.c1' * A * l + StageProblemData.c2' * y +  θ + StageProblemData.penalty * slack - 
                                                            π' * A * L )
        optimize!(F)
        result = [ JuMP.objective_value(F), - A * JuMP.value.(L) ]
    end
    return result
end



# function backward_step_f_star(StageProblemData::StageData, demand::Vector{Float64}, sum_generator::Vector{Float64}, cut_coefficient::CutCoefficient; 
#     θ_bound::Real = 0.0, interior_value ::Float64 = 0.5,
#     d::Int64 = 6, n::Int64 = 21, 
#     A::Matrix{Int64} = [1 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
#                         0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
#                         0 0 0 0 0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0;
#                         0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
#                         0 0 0 0 0 0 0 0 0 0 0 0 1 2 4 8 16 32 0 0 0;
#                         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 4;]
#     )
#     f_star = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
#                                                 "OutputFlag" => 0, 
#                                                 "Threads" => 1) 
#                     )

#     @variable(f_star, l[i = 1:n], Bin)
#     @variable(f_star, y[i = 1:d] >= 0)
#     @variable(f_star, slack >= 0 )
#     @variable(f_star, θ >= θ_bound)
#     model_info = ForwardModelInfo(f_star, l, y, θ, demand, slack, sum_generator)

#     add_generator_constraint(StageProblemData, model_info)
#     add_generator_cut(cut_coefficient, model_info)

#     @objective(f_star, Min, StageProblemData.c1'* A * l + StageProblemData.c2' * y + StageProblemData.penalty * slack +  θ)

#     optimize!(f_star)

#     return JuMP.objective_value(f_star)
# end


##  μ larger is better, λ is flexibile
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
    ######################################################################################################################
    ###############################   auxiliary function for function information   ######################################
    ######################################################################################################################
    l_interior= [interior_value for i in 1:n]

    f_star = forward_step_optimize!(StageProblemData, demand, sum_generator, cut_coefficient)
    f_star_value = f_star[3] + f_star[4]  ## here, we obtain the optimal value f^*

    """
        compute_f_G
        This function is designed to collect the info of objective and constraints, such as their values and gradients,
            and return a Dict
    """
    function compute_f_G(π::Vector{Float64}; Enhand_Cut::Bool = true )
        F_solution = backward_step_F(StageProblemData, demand, sum_generator, π, cut_coefficient)
        ##  minus F_solution because this level set method is designed to obtain minimum
        if Enhand_Cut
            function_value_info  = Dict(1 => -F_solution[1] - π' * (A * l_interior .- sum_generator),  
                                        2 => -F_solution[2] - (A * l_interior .- sum_generator),
                                        3 => Dict(1 => (1- ϵ) * f_star_value - F_solution[1]),
                                        4 => Dict(1 => - F_solution[2]),
                                        )
        else
            function_value_info  = Dict(1 => - F_solution[1] - π' *  sum_generator,
                                        2 => - F_solution[2] - sum_generator,
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
<<<<<<< HEAD
    para_oracle_bound =  α * function_info.f_his[1] + (1-α) * function_info.G_max_his[1] 
    @variable(model_oracle, z >= - 10^(ceil(log10(-para_oracle_bound))))  ## the lower bound is crucial

=======
    # para_oracle_bound =  α * function_info.f_his[1] + (1-α) * function_info.G_max_his[1] 
    # @variable(model_oracle, z >= - 10^(ceil(log10(-para_oracle_bound))))
    @variable(model_oracle, z)
>>>>>>> dev
    @variable(model_oracle, x[i = 1:d])
    @variable(model_oracle, y <= 0)

    para_oracle_bound =  abs(α * function_info.f_his[1] + (1-α) * function_info.G_max_his[1])
    z_rhs = 3 * 10^(ceil(log10(para_oracle_bound)))
    @constraint(model_oracle, oracle_bound, z >= - z_rhs)

    @objective(model_oracle, Min, z)
    oracle_info = ModelInfo(model_oracle, x, y, z)



    while true
        if true
            param_z_rhs = abs(function_info.f_his[iter])
            if z_rhs <  param_z_rhs
                z_rhs = 1.3 * z_rhs
            end

            if z_rhs > 2 * param_z_rhs
                z_rhs = 0.7 * z_rhs
            end 
            set_normalized_rhs(oracle_bound, - z_rhs)  
        end
        add_constraint(function_info, oracle_info, iter)
        optimize!(model_oracle)

        st = termination_status(model_oracle)
        if st != MOI.OPTIMAL
            @info "oracle is infeasible"
            break
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
                return [ - function_value_info[1] - function_info.x_his[iter]' * sum_generator, 
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
                return [ - function_value_info[1] - function_info.x_his[iter]' * sum_generator,  ## we only return F or L_n (see sddip)
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

