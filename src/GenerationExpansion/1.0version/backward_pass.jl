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
                                                        binaryInfo::BinaryInfo = binaryInfo)

    @constraint(model_info.model, binaryInfo.A * model_info.Lc + model_info.x .<= StageProblemData.ū )  ## no more than max num of generators
    @constraint(model_info.model, sum(model_info.y) + model_info.slack .>= model_info.demand )  # satisfy demand
    @constraint(model_info.model, StageProblemData.h * StageProblemData.N 
                            * (binaryInfo.A * model_info.Lc + model_info.x + StageProblemData.s₀ ) .>= model_info.y )  # no more than capacity

end




function add_generator_cut(cut_coefficient::CutCoefficient, model_info::BackwardModelInfo)


    iter = length(keys(cut_coefficient.v))  ## iter num
    k = length(keys(cut_coefficient.v[1]))  ## scenario num

    @constraint(model_info.model, cut[i in 1:iter-1, m in 1:k], model_info.θ >= cut_coefficient.v[i][m] + 
                                                cut_coefficient.π[i][m]' * model_info.Lt )
                                                
end






#############################################################################################
###################################  function: backward pass ################################
#############################################################################################

"""
    This is the oracle in level set method, and it will return [F, dF]
"""
function backward_step_F(StageProblemData::StageData, demand::Vector{Float64}, sum_generator::Vector{Float64}, π::Vector{Float64}, cut_coefficient::CutCoefficient; 
                        θ_bound::Real = 0.0, Enhanced_Cut::Bool = true,
                        binaryInfo::BinaryInfo = binaryInfo )

    (A, n, d) = (binaryInfo.A, binaryInfo.n, binaryInfo.d)

    F = Model(
        optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                    "OutputFlag" => 0,
                                    "Threads"=>1)
            )

    @variable(F, x[i = 1:d] >=0, Int)   # the number of generators will be built in this stage
    @variable(F, 0 <= Lc[i = 1:n]<= 1)  # auxiliary variable (copy variable)
    @variable(F, Lt[i = 1:n], Bin)      # stage variable, A * Lt is total number of generators built after this stage
    @variable(F, y[i = 1:d] >= 0)
    @variable(F, slack >= 0 )
    @variable(F, θ >= θ_bound)

    model_info = BackwardModelInfo(F, x, Lt, Lc, y, θ, demand, slack, sum_generator)
    add_generator_constraint(StageProblemData, model_info, binaryInfo = binaryInfo)
    add_generator_cut(cut_coefficient, model_info)

    @constraint(F, A * sum_generator + x .== A * Lt )  ## to ensure pass a binary variable

    if Enhanced_Cut 
        @objective(F, Min, StageProblemData.c1' * x + StageProblemData.c2' * y + θ + StageProblemData.penalty * slack + 
                                                            π' * (sum_generator .- Lc) )
        optimize!(F)
        result = [ JuMP.objective_value(F), sum_generator .- JuMP.value.(Lc) ]
    else
        @objective(F, Min, StageProblemData.c1' * x + StageProblemData.c2' * y + θ + StageProblemData.penalty * slack - 
                                                            π' * Lc )
        optimize!(F)
        result = [ JuMP.objective_value(F), - round.(JuMP.value.(Lc)) ]
    end

    return result
end



struct LevelSetMethodParam
    μ             ::Float64   ## param for adjust α
    λ             ::Float64   ## param for adjust level
    threshold     ::Float64   ## threshold for Δ
    nxt_bound     ::Float64   ## lower bound for solving next iteration point π
    max_iter      ::Int64     
    Output        ::Int64     ## Gurobi Output parameter
    Output_Gap    ::Bool      ## if True will print Δ info
    Adj           ::Bool      ## whether adjust oracle lower bound
end


# levelSetMethodParam = LevelSetMethodParam(μ, λ, threshold, nxt_bound, max_iter, Output, Output_Gap, Adj)




function LevelSetMethod_optimization!(StageProblemData::StageData, demand::Vector{Float64}, sum_generator::Vector{Float64}, cut_coefficient::CutCoefficient; 
                                        levelSetMethodParam::LevelSetMethodParam = levelSetMethodParam, 
                                        ϵ::Float64 = 1e-4, interior_value::Float64 = 0.5, Enhanced_Cut::Bool = true, binaryInfo::BinaryInfo = binaryInfo)
    
    ######################################################################################################################
    ###############################   auxiliary function for function information   ######################################
    ######################################################################################################################
    ##  μ larger is better
    (μ, λ, threshold, nxt_bound, max_iter, Output, Output_Gap, Adj) = (levelSetMethodParam.μ, levelSetMethodParam.λ, levelSetMethodParam.threshold, levelSetMethodParam.nxt_bound, levelSetMethodParam.max_iter, levelSetMethodParam.Output,levelSetMethodParam.Output_Gap, levelSetMethodParam.Adj)
    (A, n, d) = (binaryInfo.A, binaryInfo.n, binaryInfo.d)
    l_interior= [interior_value for i in 1:n]

    f_star = forward_step_optimize!(StageProblemData, demand, sum_generator, cut_coefficient, binaryInfo = binaryInfo)
    f_star_value = f_star[3] + f_star[4]


    ## collect the information from the objecive f, and constraints G
    function compute_f_G(π::Vector{Float64}; Enhanced_Cut::Bool = true, f_star_value::Float64 = f_star_value, 
        StageProblemData::StageData = StageProblemData, demand::Vector{Float64} = demand, 
        sum_generator::Vector{Float64} = sum_generator, cut_coefficient::CutCoefficient = cut_coefficient   )

        F_solution = backward_step_F(StageProblemData, demand, sum_generator, π, cut_coefficient, Enhanced_Cut = Enhanced_Cut, binaryInfo = binaryInfo)

        if Enhanced_Cut
            function_value_info  = Dict(1 => - F_solution[1] - π' * (l_interior .- sum_generator),
                                        2 => - F_solution[2] - (l_interior .- sum_generator),
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
    
    x₀ = ones(n)

    iter = 1
    α = 1/2

    ## trajectory
    function_value_info = compute_f_G(x₀, Enhanced_Cut = Enhanced_Cut)
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
        if Adj
            param_z_rhs = abs(function_info.f_his[iter])
            if z_rhs <  1.5 * param_z_rhs
                # @info "z level up $(z_rhs/param_z_rhs)"
                z_rhs = 2 * z_rhs
            end

            if z_rhs > 10 * param_z_rhs
                # @info "z level down $(z_rhs/param_z_rhs) "
                z_rhs = .1 * z_rhs
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
            if Enhanced_Cut
                return [ - function_info.f_his[iter] - function_info.x_his[iter]' * l_interior, 
                                                function_info.x_his[iter]] 
            else
                return [ - function_info.f_his[iter] - function_info.x_his[iter]' * sum_generator, 
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
            if Enhanced_Cut
                return [ - function_info.f_his[iter] - function_info.x_his[iter]' * l_interior, 
                                                function_info.x_his[iter]] 
            else
                return [ - function_info.f_his[iter] - function_info.x_his[iter]' * sum_generator, 
                                                function_info.x_his[iter]]
            end
        end
        ######################################################################################################################
        #####################################################    end   #######################################################
        ######################################################################################################################

        ## save the trajectory
        function_value_info = compute_f_G(x_nxt, Enhanced_Cut = Enhanced_Cut)
        iter = iter + 1
        function_info.x_his[iter]     = x_nxt
        function_info.G_max_his[iter] = function_value_info[3][1]
        function_info.f_his[iter]     = function_value_info[1]
        function_info.df              = function_value_info[2]
        function_info.dG              = function_value_info[4]
        function_info.G               = function_value_info[3]

    end

end

