#############################################################################################
###########################  auxiliary functions for level set method #######################
#############################################################################################

"""
This function is to constraint the model for solving gap and alpha
"""

function Δ_model_formulation(functionInfo::FunctionInfo, f_star::Float64, iter::Int64; Output::Int64 = 0)
    
    model_alpha = Model(
        optimizer_with_attributes(
                                    ()->Gurobi.Optimizer(GRB_ENV),
                                    "OutputFlag" => Output, 
                                    "Threads" => 0)
                                )

    @variable(model_alpha, z)
    @variable(model_alpha, 0 ≤ α ≤ 1)
    @constraint(model_alpha, con[j = 1:iter], z ≤  α * ( functionInfo.f_his[j] - f_star) + (1 - α) * functionInfo.G_max_his[j] )
    
    # we first compute gap Δ
    @objective(model_alpha, Max, z)
    optimize!(model_alpha)
    st = termination_status(model_alpha)
    Δ = JuMP.objective_value(model_alpha)

    
    ## then we modify above model to compute alpha
    # alpha_min
    @constraint(model_alpha, z .≥ 0)
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
function add_constraint(functionInfo::FunctionInfo, modelInfo::ModelInfo, iter::Int64)
    m = length(functionInfo.G)

    xⱼ = functionInfo.x_his[iter]
    # add constraints
    @constraint(modelInfo.model, modelInfo.z .≥ functionInfo.f_his[iter] + functionInfo.df' * (modelInfo.x - xⱼ) )
    @constraint(modelInfo.model, [k = 1:m], modelInfo.y .≥ functionInfo.G[k] + functionInfo.dG[k]' * (modelInfo.x - xⱼ) )
end



#############################################################################################
##########################    auxiliary functions for backward   ############################
#############################################################################################

function add_generator_constraint(StageProblemData::StageData, modelInfo::BackwardModelInfo;
                                                        binaryInfo::BinaryInfo = binaryInfo)

    @constraint(modelInfo.model, binaryInfo.A * modelInfo.Lc + modelInfo.x .≤ StageProblemData.ū )  ## no more than max num of generators
    @constraint(modelInfo.model, sum(modelInfo.y) + modelInfo.slack .≥ modelInfo.demand )  # satisfy demand
    @constraint(modelInfo.model, StageProblemData.h * StageProblemData.N 
                            * (binaryInfo.A * modelInfo.Lc + modelInfo.x + StageProblemData.s₀ ) .≥ modelInfo.y )  # no more than capacity

end




function add_generator_cut(cutCoefficient::CutCoefficient, modelInfo::BackwardModelInfo)


    iter = length(keys(cutCoefficient.v))  ## iter num
    k = length(keys(cutCoefficient.v[1]))  ## scenario num

    @constraint(modelInfo.model, cut[i in 1:iter-1, m in 1:k], modelInfo.θ ≥ cutCoefficient.v[i][m] + 
                                                cutCoefficient.π[i][m]' * modelInfo.Lt )
                                                
end






#############################################################################################
###################################  function: backward pass ################################
#############################################################################################

"""
    This is the oracle in level set method, and it will return [F, dF]
"""
function backwardOptimize!(StageProblemData::StageData, demand::Vector{Float64}, sum_generator::Vector{Float64}, π::Vector{Float64}, cutCoefficient::CutCoefficient; 
                        θ_bound::Real = 0.0, Enhanced_Cut::Bool = true,
                        binaryInfo::BinaryInfo = binaryInfo )

    (A, n, d) = (binaryInfo.A, binaryInfo.n, binaryInfo.d)

    F = Model(
        optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                    "OutputFlag" => 0,
                                    "Threads" => 0 )
            )

    @variable(F, x[i = 1:d] ≥ 0, Int)   # the number of generators will be built in this stage
    @variable(F, 0 ≤ Lc[i = 1:n]≤ 1)  # auxiliary variable (copy variable)
    @variable(F, Lt[i = 1:n], Bin)      # stage variable, A * Lt is total number of generators built after this stage
    @variable(F, y[i = 1:d] ≥ 0)
    @variable(F, slack ≥ 0 )
    @variable(F, θ ≥ θ_bound)

    modelInfo = BackwardModelInfo(F, x, Lt, Lc, y, θ, demand, slack, sum_generator)
    add_generator_constraint(StageProblemData, modelInfo, binaryInfo = binaryInfo)
    add_generator_cut(cutCoefficient, modelInfo)

    @constraint(F, A * sum_generator + x .== A * Lt )  ## to ensure pass a binary variable

    if Enhanced_Cut 
        @objective(F, Min, StageProblemData.c1' * x + StageProblemData.c2' * y + θ + StageProblemData.penalty * slack + 
                                                            π' * (sum_generator .- Lc) )
        optimize!(F)
        result = [ JuMP.objective_value(F), sum_generator .- round.(JuMP.value.(Lc)) ]
    else
        @objective(F, Min, StageProblemData.c1' * x + StageProblemData.c2' * y + θ + StageProblemData.penalty * slack - 
                                                            π' * Lc )
        optimize!(F)
        result = [ JuMP.objective_value(F), - round.(JuMP.value.(Lc)) ]
    end

    return result
end




function LevelSetMethod_optimization!(StageProblemData::StageData, demand::Vector{Float64}, sum_generator::Vector{Float64}, cutCoefficient::CutCoefficient; 
                                        levelSetMethodParam::LevelSetMethodParam = levelSetMethodParam, 
                                        ϵ::Float64 = 1e-4, interior_value::Float64 = 0.5, Enhanced_Cut::Bool = true, binaryInfo::BinaryInfo = binaryInfo)
    
    ######################################################################################################################
    ###############################   auxiliary function for function information   ######################################
    ######################################################################################################################
    ##  μ larger is better
    (μ, λ, threshold, nxt_bound, max_iter, Output, Output_Gap, Adj) = (levelSetMethodParam.μ, levelSetMethodParam.λ, levelSetMethodParam.threshold, levelSetMethodParam.nxt_bound, levelSetMethodParam.max_iter, levelSetMethodParam.Output,levelSetMethodParam.Output_Gap, levelSetMethodParam.Adj)
    (A, n, d) = (binaryInfo.A, binaryInfo.n, binaryInfo.d)
    l_interior = [.8 for i in 1:n] - .5 * sum_generator
    # @info "$l_interior"
    # l_interior= [interior_value for i in 1:n]

    f_star = forward_step_optimize!(StageProblemData, demand, sum_generator, cutCoefficient, binaryInfo = binaryInfo)
    f_star_value = f_star[3] + f_star[4]


    ## collect the information from the objecive f, and constraints G
    function compute_f_G(π::Vector{Float64}; 
                                            Enhanced_Cut::Bool = true, 
                                            f_star_value::Float64 = f_star_value, 
                                            StageProblemData::StageData = StageProblemData, 
                                            demand::Vector{Float64} = demand, 
                                            sum_generator::Vector{Float64} = sum_generator, 
                                            cutCoefficient::CutCoefficient = cutCoefficient, 
                                            ϵ::Float64 = ϵ )

        F_solution = backwardOptimize!(StageProblemData, demand, sum_generator, π, cutCoefficient, Enhanced_Cut = Enhanced_Cut, binaryInfo = binaryInfo)

        if Enhanced_Cut
            function_value_info  = Dict(1 => - F_solution[1] - π' * (l_interior .- sum_generator), ## obj function value
                                        2 => - F_solution[2] - (l_interior .- sum_generator),      ## gradient of obj value
                                        3 => Dict(1 => (1 - ϵ) * f_star_value - F_solution[1]),    ## constraint value
                                        4 => Dict(1 => - F_solution[2]),                           ## gradient of constraint value
                                        )
        else
            function_value_info  = Dict(1 => - F_solution[1] - π' * sum_generator,
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
    functionInfo = FunctionInfo(   Dict(1 => x₀), 
                                    Dict(1 => max(function_value_info[3][k] for k in keys(function_value_info[3]))), 
                                    Dict(1 => function_value_info[1]), 
                                    function_value_info[2], 
                                    function_value_info[4], 
                                    function_value_info[3]
                                    )

    ## model for oracle
    oracleModel = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 0)
            )



    @variable(oracleModel, z)
    @variable(oracleModel, x[i = 1:n])
    @variable(oracleModel, y ≤ 0)

    # para_oracle_bound =  abs(α * functionInfo.f_his[1] + (1-α) * functionInfo.G_max_his[1] )
    # @variable(oracleModel, z ≥ - 10^(ceil(log10(-para_oracle_bound))))
    para_oracle_bound = abs(functionInfo.f_his[1])
    z_rhs = 5.3 * 10^(ceil(log10(para_oracle_bound)))
    @constraint(oracleModel, oracle_bound, z ≥ - z_rhs)

    @objective(oracleModel, Min, z)
    oracleInfo = ModelInfo(oracleModel, x, y, z)


    nxtModel = Model(
        optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
        "OutputFlag" => Output, 
        "Threads" => 0)
        )

    @variable(nxtModel, x1[i = 1:n])
    @variable(nxtModel, z1 ≥ - nxt_bound)
    @variable(nxtModel, y1 ≥ - nxt_bound)


    while true
        add_constraint(functionInfo, oracleInfo, iter)
        optimize!(oracleModel)

        st = termination_status(oracleModel)
        if st != MOI.OPTIMAL
            @info "oracle is infeasible"
            break
        end

        f_star = JuMP.objective_value(oracleModel)

        ## formulate alpha model

        result = Δ_model_formulation(functionInfo, f_star, iter, Output = Output)
        Δ, a_min, a_max = result[1], result[2], result[3]
        
        ## update α
        if μ/2 ≤ (α-a_min)/(a_max-a_min) .≤ 1-μ/2
            α = α
        else
            α = (a_min+a_max)/2
        end

        ## update level
        w = α * f_star
        W = minimum( α * functionInfo.f_his[j] + (1-α) * functionInfo.G_max_his[j] for j in 1:iter) 
        level = w + λ * (W - w)

        if Output_Gap
            @info "Gap is $Δ, iter num is $iter, func_val is $( - function_value_info[1]), Constraint is $(functionInfo.G_max_his[iter])"
            # @info "Constraint is $(functionInfo.G_max_his[iter])"
        end
        
        ######################################################################################################################
        #########################################     next iteration point   #################################################
        ######################################################################################################################

        # obtain the next iteration point
        if iter == 1
            @constraint(nxtModel, level_constraint, α * z1 + (1 - α) * y1 ≤ level);
        else 
            delete(nxtModel, nxtModel[:level_constraint])
            unregister(nxtModel, :level_constraint)
            @constraint(nxtModel, level_constraint, α * z1 + (1 - α) * y1 ≤ level);
        end

        @constraint(nxtModel, z1 .≥ functionInfo.f_his[iter] + functionInfo.df' * (x1 - functionInfo.x_his[iter]) )
        @constraint(nxtModel, [k in keys(functionInfo.G)], y1 .≥ functionInfo.G[k] + functionInfo.dG[k]' * (x1 - functionInfo.x_his[iter]) )
        @objective(nxtModel, Min, (x1 - functionInfo.x_his[iter])' * (x1 - functionInfo.x_his[iter]))
        optimize!(nxtModel)
        st = termination_status(nxtModel)
        if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
            x_nxt = JuMP.value.(x1)
        elseif st == MOI.NUMERICAL_ERROR ## need to figure out why this case happened and fix it
            if Enhanced_Cut
                return [ - functionInfo.f_his[iter] - functionInfo.x_his[iter]' * l_interior, 
                                                functionInfo.x_his[iter]] 
            else
                return [ - functionInfo.f_his[iter] - functionInfo.x_his[iter]' * sum_generator, 
                                                functionInfo.x_his[iter]]
            end
        else
            set_normalized_rhs( level_constraint, w + 1 * (W - w))
            optimize!(nxtModel)
            x_nxt = JuMP.value.(x1)
            # break   
        end

        ## stop rule
        if Δ < threshold || iter > max_iter 
            if Enhanced_Cut
                return [ - functionInfo.f_his[iter] - functionInfo.x_his[iter]' * l_interior, 
                                                functionInfo.x_his[iter]] 
            else
                return [ - functionInfo.f_his[iter] - functionInfo.x_his[iter]' * sum_generator, 
                                                functionInfo.x_his[iter]]
            end
        end
        ######################################################################################################################
        #####################################################    end   #######################################################
        ######################################################################################################################

        ## save the trajectory
        function_value_info = compute_f_G(x_nxt, Enhanced_Cut = Enhanced_Cut)
        iter = iter + 1
        functionInfo.x_his[iter]     = x_nxt
        functionInfo.G_max_his[iter] = max(function_value_info[3][k] for k in keys(function_value_info[3]))
        functionInfo.f_his[iter]     = function_value_info[1]
        functionInfo.df              = function_value_info[2]
        functionInfo.dG              = function_value_info[4]
        functionInfo.G               = function_value_info[3]

    end

end

