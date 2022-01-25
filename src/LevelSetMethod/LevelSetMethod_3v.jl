using JuMP, Test, Statistics, StatsBase, Gurobi, Distributed, Random, Flux

const GRB_ENV = Gurobi.Env()

"""
This algorithm is aimed to minimize a convex optimization
        min f(x)
       s.t. G(x) <= 0  ## total m constraints
            x ∈ Q
"""
#############################################################################################
####################################   Data Structure   #####################################
#############################################################################################

mutable struct FunctionInfo
    x_his        :: Dict{Int64, Vector{Float64}}  ## record every x_j point
    G_max_his    :: Dict{Int64, Float64}          ## record max(g[k] for k in 1:m)(x_j)
    f_his        :: Dict{Int64, Float64}          ## record f(x_j)
    df           :: Vector{Float64}
    dG           :: Dict{Int64, Vector{Float64}}  ## actually is a matrix.  But we use dict to store it
    G            :: Dict{Int64, Float64}          
end




struct ModelInfo
    model :: Model
    x     :: Vector{VariableRef}
    y     :: VariableRef
    z     :: VariableRef
end


#############################################################################################
###############################    auxiliary functions   ####################################
#############################################################################################

"""
This function is to constraint the model for solving gap and alpha
"""

function Δ_model_formulation(function_info::FunctionInfo, f_star::Float64, iter::Int64; Output::Int64 = 1)
    
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
###############################   function information   ####################################
#############################################################################################
"""
    This function is to evaluate the obj, gradient of obj, constraints, gradient of constraint.

        Need to note is that we will return a Dict,
        f     :: Float64 
        df    ::Vector{Float64}
        G     ::Dict{Int64, Float64}  
        dG    ::Dict{Int64, Vector{Float64}}
        G_max :: Float64


"""

# function evaluate_function(x::Vector{Float64}; M::Float64 = 1e3, x̂::Float64 = 8.0, Output::Int64 = 0)
#     ## here, x : Lagrangian Multiplier
#     ## here, x̂ : the first stage variable
#     model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
#                                             "OutputFlag" => Output, 
#                                             "Threads" => 1) 
# 					)

#     @variable(model, z[i = 1:4] >= 0)
#     @variable(model, y[i = 1:4], Bin)
#     @variable(model, x̄ >= 0)
#     @variable(model, x_c >= 0)  # local copy
#     # @constraint(model, x_c == x̂)


#     @constraint(model, [i = 1:4], z[i] - x_c<= 0 )
#     @constraint(model, [i = 1:4], x_c - (1 - y[i]) * M - z[i]<= 0)
#     @constraint(model, [i = 1:4], z[i] <= M * y[i])
#     x_lower = [0, 3, 7, 12]
#     x_upper = [3, 7, 12, 15]
#     @constraint(model, [i = 1:4], x_c - x_upper[i] - M * (1 - y[i]) <= 0)
#     @constraint(model, [i = 1:4], x_lower[i] - M * (1 - y[i]) - x_c <= 0 )
#     @constraint(model, sum(y) == 1)

#     z_coef = [-1, 1, -1, 1/3]
#     y_coef = [9, 6, 19, 7]
#     # @objective(model, Min, z_coef' * z + y_coef' * y )
#     @objective(model, Min, z_coef' * z + y_coef' * y + x[1] * (x̂ - x_c) )
#     optimize!(model)

#     # evaluate obj
#     f = JuMP.objective_value(model)

#     # evaluate constraints and their gradients
#     d1 = 0.0
#     # d2 = [3,9]' * x - 15
#     Com_G = Dict{Int64, Float64}()
#     Com_G[1] = d1

#     d1 = [0.0, ]
#     # d2 = [3.0, 9.0]
#     Com_grad_G = Dict{Int64, Vector{Float64}}()
#     Com_grad_G[1] = d1
#     # Com_grad_G[2] = d2

#     result = Dict()
#     result[1] = -f
#     result[2] = [JuMP.value(x_c) - JuMP.value(x̂), ]
#     result[3] = Com_G
#     result[4] = Com_grad_G
#     result[5] = maximum(Com_G[i] for i in keys(Com_G))

#     return result
# end








function evaluate_function(x::Vector{Float64}; M::Float64 = 1e3, x̂::Vector{Float64} = [1.0,], Output::Int64 = 0, ϵ::Float64 = 1e-3, interior_value::Float64 = .5, Enhanced_Cut::Bool = true)
    ## here, x : Lagrangian Multiplier
    ## here, x̂ : the first stage variable
    n = length(x)
    x_interior = [interior_value for i in 1:n]
    ## model to compute F
    optimization_model(x, x̂ = x̂, Output = Output, Enhanced_Cut = Enhanced_Cut)
    optimize!(model)

    # evaluate obj
    f = JuMP.objective_value(model) - x' * (x_interior - x̂) 
    df = - x̂ + JuMP.value(z)


    ## model to compute f*

    if Enhanced_Cut
        @constraint(model, z .== x̂)
        @objective(model, Min, one_vec' * y)
        f_star = JuMP.objective_value(model)

    # evaluate constraints and their gradients
        g1 = (1 - ϵ) * f_star - f
        dg1 = [-df, ]
    else
        g1 = 0.0
        dg1 = [0.0, ]
    end
    Com_G = Dict{Int64, Float64}()
    Com_G[1] = g1

    Com_grad_G = Dict{Int64, Vector{Float64}}()
    Com_grad_G[1] = dg1


    result = Dict()
    result[1] = -f                                    ## -f
    result[2] = [ - df, ]                             ## - df
    result[3] = Com_G
    result[4] = Com_grad_G
    result[5] = maximum(Com_G[i] for i in keys(Com_G))

    return result
end


function optimization_model(x::Vector{Float64};Output::Int64 = 0, x̂::Vector{Float64} = [1.0,], M::Float64 = 1e3, Enhanced_Cut::Bool = true)
    n = length(x)
    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                            "OutputFlag" => Output, 
                                            "Threads" => 1) 
					)

    @variable(model, 0 <= y[1:n] <= 2, Int)
    @variable(model, 0 <= z[1:n] <= 2)

    @constraint(model, 1.5 * y >= z)
    
    one_vec = ones(n)
    if Enhanced_Cut
        @objective(model, Min, one_vec' * y + x' * (z - x̂) )
    else
        @constraint(model, z .== x̂)
        @objective(model, Min, one_vec' * y)
    end


    return model
end


#############################################################################################
####################################    main function   #####################################
#############################################################################################

function LevelSetMethod(x₀::Vector{Float64}; μ::Float64 = 0.5, λ::Float64 = 0.1, Output::Int64 = 0, threshold::Float64 = 1e-5)

    n = length(x₀)

    iter = 1
    α = 1/2

    ## trajectory
    evaluate_function_value = evaluate_function(x₀)
    function_info = FunctionInfo(   Dict(1 => x₀), 
                                    Dict(1 => evaluate_function_value[5]), 
                                    Dict(1 => evaluate_function_value[1]), 
                                    evaluate_function_value[2], 
                                    evaluate_function_value[4], 
                                    evaluate_function_value[3]
                                    )

    ## model for oracle
    model_oracle = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 1)
            )

    @variable(model_oracle, z >= - 50)
    @variable(model_oracle, x[i = 1:n])
    @variable(model_oracle, y <= 0)
    @objective(model_oracle, Min, z)
    oracle_info = ModelInfo(model_oracle, x, y, z)


    # Method 1:  model for next point
    # model_nxt = Model(
    #         optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
    #         "OutputFlag" => Output, 
    #         "Threads" => 1)
    #         )

    # @variable(model_nxt, x1[i = 1:n])
    # @variable(model_nxt, z1 >= -10)
    # @variable(model_nxt, y1)
    # @constraint(model_nxt, level_constraint, α * z1 + (1 - α) * y1 <= 0)
    # nxt_info = ModelInfo(model_nxt, x1, y1, z1)

    while true
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
        if Δ < threshold
            return function_info
        end

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
        
        ######################################################################################################################
        #########################################     next iteration point   #################################################
        ######################################################################################################################

        ######################################################################################################################
        ##  Method 1: obtain the next iteration point Method 1
        # add_constraint(function_info, nxt_info, iter)
        # set_normalized_rhs( level_constraint, level)
        # @objective(model_nxt, Min, (x1 - function_info.x_his[iter])' * (x1 - function_info.x_his[iter]))
        # @info "$model_nxt"
        # optimize!(model_nxt)
        # st = termination_status(model_nxt)
        # if st != MOI.OPTIMAL
        #     @info "Next point model is infeasible"
        #     break
        # end
        # x_nxt = JuMP.value.(x1)

        # but something worry here
        ######################################################################################################################



        ## Method 2:  obtain the next iteration point
        model_nxt = Model(
            optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 1)
            )

        @variable(model_nxt, x1[i = 1:n])
        @variable(model_nxt, z1 >= -1e3)
        @variable(model_nxt, y1)

        @constraint(model_nxt, α * z1 + (1 - α) * y1 <= level)
        @constraint(model_nxt, z1 .>= function_info.f_his[iter] + function_info.df' * (x1 - function_info.x_his[iter]) )
        @constraint(model_nxt, [k in keys(function_info.G)], y1 .>= function_info.G[k] + function_info.dG[k]' * (x1 - function_info.x_his[iter]) )
        @objective(model_nxt, Min, (x1 - function_info.x_his[iter])' * (x1 - function_info.x_his[iter]))
        optimize!(model_nxt)
        st = termination_status(model_nxt)
        if st != MOI.OPTIMAL
            @info "Next point model is infeasible"
            # break
        end
        x_nxt = JuMP.value.(x1)
        ######################################################################################################################
        #####################################################    end   #######################################################
        ######################################################################################################################

        ## save the trajectory
        evaluate_function_value = evaluate_function(x_nxt)

        iter = iter + 1
        function_info.x_his[iter]     = x_nxt
        function_info.G_max_his[iter] = evaluate_function_value[5]
        function_info.f_his[iter]     = evaluate_function_value[1]
        function_info.df              = evaluate_function_value[2]
        function_info.dG              = evaluate_function_value[4]
        function_info.G               = evaluate_function_value[3]

        @info "$Δ, $x_nxt"

    end

end

x₀ = [10.0]
a = LevelSetMethod(x₀, λ = .1, Output = 0)



















