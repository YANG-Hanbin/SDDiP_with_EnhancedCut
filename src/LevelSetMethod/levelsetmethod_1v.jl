using JuMP, Test, Statistics, StatsBase, Gurobi, Distributed, Random, Flux

const GRB_ENV = Gurobi.Env()

```
This algorithm is aimed to minimize a convex optimization
        min f(x)
       s.t. G(x) <= 0  ## total m constraints
            x ∈ Q
```
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
###############################   function information   ####################################
#############################################################################################
```
This function is to compute the objective function, which should be a convex one.
    return its value :: Float64
```

# function Com_f( x::Vector{Float64} )

#     return x[1]^2 + x[2]^2
# end


# function Com_f( x::Vector{Float64} )

#     return x[1]^2
# end

x = [-18.0]  # -96
x = [-21.0]  # -114
x = [-2.12]  # -0.72
x = [0.0]    # 6  optimal
x = [-4.53]  # -15.19
Com_f(x)


function Com_f( x::Vector{Float64}; M::Float64 = 1e2, x̂::Float64 = 9.0 )
    model = Model(
    optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0)
                  )

    @variable(model, z[i = 1:3] >= 0)
    @variable(model, y[i = 1:3], Bin)
    @variable(model, x̄ >= 0)
    @variable(model, x_c >= 0)  # local copy
    # @constraint(model, x_c == x̂)


    @constraint(model, [i = 1:3], z[i] - x_c<= 0 )
    @constraint(model, [i = 1:3], x_c - (1 - y[i]) * M - z[i]<= 0)
    @constraint(model, [i = 1:3], z[i] <= M * y[i])
    x_lower = [0, 3, 12]
    x_upper = [3, 12, 15]
    @constraint(model, [i = 1:3], x_c - x_upper[i] - M * (1 - y[i]) <= 0)
    @constraint(model, [i = 1:3], x_lower[i] - M * (1 - y[i]) - x_c <= 0 )
    @constraint(model, sum(y) == 1)

    @objective(model, Min, - z[1] + 9 * y[1] + 7 * y[2] + 2 * z[3] - 18 * y[3] + x[1] * (x̂ - x_c) )
    # @objective(model, Min, 0)
    # @objective(model, Min, - z[1] + 9 * y[1] + 7 * y[2] + 2 * z[3] - 18 * y[3] )
    optimize!(model)

    f = JuMP.objective_value(model)

    return -f
end


```
This function is to compute the gradient of the objective function.
    return its gradient :: Vector{Float64}
```

function Com_grad_f( x::Vector{Float64}; M::Float64 = 1e2, x̂::Float64 = 9.0)
    model = Model(
    optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0)
                  )

    @variable(model, z[i = 1:3] >= 0)
    @variable(model, y[i = 1:3], Bin)
    @variable(model, x̄ >= 0)
    @variable(model, x_c >= 0)  # local copy
    # @constraint(model, x_c == x̂)


    @constraint(model, [i = 1:3], z[i] - x_c<= 0 )
    @constraint(model, [i = 1:3], x_c - (1 - y[i]) * M - z[i]<= 0)
    @constraint(model, [i = 1:3], z[i] <= M * y[i])
    x_lower = [0, 3, 12]
    x_upper = [3, 12, 15]
    @constraint(model, [i = 1:3], x_c - x_upper[i] - M * (1 - y[i]) <= 0)
    @constraint(model, [i = 1:3], x_lower[i] - M * (1 - y[i]) - x_c <= 0 )
    @constraint(model, sum(y) == 1)

    @objective(model, Min, - z[1] + 9 * y[1] + 7 * y[2] + 2 * z[3] - 18 * y[3] + x[1] * (x̂ - x_c) )
    # @objective(model, Min, 0)
    # @objective(model, Min, - z[1] + 9 * y[1] + 7 * y[2] + 2 * z[3] - 18 * y[3] )
    optimize!(model)

    return [JuMP.value(x_c) - JuMP.value(x̂), ]
end

# function Com_grad_f( x::Vector{Float64} )
#     return [2*x[1]]
# end




```
This function is to compute values for the constraint functions.
    return their gradients ::Dict{Int64, Float64} index for every components
    G(x) <=> g_k(x) for k = 1:m
```
## basically, this will output a vector
# function Com_G( x::Vector{Float64} )
#     d1 = [5,6]' * x - 10
#     d2 = [3,9]' * x - 15
#     result = Dict{Int64, Float64}()
#     result[1] = d1
#     result[2] = d2
#     return result
# end

# function Com_G( x::Vector{Float64} )
#     d1 = x[1] - 5
#     d2 = 3 - x[1]
#     result = Dict{Int64, Float64}()
#     result[1] = d1
#     result[2] = d2
#     return result
# end
function Com_G( x::Vector{Float64} )
    d1 = 0.0
    # d2 = [3,9]' * x - 15
    result = Dict{Int64, Float64}()
    result[1] = d1
    # result[2] = d2
    return result
end

```
This function is to compute the gradients for every the constraint function.
    return their gradients ::Dict{Int64, Vector{Float64}} index for every components
```
function Com_grad_G( x::Vector{Float64} )
    d1 = [0.0, ]
    # d2 = [3.0, 9.0]
    result = Dict{Int64, Vector{Float64}}()
    result[1] = d1
    # result[2] = d2
    return result
end

# function Com_grad_G( x::Vector{Float64} )
#     d1 = [ 1.0 ]
#     d2 = [-1.0 ]
#     result = Dict{Int64, Vector{Float64}}()
#     result[1] = d1
#     result[2] = d2
#     return result
# end



```
This function is to compute maximal value of all constraints.
    return its value ::Float64 
```
function Com_max_g( x::Vector{Float64} )
    g = Com_G(x)
    return maximum(g[i] for i in keys(g))
end





#############################################################################################
###############################    auxiliary functions   ####################################
#############################################################################################

```
This function is to constraint the model for solving gap and alpha
```

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


```
    This function is to add constraints for the model f_star and nxt pt.
```
function add_constraint(function_info::FunctionInfo, model_info::ModelInfo, iter::Int64)
    m = length(function_info.G)

    xⱼ = function_info.x_his[iter]
    # add constraints
    @constraint(model_info.model, model_info.z .>= function_info.f_his[iter] + function_info.df' * (model_info.x - xⱼ) )
    @constraint(model_info.model, [k = 1:m], model_info.y .>= function_info.G[k] + function_info.dG[k]' * (model_info.x - xⱼ) )
end






x₀ = [1.0,3.0]
x₀ = [3.0]
Output = 0
μ= 0.8
λ= 0.1 ## less is better

threshold = 1e-3
#############################################################################################
####################################    main function   #####################################
#############################################################################################

function LevelSetMethod(x₀::Vector{Float64}; μ::Float64 = 0.5, λ::Float64 = 0.1, Output::Int64 = 0, threshold::Float64 = 1e-5)

    n = length(x₀)

    iter = 1
    α = 1/2

    ## trajectory
    function_info = FunctionInfo(   Dict(1 => x₀), 
                                    Dict(1 => Com_max_g(x₀)), 
                                    Dict(1 => Com_f(x₀)), 
                                    Com_grad_f(x₀), 
                                    Com_grad_G(x₀), 
                                    Com_G(x₀)
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

        result = Δ_model_formulation(function_info, f_star, iter, Output = 1)
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
        iter = iter + 1
        function_info.x_his[iter]     = x_nxt
        function_info.G_max_his[iter] = Com_max_g(x_nxt)
        function_info.f_his[iter]     = Com_f(x_nxt)
        function_info.df              = Com_grad_f(x_nxt)
        function_info.dG              = Com_grad_G(x_nxt)
        function_info.G               = Com_G(x_nxt)

        @info "$Δ, $x_nxt"

    end

end

x₀ = [40.0]
a = LevelSetMethod(x₀, λ = .1, Output = 0)



## test case
x₀ = [1.0,3.0]
x₀ = [9.0]
Output = 1
μ= 0.8
λ= 0.1 ## less is better

threshold = 1e-3
a = LevelSetMethod(x₀, λ = .1, Output = 0)

a.x_his

















function function_information( x::Vector{Float64})
    f(x) = x[1]^2 + x[2]^2
 
    g1(x) = [5,6]' * x - 10
    g2(x) = [3,9]' * x - 15
    G_value = [g1(x),g2(x)]
    

    f_value = f(x)
    df_value = df(x)
    G_max_value = maximum(G_value[i] for i in keys(G_value))
    G_value_dict = Dict(i => G_value[i] for i in 1:length(G_value))
    dG_value_dict = Dict(1 => gradient(g1, x)[1], 2 => gradient(g2, x)[1])

end









