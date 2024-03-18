using JuMP, Test, Statistics, StatsBase, Gurobi, Distributed, Random, Flux, Printf, LinearAlgebra

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
    x_his        :: Dict{Int64, Union{Float64, Vector{Float64}}}  ## record every x_j point
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
"""
    This function is to evaluate the obj, gradient of obj, constraints, gradient of constraint.

        Need to note is that we will return a Dict,
        f     :: Float64 
        df    ::Vector{Float64}
        G     ::Dict{Int64, Float64}  
        dG    ::Dict{Int64, Vector{Float64}}
        G_max :: Float64


"""

struct FunctionCoefficientInfo
    c :: Union{Float64, Vector{Float64}}
    f :: Union{Float64, Vector{Float64}}
    A :: Union{Matrix{Float64}, LinearAlgebra.Adjoint{Float64, Vector{Float64}}}
    b :: Union{Float64, Vector{Float64}}
    B :: Union{Matrix{Float64}, LinearAlgebra.Adjoint{Float64, Vector{Float64}}}
end
c = [1/2, 0, 2, 1];
A = [1, -1.5, 1, -1];
f = 1.0;
b = 0.0;
B = [-1.0];
functionCoefficientInfo = FunctionCoefficientInfo(c, f, A', b, B');


"""
This algorithm is aimed to minimize a convex optimization
   max min f(x) + λ (x -  x̂)
       s.t. G(x) <= 0  ## total m constraints
            x ∈ Q
"""
## function to evaluate the Lagrangian
function evaluate_function(λ::Union{Float64, Vector{Float64}}; 
                            functionCoefficientInfo::FunctionCoefficientInfo = functionCoefficientInfo, 
                            ŷ::Union{Float64, Vector{Float64}} = [3.0, 0.0], Output::Int64 = Output)
    ## here, x : Lagrangian Multiplier
    ## here, x̂ : the first stage variable
    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                            "OutputFlag" => Output, 
                                            "Threads" => 0) 
					)


    @variable(model, x[1:length(functionCoefficientInfo.c)] .≥ 0.0)
    set_integer(x[1])
    set_integer(x[2])
    if length(functionCoefficientInfo.f) == 1
        @variable(model, y)
    else 
        @variable(model, y[1:length(functionCoefficientInfo.f)])
    end

    @constraint(model, functionCoefficientInfo.A * x .== functionCoefficientInfo.b .- functionCoefficientInfo.B * y)
    # @constraint(model, y == - 1.0)

    @objective(model, Min, functionCoefficientInfo.c' * x .- λ' * (y .- ŷ) ) 
    optimize!(model)

    # evaluate obj
    f = JuMP.objective_value(model)
    df = ŷ - JuMP.value.(y)
    # evaluate constraints and their gradients
    d1 = 0.0
    # d2 = [3,9]' * x - 15
    Com_G = Dict{Int64, Float64}()
    Com_G[1] = d1

    d1 = [0.0 for i in 1:length(functionCoefficientInfo.f)]
    # d2 = [3.0, 9.0]
    Com_grad_G = Dict{Int64, Vector{Float64}}()
    Com_grad_G[1] = d1
    # Com_grad_G[2] = d2

    return (    f = -f, 
                df = - [df], 
                G = Com_G, 
                dG = Com_grad_G,
                maxG = maximum(Com_G[i] for i in keys(Com_G)))
end





#############################################################################################
###############################    auxiliary functions   ####################################
#############################################################################################

"""
This function is to constraint the model for solving gap and alpha
"""

function Δ_model_formulation(functionInfo::FunctionInfo, f_star::Float64, iter::Int64; Output::Int64 = Output)
    
    model_alpha = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV),
            "OutputFlag" => Output, 
            "Threads" => 0)
            )

    @variable(model_alpha, z)
    @variable(model_alpha, 0 <= α <= 1)
    @constraint(model_alpha, con[j = 1:iter], z <=  α * ( functionInfo.f_his[j] - f_star) + (1 - α) * functionInfo.G_max_his[j] )
    
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
function add_constraint(functionInfo::FunctionInfo, modelInfo::ModelInfo, iter::Int64)
    m = length(functionInfo.G)

    xⱼ = functionInfo.x_his[iter]
    # add constraints
    @constraint(modelInfo.model, modelInfo.z .≥ functionInfo.f_his[iter] .+ functionInfo.df' * (modelInfo.x .- xⱼ) )
    @constraint(modelInfo.model, [k = 1:m], modelInfo.y .≥ functionInfo.G[k] .+ functionInfo.dG[k]' * (modelInfo.x .- xⱼ) )
end





struct LevelSetMethodParam
    μ             ::Float64                     ## param for adjust α
    λ             ::Union{Float64, Nothing}     ## param for adjust level
    threshold     ::Float64                     ## threshold for Δ
    nxt_bound     ::Float64                     ## lower bound for solving next iteration point π
    max_iter      ::Int64     
    Output        ::Int64                       ## Gurobi Output parameter
    Output_Gap    ::Bool                        ## if True will print Δ info
end
λ_value = nothing; Output = 0; Output_Gap = false; Enhanced_Cut = true; threshold = 1e-2; 
levelSetMethodParam = LevelSetMethodParam(0.9, λ_value, threshold, 1e13, 60, Output, Output_Gap);



#############################################################################################
####################################    main function   #####################################
#############################################################################################
function LevelSetMethod(ŷ::Union{Float64, Vector{Float64}}, x₀::Union{Float64, Vector{Float64}}; 
                        levelSetMethodParam::LevelSetMethodParam = levelSetMethodParam,
                        functionCoefficientInfo::FunctionCoefficientInfo = functionCoefficientInfo)
    (μ, λ, threshold, nxt_bound, max_iter, Output, Output_Gap) = (levelSetMethodParam.μ, levelSetMethodParam.λ, levelSetMethodParam.threshold, levelSetMethodParam.nxt_bound, levelSetMethodParam.max_iter, levelSetMethodParam.Output,levelSetMethodParam.Output_Gap);
    n = length(ŷ)
    iter = 1
    α = 1/2

    ## trajectory
    lagrangianInfo = evaluate_function(x₀; ŷ = ŷ)
    functionInfo = FunctionInfo(   Dict(1 => x₀), 
                                    Dict(1 => lagrangianInfo.maxG), 
                                    Dict(1 => lagrangianInfo.f), 
                                    lagrangianInfo.df, 
                                    lagrangianInfo.dG, 
                                    lagrangianInfo.G
                                    )

    ## model for oracle
    oracleModel = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 0)
            )

    @variable(oracleModel, z ≥ - nxt_bound)
    @variable(oracleModel, x[i = 1:n])
    @variable(oracleModel, y ≤ 0)
    @objective(oracleModel, Min, z)
    oracleInfo = ModelInfo(oracleModel, x, y, z)

    ## model for next point
    nxtModel = Model(
        optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
        "OutputFlag" => Output, 
        "Threads" => 0)
        )

    @variable(nxtModel, x1[i = 1:n])
    @variable(nxtModel, z1 ≥ -nxt_bound)
    @variable(nxtModel, y1)

    while true
        add_constraint(functionInfo, oracleInfo, iter)
        optimize!(oracleModel)

        st = termination_status(oracleModel)
        if st != MOI.OPTIMAL
            @info "oracle is infeasible"
            # break
        end

        f_star = JuMP.objective_value(oracleModel)

        ## formulate alpha model
        gapInfo = Δ_model_formulation(functionInfo, f_star, iter, Output = Output)
        Δ, a_min, a_max = gapInfo[1], gapInfo[2], gapInfo[3]
        if Δ < threshold
            return functionInfo
        end

        if Output_Gap # && (iter % 30 == 0)
            if iter == 1
                println("------------------------------------ Iteration Info --------------------------------------")
                println("Iter |   Gap                              Objective                             Constraint")
            end
            @printf("%3d  |   %5.3g                         %5.3g                              %5.3g\n", iter, Δ,  - functionInfo.f_his[iter], functionInfo.G_max_his[iter])
        end

        ## update α
        if μ/2 <= (α-a_min)/(a_max-a_min) .≤ 1-μ/2
            α = α
        else
            α = (a_min+a_max)/2
        end

        ## update level
        w = α * f_star
        W = minimum( α * functionInfo.f_his[j] + (1-α) * functionInfo.G_max_his[j] for j in 1:iter) 

        λ = iter ≤ 10 ? 0.05 : 0.15
        λ = iter ≥ 20 ? 0.25 : λ
        λ = iter ≥ 30 ? 0.4 : λ
        λ = iter ≥ 40 ? 0.6 : λ
        λ = iter ≥ 50 ? 0.7 : λ
        λ = iter ≥ 55 ? 0.8 : λ

        level = w + λ * (W - w)

        ## obtain the next iteration point
        if iter == 1
            @constraint(nxtModel, levelConstraint, α * z1 + (1 - α) * y1 ≤ level);
        else 
            delete(nxtModel, nxtModel[:levelConstraint]);
            unregister(nxtModel, :levelConstraint);
            @constraint(nxtModel, levelConstraint, α * z1 + (1 - α) * y1 ≤ level);
        end
        @constraint(nxtModel, z1 .≥ functionInfo.f_his[iter] .+ functionInfo.df' * (x1 .- functionInfo.x_his[iter]) )
        @constraint(nxtModel, [k in keys(functionInfo.G)], y1 .≥ functionInfo.G[k] .+ functionInfo.dG[k]' * (x1 .- functionInfo.x_his[iter]) )
        @objective(nxtModel, Min, (x1 .- functionInfo.x_his[iter])' * (x1 .- functionInfo.x_his[iter]))
        optimize!(nxtModel)
        st = termination_status(nxtModel)
        if st != MOI.OPTIMAL
            @info "Next point model is infeasible"
            # break
        end
        x_nxt = JuMP.value.(x1)

        ## =============================================== End ============================================= ##

        ## save the trajectory
        lagrangianInfo = evaluate_function(x_nxt[1]; ŷ = ŷ)

        iter = iter + 1
        functionInfo.x_his[iter]     = x_nxt
        functionInfo.G_max_his[iter] = lagrangianInfo.maxG
        functionInfo.f_his[iter]     = lagrangianInfo.f
        functionInfo.df              = lagrangianInfo.df
        functionInfo.dG              = lagrangianInfo.dG
        functionInfo.G               = lagrangianInfo.G
    end

end

ŷ = -1.0
x₀ = 0.0
a = LevelSetMethod(ŷ, x₀)

a.x_his
a.f_his


λ = a.x_his[length(a.x_his)]
v = - a.f_his[length(a.f_his)] .- π' * ŷ
## Therefore, the cut is of the form θ ≥ λ' * y  + v










