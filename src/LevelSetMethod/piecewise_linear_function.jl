"""
    This script is used to test the Lagrangian cut for piece-wise linear function using the level method:


        min f(x) s.t. x̲ <= x <= x̄     
    where f(x) is a piecewise linear function with m pieces, and x is a vector of dimension n,
    f(x) = min ∑ᵢ yᵢ(aᵢx + bᵢ) 
          s.t. ∑ᵢ yᵢ = 1
               yᵢ ∈ {0, 1}
               x̲_i - M (1 - yᵢ)<= x <= x̄_i - M (1 - yᵢ)
"""
struct FunctionCoefficientInfo
    n       :: Int64                        ## dim of x
    a       :: Vector{Vector{Float64}}
    b       :: Vector{Vector{Float64}}
    x̲       :: Vector{Vector{Float64}}
    x̄       :: Vector{Vector{Float64}}
end


function min_fx(functionCoefficientInfo::FunctionCoefficientInfo; 
                        M::Float64 = 1e5, x̂::Vector{Float64} = [3.0, 0.0], Output::Int64 = 0, nonancipativity::Bool = true
                )

    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                            "OutputFlag" => Output, 
                                            "Threads" => 1) 
					)

    ## define variables
    @variable(model, z[i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.a[i])])                # zᵢ = xyᵢ
    @variable(model, y[i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.a[i])], Bin)           # yᵢ ∈ {0, 1}
    @variable(model, x_c[i = 1:functionCoefficientInfo.n])                                                          # local copy

    # nonancipativity
    if true 
        @constraint(model, x_c .== x̂)
    end

    # McCormick relaxation
    @constraint(model, [i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.a[i])], 
                            z[i, j] - x_c[i] ≤ 0 )

    @constraint(model, [i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.a[i])], 
                            x_c[i] - (1 - y[i, j]) * M - z[i, j] ≤ 0)

    @constraint(model, [i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.a[i])], 
                            z[i, j] ≤ M * y[i, j])
    @constraint(model, [i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.a[i])], 
                            z[i, j] ≥ - M * y[i, j])
                            
    # choose x region
    @constraint(model, [i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.a[i])], 
                            x_c[i] - functionCoefficientInfo.x̄[i][j] - M * (1 - y[i, j]) ≤ 0)

    @constraint(model, [i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.a[i])], 
                            functionCoefficientInfo.x̲[i][j] - M * (1 - y[i, j]) - x_c[i] ≤ 0 )

    @constraint(model, [i = 1:functionCoefficientInfo.n], sum(y[i, j] for j in 1:length(functionCoefficientInfo.a[i])) == 1)


    ## objective function
    @objective(model, Min, sum(sum(functionCoefficientInfo.a[i][j] * z[i,j] for j in 1:length(functionCoefficientInfo.a[i])) + 
    sum(functionCoefficientInfo.b[i][j] * y[i, j] for j in 1:length(functionCoefficientInfo.a[i])) for i in 1:functionCoefficientInfo.n)) 

    optimize!(model)

    # evaluate obj
    f_star = JuMP.objective_value(model)
    x_star = JuMP.value.(x_c)

    return [f_star, x_star]
end








