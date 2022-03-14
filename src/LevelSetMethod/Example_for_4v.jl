struct FunctionCoefficientInfo
    n       :: Int64  ## dim of x
    z_coef  :: Vector{Vector{Float64}}
    y_coef  :: Vector{Vector{Float64}}
    x_lower :: Vector{Vector{Float64}}
    x_upper :: Vector{Vector{Float64}}
end

n = 2
x_lower = [[0., 3, 7, 12], [0, 4, 8]]
x_upper = [[3., 7, 12, 15], [4, 8, 15]]

z_coef = [[-1, 1, -1, 1/3], [1, -1, 5]]
y_coef = [[9., 6, 19, 7], [3, 8, 10]]

functionCoefficientInfo = FunctionCoefficientInfo(n ,z_coef, y_coef, x_lower, x_upper)


function min_fx(functionCoefficientInfo::FunctionCoefficientInfo; M::Float64 = 1e5, x̂::Vector{Float64} = [3.0, 0.0], Output::Int64 = 0, specify::Bool = true)
    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                            "OutputFlag" => Output, 
                                            "Threads" => 1) 
					)

    ## define variables
    @variable(model, z[i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.z_coef[i])] >= 0)
    @variable(model, y[i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.z_coef[i])], Bin)

    @variable(model, x_c[i = 1:functionCoefficientInfo.n] >= 0)  # local copy
    if specify
        @constraint(model, x_c .== x̂)
    end

    ## constraints
    @constraint(model, [i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.z_coef[i])], 
                            z[i, j] - x_c[i]<= 0 )

    @constraint(model, [i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.z_coef[i])], 
                            x_c[i] - (1 - y[i, j]) * M - z[i, j]<= 0)

    @constraint(model, [i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.z_coef[i])], 
                            z[i, j] <= M * y[i, j])

    @constraint(model, [i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.z_coef[i])], 
                            x_c[i] - functionCoefficientInfo.x_upper[i][j] - M * (1 - y[i, j]) <= 0)

    @constraint(model, [i = 1:functionCoefficientInfo.n, j = 1:length(functionCoefficientInfo.z_coef[i])], 
                            functionCoefficientInfo.x_lower[i][j] - M * (1 - y[i, j]) - x_c[i] <= 0 )

    @constraint(model, [i = 1:functionCoefficientInfo.n], sum(y[i, j] for j in 1:length(functionCoefficientInfo.z_coef[i])) == 1)


    ## objective function
    @objective(model, Min, sum(sum(functionCoefficientInfo.z_coef[i][j] * z[i,j] for j in 1:length(functionCoefficientInfo.z_coef[i])) + 
    sum(functionCoefficientInfo.y_coef[i][j] * y[i, j] for j in 1:length(functionCoefficientInfo.z_coef[i])) for i in 1:functionCoefficientInfo.n)) 

    optimize!(model)

    # evaluate obj
    f_star = JuMP.objective_value(model)
    x_star = JuMP.value.(x_c)

    return [f_star, x_star]
end


x̂ = [3.0, 0.0]
min_fx(functionCoefficientInfo, x̂ = x̂, specify = true)





