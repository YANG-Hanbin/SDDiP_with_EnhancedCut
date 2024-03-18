using Pkg
Pkg.activate(".")
using JuMP, Test, Statistics, StatsBase, Gurobi, Distributed, Random, Flux

const GRB_ENV = Gurobi.Env()
include("src/LevelSetMethod/levelmethod.jl")

## test 1: 1-dimension cases
n = 1
x̲ = [[0, 1, 2]]
x̄ = [[1, 2, 3]]

a = [[-2, 0, 2]]
b = [[2, 2, -4]]
functionCoefficientInfo = FunctionCoefficientInfo(n ,a, b, x̲, x̄)
x̂ = [1.5]
x₀ = [-1.]
a = LevelSetMethod(x̂, x₀, λ = .1, Output = 0, functionCoefficientInfo = functionCoefficientInfo, )

a.x_his
a.f_his


π = a.x_his[length(a.x_his)]
v = - a.f_his[length(a.f_his)] - π' * x̂
## Therefore, the cut is of the form θ ≥ π' * x + v




## test 2: 2-dimension cases
# n = 2
# x̲ = [[0., 3, 7, 12], [0, 4, 8]]
# x̄ = [[3., 7, 12, 15], [4, 8, 15]]

# a = [[-1, 1, -1, 1/3], [1, -1, 5]]
# b = [[9., 6, 19, 7], [3, 8, 10]]

# functionCoefficientInfo = FunctionCoefficientInfo(n ,a, b, x̲, x̄)
# x̂ = [3.0, 0.0]
# x₀ = [10.0 for i in 1:n]
# a = LevelSetMethod(x̂, x₀, λ = .1, Output = 0, functionCoefficientInfo = functionCoefficientInfo, )

# a.x_his
# a.f_his


# π = a.x_his[length(a.x_his)]
# v = - a.f_his[length(a.f_his)] - π' * x̂

# θ_ = π' * x̂ + v
# θ_star = π' * [3.0, 8.0] + v
# f_true = min_fx(functionCoefficientInfo, x̂ = x̂, specify = true)

# gap = f_true[1] - θ_
# ## Therefore, the cut is of the form θ ≥ π' * x + v
