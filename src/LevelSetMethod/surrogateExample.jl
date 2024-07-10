using Pkg
Pkg.activate(".")
using JuMP, Test, Statistics, StatsBase, Gurobi, Distributed, Random, Flux

const GRB_ENV = Gurobi.Env()
include("src/LevelSetMethod/levelmethod.jl")
include("src/LevelSetMethod/augmented_levelmethod.jl")


## Without Lifting
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


## Litfing 
n = 1
x̲ = [[0, 1, 2]]
x̄ = [[1, 2, 3]]

a = [[-2, 0, 2]]
b = [[2, 2, -4]]
functionCoefficientInfo = FunctionCoefficientInfo(n ,a, b, x̲, x̄)
x̂ = [1.5]
x₀ = [-1.]
surhat = 0.
sur_x = 0.
a = LevelSetMethod(x̂, x₀, λ = .1, Output = 0, functionCoefficientInfo = functionCoefficientInfo, surhat = surhat, sur_x = sur_x)

a.x_his
a.f_his


π = a.x_his[length(a.x_his)]
v = - a.f_his[length(a.f_his)] - π' * [x̂..., surhat]

π' * x + v
## Therefore, the cut is of the form θ ≥ π' * x + v







