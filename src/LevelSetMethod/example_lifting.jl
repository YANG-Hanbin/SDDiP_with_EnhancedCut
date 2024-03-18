using Pkg
Pkg.activate(".")
using JuMP, Test, Statistics, StatsBase, Gurobi, Distributed, Random, Flux

const GRB_ENV = Gurobi.Env()
include("src/LevelSetMethod/levelmethod_lifting_twobinaryvariable.jl")
n = 1
x̲ = [[0, 1, 2]]
x̄ = [[1, 2, 3]]

a = [[-2, 0, 2]]
b = [[2, 2, -4]]
functionCoefficientInfo = FunctionCoefficientInfo(n ,a, b, x̲, x̄)

## test 1: two binary cases
x̂ = Dict{Any, Any}(:x => [1.2], :λ => [1.0, 0.0])
x₀ = Dict{Any, Any}(:x => [1.2], :λ => [1.2, -.2])
a = LevelSetMethod(x̂, x₀; Output = 0, functionCoefficientInfo = functionCoefficientInfo, maxIter = 20)

a.x_his
a.f_his


π = [a.x_his[length(a.x_his)][:x]; a.x_his[length(a.x_his)][:λ]]
v = - a.f_his[length(a.f_his)] - π' * [x̂[:x]; x̂[:λ]]
# the cut value at the current point is: π' * [x̂[:x]; x̂[:λ]] + v
## Therefore, the cut is of the form θ ≥ π' * x + v


## test 2: single binary cases
include("src/LevelSetMethod/levelmethod_lifting_onebinaryvariable.jl")
x̂ = Dict{Any, Any}(:x => [1.5], :λ => [0.0])
x₀ = Dict{Any, Any}(:x => [1.2], :λ => [1.0])
a = LevelSetMethod(x̂, x₀; Output = 0, functionCoefficientInfo = functionCoefficientInfo, maxIter = 20)

a.x_his
a.f_his


π = [a.x_his[length(a.x_his)][:x]; a.x_his[length(a.x_his)][:λ]]
v = - a.f_his[length(a.f_his)] - π' * [x̂[:x]; x̂[:λ]]
π' * [x̂[:x]; x̂[:λ]] + v
# the cut value at the current point is: π' * [x̂[:x]; x̂[:λ]] + v
