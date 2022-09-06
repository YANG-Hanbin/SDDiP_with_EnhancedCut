"""
    The data is provided by Jin2011.

"""
using JuMP, Statistics, StatsBase, Gurobi, Distributions, Random, Dates
using CSV, DataFrames, Printf
using JLD2, FileIO


include("src/GenerationExpansion/3.0version/def.jl")
include("src/GenerationExpansion/3.0version/setting.jl")



## the annualized interest rate
r = 0.08
T = 3
## Generator rating
N = [
    1130.0 0 0 0 0 0;
    0 390 0 0 0 0; 
    0 0 380 0 0 0; 
    0 0 0 1180 0 0;
    0 0 0 0 175 0; 
    0 0 0 0 0 560]


ū = [4.0, 10, 10, 1, 45, 4]                             ## maximum number of each type of generators
c = [1.445, 0.795, 0.575, 1.613, 1.650, 1.671] * 10^6   # c_g from table 4, cost/MW to build a generator of type g
mg = [1200, 400, 400, 1200, 500, 600]                   # Install capacity
fuel_price = [3.37, 9.11, 9.11, 9.3e-4, 0, 3.7]
heat_rate = [8844, 7196, 10842, 10400, 1, 8613]
eff = [0.4, 0.56, 0.4, 0.45, 0, 0.48]
om_cost = [4.7, 2.11, 3.66, 0.51, 5.00, 2.98]
s₀ = [1,1,1,1,1,1]                                      # initial number of generators


(probList, 
            stageDataList, 
            Ω, 
            binaryInfo) = dataGeneration(; T = T , initial_demand = 1e8, seed = 1234567)



scenario_sequence = Dict{Int64, Dict{Int64, Any}}()  ## the first index is for scenario index, the second one is for stage
pathList = Vector{Int64}()
push!(pathList, 1)
P = 1.0

recursion_scenario_tree(pathList, P, scenario_sequence, 2, T = T, prob = probList)
scenario_tree = scenario_sequence

save("src/GenerationExpansion/3.0version/testData_3/stageDataList.jld2", "stageDataList", stageDataList)
save("src/GenerationExpansion/3.0version/testData_3/Ω.jld2", "Ω", Ω)
save("src/GenerationExpansion/3.0version/testData_3/binaryInfo.jld2", "binaryInfo", binaryInfo)
save("src/GenerationExpansion/3.0version/testData_3/scenario_sequence.jld2", "scenario_sequence", scenario_sequence)
save("src/GenerationExpansion/3.0version/testData_3/probList.jld2", "probList", probList)
