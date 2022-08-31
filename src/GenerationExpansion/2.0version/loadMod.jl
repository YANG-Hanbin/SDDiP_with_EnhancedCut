using JuMP, Test, Statistics, StatsBase, Gurobi, Distributed, ParallelDataTransfer, Random, DataFrames, Dates, Printf

const GRB_ENV = Gurobi.Env()


include("src/GenerationExpansion/2.0version/def.jl")
include("src/GenerationExpansion/2.0version/backwardPass.jl")
include("src/GenerationExpansion/2.0version/forwardPass.jl")
include("src/GenerationExpansion/2.0version/extFormGurobi.jl")
include("src/GenerationExpansion/2.0version/setting.jl")
include("src/GenerationExpansion/2.0version/generationTest.jl")
include("src/GenerationExpansion/2.0version/runtests_small2.jl")  ## M = 4

#############################################################################################
####################################    main function   #####################################
#############################################################################################
max_iter = 200; ϵ = 1e-2; Enhanced_Cut = true