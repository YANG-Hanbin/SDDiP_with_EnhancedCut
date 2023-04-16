using JuMP, Gurobi, Random

const GRB_ENV = Gurobi.Env()
include("def.jl")

## input data
# include("generationTest.jl")


################################################################################################################################################
############################################################     Gurobi function   #############################################################
################################################################################################################################################
function gurobiOptimize!(Ω::Dict{Int64,ScenarioData}, 
                        indexSets::IndexSets = indexSets,
                        prob::Prob, 
                        stageData::StageData; 
                        mipGap::Float64 = 1e-3, timeLimit::Float64 = 1e3, outputFlag::Int64 = 1)

    ## construct forward second-stage problem 
    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => outputFlag, 
                                          "Threads" => 0,
                                          "MIPGap" => mipGap, 
                                          "TimeLimit" => timeLimit) 
                                          )

    @variable(model, π[indexSets.N, indexSets.Ω] ≥ 0)  
    @variable(model, x[indexSets.D], Bin) 

    @constraint(model, sum(stageData.c[(i,j)] * x[(i,j)] for (i,j) in indexSet.D) ≤ stageData.b ) 
    @constraint(model, [ω in indexSets.Ω], π[Ω[ω].t, ω] == 1)
    @constraint(model, [(i,j) in indexSets.D, ω in indexSets.Ω], π[i, ω] - prob.q[(i,j)] * π[j, ω] ≥ 0)
    @constraint(model, [(i,j) in indexSets.Dᶜ,ω in indexSets.Ω], π[i, ω] - prob.r[(i,j)] * π[j, ω] ≥ 0)
    @constraint(model, [(i,j) in indexSets.D, ω in indexSets.Ω], π[i, ω] - prob.r[(i,j)] * π[j, ω] ≥ 
                                                    (prob.q[(i,j)] - prob.r[(i,j)]) * Ω[ω].ϕ[j] * x[(i,j)]
                    )

    @objective(model, Min, sum(π[Ω[ω].s, ω] * prob.p[ω] for ω in indexSets.Ω) )


    optimize!(model)
    ####################################################### solve the model and display the result ###########################################################
    gurobiResult = (OPT = JuMP.objective_value(model), 
                        statevariable_x = JuMP.value.(x)
                        )

    return gurobiResult
end
