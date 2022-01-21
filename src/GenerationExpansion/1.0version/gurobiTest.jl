using JuMP, Gurobi, Random

const GRB_ENV = Gurobi.Env()


include("data_struct.jl")

## input data

# include("generationTest.jl")
# include("runtests_small2.jl")
include("runtests_small3.jl")


################################################################################################################################################
###############################################################     Gurobi Model   #############################################################
################################################################################################################################################


W = num_Ω^(T-1) # number of scenarios


model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                            "OutputFlag" => 1, 
                                            "Threads" => 1) 
                                            )
@variable(model, x[i = 1:d, t = 1:T, ω in 1:W] >= 0, Int)   ## for current state, x is the number of generators will be built in this stage
@variable(model, y[i = 1:d, t = 1:T, ω in 1:W] >= 0)        ## amount of electricity
@variable(model, slack[t = 1:T, ω in 1:W] >=0 )
@variable(model, S[i = 1:d, t = 1:T, ω in 1:W] >=0 )

@constraint(model, [t in 1:T, ω in 1:W], S[:, t, ω] .== sum(x[:, j, ω] for j in 1:t ) )
@constraint(model, S .<= ū)

@constraint(model, [t in 1:T, ω in 1:W], y[:,t, ω] .<= h * N * (S[:, t, ω] + s₀))

@constraint(model, [i in 1:W, j in 1:W], x[:, 1, i] .== x[:, 1, j])  ## nonanticipativity for 2-stage problem


########################################################################################################################################################
##################################################### nonanticipativity for multistage problem ########################################################
########################################################################################################################################################

function recursion_scenario_constraint(pathList::Vector{Int64}, P::Float64, scenario_sequence::Dict{Int64, Dict{Int64, Any}}, t::Int64;   
                    Ω::Dict{Int64,Dict{Int64,RandomVariables}} = Ω, prob::Dict{Int64,Vector{Float64}} = prob, T::Int64 = 2)

    if t <= T
        for ω_key in keys(Ω[t])

            pathList_copy = copy(pathList)
            P_copy = copy(P)

            push!(pathList_copy, ω_key)
            P_copy = P_copy * prob[t][ω_key]

            ## nonanticipativity for multi-stage problem
            if t < T
                if haskey(scenario_sequence, 1)
                    first =  maximum(keys(scenario_sequence)) + 1
                    last  =  maximum(keys(scenario_sequence)) + num_Ω^(T-t)
                    @constraint(model, [i = first, j in (first + 1): last], x[:, t, i] .== x[:, t, j]) 
                    @constraint(model, [i = first, j in (first + 1): last], y[:, t, i] .== y[:, t, j]) 
                    @constraint(model, [i = first, j in (first + 1): last], slack[t, i] .== slack[t, j]) 

                else
                    @constraint(model, [i = 1, j in 2:num_Ω^(T-t)], x[:, t, i] .== x[:, t, j]) 
                    @constraint(model, [i = 1, j in 2:num_Ω^(T-t)], y[:, t, i] .== y[:, t, j]) 
                    @constraint(model, [i = 1, j in 2:num_Ω^(T-t)], slack[t, i] .== slack[t, j]) 
                end
            end
            recursion_scenario_constraint(pathList_copy, P_copy, scenario_sequence, t+1, Ω = Ω, prob = prob, T = T)
        end
    else
        if haskey(scenario_sequence, 1)
            scenario_sequence[maximum(keys(scenario_sequence))+1] = Dict(1 => pathList, 2 => P)
        else
            scenario_sequence[1] = Dict(1 => pathList, 2 => P)
        end
        return scenario_sequence
    end

end

scenario_sequence = Dict{Int64, Dict{Int64, Any}}()  ## the first index is for scenario index, the second one is for stage
pathList = Vector{Int64}()
push!(pathList, 1)

recursion_scenario_constraint(pathList, 1.0, scenario_sequence, 2, T = T)
scenario_tree = scenario_sequence
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################



@constraint(model, [t in 1:T, ω in 1:W], sum(y[i, t, ω] for i in 1:d) + slack[t, ω] >= Ω[t][scenario_tree[ω][1][t]].d[1] )

@objective(model, Min, sum( sum( scenario_tree[ω][2] * (c1[t]' * x[:, t, ω] + c2[t]' * y[:, t, ω] + penalty * slack[t, ω]) for t in 1:T ) for ω in 1:W) )

@info "$model"



####################################################### solve the model and display the result ###########################################################
optimize!(model)

JuMP.objective_value(model)
JuMP.value.(x[:, 1, :])
JuMP.value.(x[:, 2, :])
JuMP.value.(x[:, 3, :])

JuMP.value.(y[:, 2, :])

JuMP.value.(slack)