##########################################################################################
############################  To generate Stage Data  ####################################
##########################################################################################
ū = [3.] # ū = [4.,10.,10.,1.,45.,4.]

binaryDict = binarize_gen(ū)
(A, n, d) = (binaryDict[1], binaryDict[2], binaryDict[3])


c1 = [[30.],[20.]]
c2 = [[14.],[12.5]]

StageCoefficient = Dict{Int64,StageData}()

s₀ = [1]
penalty = 1e2
N = Array{Float64,2}(undef,1,1)
N[1,1] = 2.0
h = 5
for t in 1:T 
    StageCoefficient[t] = StageData(c1[t], c2[t], ū, h, N, s₀, penalty)
end




##########################################################################################
############################  To generate random variable  ###############################
##########################################################################################
T = 2
N_rv = Vector{Int64}()  # the number of realization of each stage
num_Ω = 3
N_rv = [num_Ω for t in 1:T]  ## xxxx 需要update


Random.seed!(12345)

Ω = Dict{Int64,Dict{Int64,RandomVariables}}()   # each stage t, its node are included in Ω[t]
initial_demand = 50  #  5.685e8

for t in 1:T 
    Ω[t] = Dict{Int64,RandomVariables}()
    for i in 1:N_rv[t]
        if t == 1
            Ω[t][i]= RandomVariables([initial_demand])
        else
            Ω[t][i]= RandomVariables( (rand(1)[1]/5+1)*Ω[t-1][i].d )
        end
    end
end




prob = Dict{Int64,Vector{Float64}}()  # P(node in t-1 --> node in t ) = prob[t]
for t in 1:T 
    prob[t] = [0.333 for i in 1:N_rv[t]]
end

# prob[2][1] = .3
# prob[2][2] = .7

















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

################### nonanticipativity for multistage problem #######################

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

                else
                    @constraint(model, [i = 1, j in 2:num_Ω^(T-t)], x[:, t, i] .== x[:, t, j]) 
                end
            end
            recursion_scenario(pathList_copy, P_copy, scenario_sequence, t+1, Ω = Ω, prob = prob, T = T)
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
P = 1.0

recursion_scenario_constraint(pathList, P, scenario_sequence, 2, T = T)
scenario_tree = scenario_sequence
#####################################################################################

@constraint(model, [t in 1:T, ω in 1:W], sum(y[i, t, ω] for i in 1:d) + slack[t, ω] >= Ω[t][scenario_tree[ω][1][t]].d[1] )

@objective(model, Min, sum( sum( scenario_tree[ω][2][t] * (c1[t]' * x[:, t, ω] + c2[t]' * y[:, t, ω] + penalty * slack[t, ω]) for t in 1:T ) for ω in 1:W) )

@info "$model"
optimize!(model)