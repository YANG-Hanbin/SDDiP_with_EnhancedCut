## this file is to generate data

#  the matrix A is to help to form decision variables X
# A = [   1 2 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0 1 2 4 8 16 32 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 4;  ]


ū = [3] # ū = [4,10,10,1,45,4]

row_num = size(ū); row_num = row_num[1];

var_num = floor.(Int, log.(2,ū) ) .+ 1; col_num = sum(Int, var_num);

A = zeros(row_num, col_num)
for i in 1:row_num
    l = sum(var_num[l] for l in 1:i)
    for j in (l+1-var_num[i]):(l+var_num[i]-var_num[i])
        A[i,j] = 2^(j-(l+1-var_num[i]))
        # println(1)
    end
end

A

(d,n) = size(A)

struct RandomVariables
    d::Vector{Float64}
end


struct StageData ## with the assumption that only b has stochasticity
    c1       ::Vector{Float64}
    c2       ::Vector{Float64}
    ū        ::Vector{Float64}
    h        ::Float64
    N        ::Matrix{Float64}
    s₀       ::Vector{Float64}
    penalty  ::Float64
end



##########################################################################################
############################  To generate Stage Data  ####################################
##########################################################################################
T = 2



ū = [3]

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
M =4
N = Vector{Int64}()  # the number of realization of each stage
N_rv = [2, 2]  ## xxxx 需要update


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
    prob[t] = [0.5 for i in 1:N_rv[t]]
end

# prob[2][1] = .3
# prob[2][2] = .7

















################################################################################################################################################
###############################################################     Gurobi Model   #############################################################
################################################################################################################################################


W = 2^(T-1) # number of scenarios


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
@constraint(model, [t in 1:T, ω in 1:W], sum(y[i, t, ω] for i in 1:d) + slack[t, ω] >= Ω[t][ω].d[1] )
@constraint(model, [t in 1:T, ω in 1:W], y[:,t, ω] .<= h * N * (S[:, t, ω] + s₀))

@constraint(model, [i in 1:W, j in 1:W], x[:, 1, i] .== x[:, 1, j])  ## nonanticipativity for 2-stage problem

@objective(model, Min, sum( sum( prob[t][ω] * (c1[t]' * x[:, t, ω] + c2[t]' * y[:, t, ω] + penalty * slack[t, ω]) for t in 1:T ) for ω in 1:W) )
@info "$model"
optimize!(model)

































