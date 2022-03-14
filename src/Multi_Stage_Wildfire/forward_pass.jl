#############################################################################################
##########################    auxiliary functions for forward    ############################
#############################################################################################

function add_generator_constraint(indexSets::IndexSets, 
                                  paramDemand::ParamDemand, 
                                  paramOPF::ParamOPF, 
                                  model_info::ForwardModelInfo
                                  )   

    (D, G, L, B, T) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T) 
    (_D, _G, in_L, out_L) = (indexSets._D, indexSets._G, indexSets.in_L, indexSets.out_L) 

    ## constraint 1b 1c
    for l in L
    i = l[1]
    j = l[2]
    @constraint(model_info.model, [t in T], model_info.P[l, t, ω] <= - paramOPF.b[l] * (model_info.θ_angle[i, t, ω] - model_info.θ_angle[j, t, ω] + paramOPF.θmax * (1 - model_info.zl[l, t, ω] ) ) )
    @constraint(model_info.model, [t in T], model_info.P[l, t, ω] >= - paramOPF.b[l] * (model_info.θ_angle[i, t, ω] - model_info.θ_angle[j, t, ω] + paramOPF.θmin * (1 - model_info.zl[l, t, ω] ) ) )
    end

    ## constraint 1d
    @constraint(model_info.model, [l in L, t in 1:T], - paramOPF.W[l, ω] * model_info.zl[l, t, ω] <= model_info.P[l, t, ω] <= paramOPF.W[l, ω] * model_info.zl[l, t, ω] )

    ## constraint 1e
    @constraint(model_info.model, [i in B, t in 1:T], sum(model_info.s[g, t, ω] for g in _G[i]) + sum(model_info.P[(i, j), t, ω] for j in out_L[i] ) .== sum(paramDemand.demand[t][d] * model_info.x[t, d, ω] for d in _D[i]) )

    ## constraint 1f
    @constraint(model_info.model, [g in G, t in 1:T], paramOPF.smin * model_info.zg[g, t, ω] <= model_info.s[g, t, ω] <= paramOPF.smax * model_info.zg[g, t, ω] )

    ## constraint g h i 
    @constraint(model_info.model, [i in B, t in 1:T, d in _D[i]], model_info.zb[i, t, ω] >= model_info.x[t, d, ω] )
    @constraint(model_info.model, [i in B, t in 1:T, d in _G[i]], model_info.zb[i, t, ω] >= model_info.zg[g, t, ω])
    @constraint(model_info.model, [i in B, t in 1:T, j in out_L[i]], model_info.zb[i, t, ω] >= model_info.zl[(i, j), t, ω] )

    ## constraint j k l m
    @constraint(model_info.model, [i in B, t in 1:T, j in in_L[i]], model_info.zb[i, t, ω] >= model_info.zl[(j, i), t, ω] )
    @constraint(model_info.model, [i in B, t in 1:T-1], model_info.zb[i, t, ω] >= model_info.zb[i, t+1, ω] )
    @constraint(model_info.model, [g in G, t in 1:T-1], model_info.zg[g, t, ω] >= model_info.zg[g, t+1, ω] )
    @constraint(model_info.model, [l in L, t in 1:T-1], model_info.zl[l, t, ω] >= model_info.zl[l, t+1, ω] )

end




function add_generator_cut(cut_coefficient::CutCoefficient, model_info::ForwardModelInfo)

    iter = length(keys(cut_coefficient.v))  ## iter num
    k = length(keys(cut_coefficient.v[1]))  ## scenario num

    @constraint(model_info.model, cut[i in 1:iter-1, m in 1:k, ω in Ω], model_info.θ[ω] >= cut_coefficient.v[i][m] + 
                                              cut_coefficient.π[i][m]' * model_info.Lt )

end






#############################################################################################
###################################  function: forward pass #################################
#############################################################################################

"""     forward_step_optimize(StageCoefficient[t], b, x_ancestor, cut_collection, t)
1. Solve the forward problem and return the useful info

cut_coefficient: is the cut info given stage t
StageProblemData: is the param info given stage t
"""


function forward_step_optimize!(indexSets::IndexSets, paramDemand::ParamDemand, paramOPF::ParamOPF, cut_coefficient::CutCoefficient;
                                θ_bound::Real = 0.0)

    (A, n, d) = (binaryInfo.A, binaryInfo.n, binaryInfo.d)                            
    ## construct forward problem (3.1)
    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => 0, 
                                          "Threads" => 1) 
                                          )
    @variable(Q, x[i = 1:d] >= 0, Int)   ## for current state, x is the number of generators will be built in this stage
    @variable(Q, y[i = 1:d] >= 0)        ## amount of electricity
    @variable(Q, Lt[i = 1:n], Bin)       ## stage variable, A * Lt is total number of generators built after this stage
    @variable(Q, slack >=0 )
    @variable(Q, θ >= θ_bound)
    model_info = ForwardModelInfo(Q, x, Lt, y, θ, demand, slack, sum_generator)

    add_generator_constraint(StageProblemData, model_info, binaryInfo = binaryInfo)
    add_generator_cut(cut_coefficient, model_info)
    @constraint(Q, A * sum_generator + x.== A * Lt)

    @objective(Q, Min, StageProblemData.c1'* x + StageProblemData.c2' * y + StageProblemData.penalty * slack + θ )
    optimize!(Q)




    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => 0, 
                                          "Threads" => 1) 
                                          )

    @variable(Q, θ_angle[B, 1:T, Ω]) 
    @variable(Q, P[L, 1:T, Ω]) ## elements in L is Tuple (i, j)
    @variable(Q, s[G, 1:T, Ω] )
    @variable(Q, x[1:T, D, Ω])

    @variable(Q, zg[G, 1:T, Ω] )
    @variable(Q, zb[B, 1:T, Ω] )
    @variable(Q, zl[L, 1:T, Ω] )

    @variable(Q, yb[B, 1:T, Ω])
    @variable(Q, yg[G, 1:T, Ω])
    @variable(Q, yl[L, 1:T, Ω])

    @variable(Q, νb[B, 1:T, Ω ])
    @variable(Q, νg[G, 1:T, Ω ])
    @variable(Q, νl[L, 1:T, Ω ])

    @variable(Q, θ[Ω] >= θ_bound)

    model_info = ForwardModelInfo(Q, θ_angle, P, s, x, zg, zb, zl, yg, yb, yl, νb, νg, νl, θ)

    add_generator_constraint(indexSets, 
                              paramDemand, 
                              paramOPF, 
                              model_info
                                  )   
    
    ## objective function
    @objecive(Q, Min, sum( prob[ω] * ( sum( sum(paramDemand.w[d] * paramDemand.demand[t, d] * (1 - x[t, d, ω]) for d in D ) for t in 1:Ω_rv[ω].τ - 1 ) + θ[ω] ) for ω in Ω)  )

    optimize!(Q)


    return [round.(JuMP.value.(Lt)), JuMP.value.(y), JuMP.value(θ), JuMP.objective_value(Q) - JuMP.value(θ)]  ## returen [Lt, y, θ, f]
end



    # ## constraint 1b 1c
    # for l in L
    #   i = l[1]
    #   j = l[2]
    #   @constraint(Q, [t in T], P[l, t, ω] <= - paramOPF.b[l] * (θ_angle[i, t, ω] - θ_angle[j, t, ω] + ParamOPF.θmax * (1 - zl[l, t, ω] ) ) )
    #   @constraint(Q, [t in T], P[l, t, ω] >= - paramOPF.b[l] * (θ_angle[i, t, ω] - θ_angle[j, t, ω] + ParamOPF.θmin * (1 - zl[l, t, ω] ) ) )
    # end

    # ## constraint 1d
    # @constraint(Q, [l in L, t in 1:T], - ParamOPF.W[l, ω] * zl[l, t, ω] <= P[l, t, ω] <= ParamOPF.W[l, ω] * zl[l, t, ω] )

    # ## constraint 1e
    # @constraint(Q, [i in B, t in 1:T], sum(s[g, t, ω] for g in _G[i]) + sum(P[(i, j), t, ω] for j in out_L[i] ) .== sum(paramDemand.demand[t][d] * x[t, d, ω] for d in _D[i]) )
    
    # ## constraint 1f
    # @constraint(Q, [g in G, t in 1:T], ParamOPF.smin * zg[g, t, ω] <= s[g, t, ω] <= ParamOPF.smax * zg[g, t, ω] )

    # ## constraint g h i 
    # @constraint(Q, [i in B, t in 1:T, d in _D[i]], zb[i, t, ω] >= x[t, d, ω] )
    # @constraint(Q, [i in B, t in 1:T, d in _G[i]], zb[i, t, ω] >= zg[g, t, ω])
    # @constraint(Q, [i in B, t in 1:T, j in out_L[i]], zb[i, t, ω] >= zl[(i, j), t, ω] )

    # ## constraint j k l m
    # @constraint(Q, [i in B, t in 1:T, j in in_L[i]], zb[i, t, ω] >= zl[(j, i), t, ω] )
    # @constraint(Q, [i in B, t in 1:T-1], zb[i, t, ω] >= zb[i, t+1, ω] )
    # @constraint(Q, [g in G, t in 1:T-1], zg[g, t, ω] >= zg[g, t+1, ω] )
    # @constraint(Q, [l in L, t in 1:T-1], zl[l, t, ω] >= zl[l, t+1, ω] )





#     S = [(1, 1, 1), (N, N, N)]
# @variable(model, x1[i=1:N, j=1:N, k=1:N; (i, j, k) in S])
