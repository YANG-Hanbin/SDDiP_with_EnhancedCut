#############################################################################################
##########################    auxiliary functions for forward    ############################
#############################################################################################

function add_generator_constraint(StageProblemData::StageData, model_info::ForwardModelInfo;
  A::Matrix{Int64} = [1 2;])    

  @constraint(model_info.model, A * model_info.sum_generator + model_info.x .<= StageProblemData.ū )  ## no more than max num of generators
  @constraint(model_info.model, sum(model_info.y) + model_info.slack .>= model_info.demand )  # satisfy demand
  @constraint(model_info.model, StageProblemData.h * StageProblemData.N 
                          * (A * model_info.sum_generator + model_info.x + StageProblemData.s₀ ) .>= model_info.y )  # no more than capacity
  
end



function add_generator_cut(cut_coefficient::CutCoefficient, model_info::ForwardModelInfo; 
                                                  d::Int64 = 1, n::Int64 = 2, Enhand_Cut::Bool = true,
                                                  A::Matrix{Int64} = [1 2;] )

  l_interior= [.5 for i in 1:n]

  iter = length(keys(cut_coefficient.v))  ## iter num
  k = length(keys(cut_coefficient.v[1]))  ## scenario num

  if Enhand_Cut 
      @constraint(model_info.model, cut[i in 1:iter-1, m in 1:k], model_info.θ >= cut_coefficient.v[i][m] + 
                                            cut_coefficient.π[i][m]' * (model_info.Lt .- l_interior))
  else
      @constraint(model_info.model, cut[i in 1:iter-1, m in 1:k], model_info.θ >= cut_coefficient.v[i][m] + 
                                            cut_coefficient.π[i][m]' * model_info.Lt )
  end
end






#############################################################################################
###################################  function: forward pass #################################
#############################################################################################

"""     forward_step_optimize(StageCoefficient[t], b, x_ancestor, cut_collection, t)
1. Solve the forward problem and return the useful info

cut_coefficient: is the cut info given stage t
StageProblemData: is the param info given stage t
"""


function forward_step_optimize!(StageProblemData::StageData, demand::Vector{Float64}, sum_generator::Vector{Float64}, 
                                cut_coefficient::CutCoefficient; Enhand_Cut::Bool = true,
                                θ_bound::Real = 0.0, d::Int64 = 1, n::Int64 = 2, 
                                A::Matrix{Int64} = [1 2;]  )

                                
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

    add_generator_constraint(StageProblemData, model_info, A = A)
    add_generator_cut(cut_coefficient, model_info, Enhand_Cut = Enhand_Cut, A = A, d = d, n = n)
    @constraint(Q, A * sum_generator + x.== A * Lt)

    @objective(Q, Min, StageProblemData.c1'* x + StageProblemData.c2' * y + StageProblemData.penalty * slack + θ )
    optimize!(Q)

    return [round.(JuMP.value.(Lt)), JuMP.value.(y), JuMP.value(θ), JuMP.objective_value(Q) - JuMP.value(θ)]  ## returen [Lt, y, θ, f]
end

