#############################################################################################
##########################    auxiliary functions for forward    ############################
#############################################################################################

function add_generator_constraint(StageProblemData::StageData, model_info::ForwardModelInfo;
  A::Matrix{Int64} = [1 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                      0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 1 2 4 8 16 32 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 4;]
                  )    
  @constraint(model_info.model, A * model_info.l + model_info.sum_generator .<= StageProblemData.ū )  ## no more than max num of generators
  @constraint(model_info.model, sum(model_info.y) + model_info.slack .>= model_info.demand )  # satisfy demand
  @constraint(model_info.model, StageProblemData.h * StageProblemData.N 
                          * (A * model_info.l + model_info.sum_generator + StageProblemData.s₀ ) .>= model_info.y )  # no more than capacity

end




function add_generator_cut(cut_coefficient::CutCoefficient, model_info::ForwardModelInfo; 
                                                  d::Int64 = 8, n::Int64 = 21, Enhand_Cut::Bool = true,
                                                  A::Matrix{Int64} = [1 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                                                      0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                                                      0 0 0 0 0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0;
                                                                      0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
                                                                      0 0 0 0 0 0 0 0 0 0 0 0 1 2 4 8 16 32 0 0 0;
                                                                      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 4;]
                                                  )
  l_interior= [.5 for i in 1:n]

  iter = length(keys(cut_coefficient.v))  ## iter num
  k = length(keys(cut_coefficient.v[1]))  ## scenario num

  if Enhand_Cut 
      @constraint(model_info.model, cut[i in 1:iter-1, m in 1:k], model_info.θ >= cut_coefficient.v[i][m] + 
                                            cut_coefficient.π[i][m]' * A * (model_info.l .- l_interior))
  else
      @constraint(model_info.model, cut[i in 1:iter-1, m in 1:k], model_info.θ >= cut_coefficient.v[i][m] + 
                                            cut_coefficient.π[i][m]' * A * model_info.l )
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
                                θ_bound::Real = 0.0, d::Int64 = 6, n::Int64 = 21, 
                                A::Matrix{Int64} = [1 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                                    0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                                    0 0 0 0 0 0 0 1 2 4 8 0 0 0 0 0 0 0 0 0 0;
                                                    0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
                                                    0 0 0 0 0 0 0 0 0 0 0 0 1 2 4 8 16 32 0 0 0;
                                                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 4;]
                                )

    ## construct forward problem (3.1)
    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                "OutputFlag" => 0, 
                "Threads" => 1) 
              )
    @variable(Q, l[i = 1:n], Bin)
    @variable(Q, y[i = 1:d] >= 0)
    @variable(Q, slack >=0 )
    @variable(Q, θ >= θ_bound)
    model_info = ForwardModelInfo(Q, l, y, θ, demand, slack, sum_generator)

    add_generator_constraint(StageProblemData, model_info)
    add_generator_cut(cut_coefficient, model_info, Enhand_Cut = Enhand_Cut)
    # @info "$Q"
    @objective(Q, Min, StageProblemData.c1'* A * l + StageProblemData.c2' * y + StageProblemData.penalty * slack + θ )
    optimize!(Q)
    return [A * JuMP.value.(l) + sum_generator, JuMP.value.(y), JuMP.value(θ), JuMP.objective_value(Q) - JuMP.value(θ)]  ## returen [x,y,θ,f]
end
