function forward_step_optimize!(StageProblemData::StageData, demand::Vector{Float64}, sum_generator::Vector{Float64}, 
    cut_coefficient::CutCoefficient; Enhand_Cut::Bool = true,
    θ_bound::Real = 0.0, d::Int64 = 1, n::Int64 = 2, 
    A::Matrix{Int64} = [1 2;]  )



    sum_generator= [0.0 for i in 1:n]



    
    ## construct forward problem (3.1)
    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                "OutputFlag" => 1, 
                "Threads" => 1) 
                )
    @variable(Q, x[i = 1:d] >= 0, Int)   ## for current state, x is the number of generators will be built in this stage
    @variable(Q, y[i = 1:d] >= 0)        ## amount of electricity
    @variable(Q, Lt[i = 1:n], Bin)       ## stage variable, A * Lt is total number of generators built after this stage
    @variable(Q, slack >=0 )
    @variable(Q, θ >= 0)
    model_info = ForwardModelInfo(Q, x, Lt, y, θ, Ω[t][j].d, slack, sum_generator)

    add_generator_constraint(StageCoefficient[t], model_info, A = A)
    add_generator_cut(cut_collection[t], model_info, Enhand_Cut = Enhand_Cut, A = A, d = d, n = n)
    @constraint(Q, A * sum_generator + x.== A * Lt)

    @objective(Q, Min, StageCoefficient[t].c1'* x + StageCoefficient[t].c2' * y + StageCoefficient[t].penalty * slack + θ )
    optimize!(Q)

    return [JuMP.value.(Lt), JuMP.value.(y), JuMP.value(θ), JuMP.objective_value(Q) - JuMP.value(θ)]  ## returen [Lt, y, θ, f]
