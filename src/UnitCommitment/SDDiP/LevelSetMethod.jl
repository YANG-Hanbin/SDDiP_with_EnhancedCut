"""
    This function is to constraint the model for solving gap and alpha
"""
function Δ_model_formulation(functionHistory::FunctionHistory, f_star::Float64, iter::Int64; Output::Int64 = 0, mipGap::Float64 = 1e-3, timelimit::Int64 = 3)
    
    alphaModel = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                                    "OutputFlag" => Output, 
                                                    "Threads" => 0)); 
    MOI.set(alphaModel, MOI.Silent(), true);
    set_optimizer_attribute(alphaModel, "MIPGap", mipGap);
    set_optimizer_attribute(alphaModel, "TimeLimit", timelimit);

    @variable(alphaModel, z)
    @variable(alphaModel, 0 ≤ α ≤ 1)
    @constraint(alphaModel, con[j = 1:iter], z ≤  α * ( functionHistory.f_his[j] - f_star) + (1 - α) * functionHistory.G_max_his[j] )
    
    # we first compute gap Δ
    @objective(alphaModel, Max, z)
    optimize!(alphaModel)
    st = termination_status(alphaModel)
    Δ = JuMP.objective_value(alphaModel)

    
    ## then we modify above model to compute alpha
    # alpha_min
    @constraint(alphaModel, z .≥ 0)
    @objective(alphaModel, Min, α)
    optimize!(alphaModel)
    a_min = JuMP.value(α)

    # alpha_max
    @objective(alphaModel, Max, α)
    optimize!(alphaModel)
    a_max = JuMP.value(α)

    return Dict(1 => Δ, 2 => a_min, 3 => a_max)

end


"""
    This function is to add constraints for the model f_star and nxt pt.
"""

function add_constraint(currentInfo::CurrentInfo, modelInfo::ModelInfo; κ::Dict = κ)
    m = length(currentInfo.G)

    xⱼ = currentInfo.x
    # add constraints     
    @constraint(modelInfo.model, modelInfo.z .≥ currentInfo.f + sum(sum(currentInfo.df[:λ][g][i] * (modelInfo.xλ[g, i] - xⱼ[:λ][g][i]) for i in 1:κ[g]) for g in keys(currentInfo.df[:λ]))
                                                                    + sum(currentInfo.df[:y][g] * (modelInfo.xy[g] - xⱼ[:y][g]) for g in keys(currentInfo.df[:y])))

    @constraint(modelInfo.model, [k = 1:m], modelInfo.y .≥ currentInfo.G[k] + sum(sum(currentInfo.dG[k][:λ][g][i] * (modelInfo.xλ[g, i] .- xⱼ[:λ][g][i]) for i in 1:κ[g]) for g in keys(currentInfo.df[:λ]))
                                                                                 + sum(currentInfo.dG[k][:y][g] * (modelInfo.xy[g] .- xⱼ[:y][g]) for g in keys(currentInfo.df[:y])) 
                                                                                 )                                                                        
end


"""
    This function is to setup the levelset method.
"""
function setupLevelSetMethod(; stageDecision::Dict = stageDecision, 
                                f_star_value::Float64 = f_star_value,
                                    cutSelection::String = cutSelection, κ::Dict = κ,
                                        Output_Gap::Bool = false, max_iter::Int64 =100, ℓ::Real = .0, λ::Union{Real, Nothing} = .1
                            )
    if cutSelection == "SMC" 
        λ_value = λ; Output = 0; threshold = 1e-4;
        levelSetMethodParam = LevelSetMethodParam(0.9, λ_value, threshold, 1e13, max_iter, Output, Output_Gap, f_star_value);
        x_interior = nothing;
    elseif cutSelection == "ELC"
        λ_value = λ; Output = 0; threshold = 1e-4;
        levelSetMethodParam = LevelSetMethodParam(0.9, λ_value, threshold, 1e13, max_iter, Output, Output_Gap, f_star_value);
        x_interior = Dict(:λ => Dict(g => Dict(i => stageDecision[:λ][g][i] * ℓ .+ (1 - ℓ)/2 for i in 1:κ[g]) for g in keys(stageDecision[:y])), 
                            :y => Dict( g => stageDecision[:y][g] * ℓ .+ (1 - ℓ)/2 for g in keys(stageDecision[:y])));
    elseif cutSelection == "LC"
        λ_value = λ; Output = 0; threshold = 1e-4;
        levelSetMethodParam = LevelSetMethodParam(0.95, λ_value, threshold, 1e15, max_iter, Output, Output_Gap, f_star_value);
        x_interior = nothing;
    end
    x₀ = Dict(:λ => Dict(g => Dict(i => 0.0 for i in 1:κ[g]) for g in keys(stageDecision[:y])), 
                            :y => Dict( g => 0.0 for g in keys(stageDecision[:y])));
    return (x_interior = x_interior, 
            levelSetMethodParam = levelSetMethodParam,
            x₀ = x₀
            )
end

"""
function_info(; x₀::Dict{Symbol, Dict{Int64, Float64}} = x₀,
                model::Model = model,   
                f_star_value::Float64 = f_star_value,
                stageDecision::Dict{Symbol, Dict{Int64, Float64}} = stageDecision,
                cutSelection::String = "LC"
                            )

# Arguments

    1. `x₀::Dict{Symbol, Dict{Int64, Float64}}` : the initial point of the lagrangian dual variables
    2. `model::Model` : the backward model
    3. `f_star_value::Float64` : the optimal value of the current approximate value function
    4. `stageDecision::Dict{Symbol, Dict{Int64, Float64}}` : the decision of the last stage
    5. `cutSelection::String` : the cut selection
  
# Returns
    1. `currentInfo::CurrentInfo` : the current information
  
"""
function function_info(; x₀::Dict = x₀, 
                        model::Model = model, 
                        f_star_value::Float64 = f_star_value, x_interior::Union{Dict, Nothing} = nothing, 
                        stageDecision::Dict = stageDecision, 
                        cutSelection::String = cutSelection, κ::Dict = κ,
                        paramDemand::ParamDemand = paramDemand, paramOPF::ParamOPF = paramOPF, indexSets::IndexSets = indexSets, 
                        ϵ::Float64 = 1e-4, δ::Float64 = 50.)
    if cutSelection == "ELC"
        # objective function
        @objective(model, Min,  sum(model[:h][g] +
                                        paramOPF.C_start[g] * model[:v][g] + 
                                            paramOPF.C_down[g] * model[:w][g] for g in indexSets.G) + 
                                                sum(paramDemand.w[d] * (1 - model[:x][d]) for d in indexSets.D) + sum(model[:θ]) +
                                                    sum(sum(x₀[:λ][g][i] * (stageDecision[:λ][g][i] - model[:λ_copy][g, i]) for i in 1:κ[g]) + x₀[:y][g] * (stageDecision[:y][g] - model[:y_copy][g]) for g in indexSets.G) 
                    );
        ## ==================================================== solve the model and display the result ==================================================== ##
        optimize!(model);
        F  = JuMP.objective_value(model);
        negative_∇F = Dict( :y => Dict(g => JuMP.value(model[:y_copy][g]) - stageDecision[:y][g] for g in indexSets.G),
                            :λ => Dict(g => Dict(i => JuMP.value(model[:λ_copy][g, i]) - stageDecision[:λ][g][i] for i in 1:κ[g]) for g in indexSets.G));

        currentInfo = CurrentInfo( x₀,                                                                                                                                                                                               ## current point
                                    - F - sum(sum(x₀[:λ][g][i] * (x_interior[:λ][g][i] - stageDecision[:λ][g][i]) for i in 1:κ[g]) + x₀[:y][g] * (x_interior[:y][g] - stageDecision[:y][g]) for g in indexSets.G),                   ## obj function value
                                    Dict(1 => f_star_value - F - δ),                                                                                                                                                                 ## constraint value
                                    Dict( :λ => Dict(g => Dict(i => JuMP.value(model[:λ_copy][g, i]) - x_interior[:λ][g][i] for i in 1:κ[g]) for g in indexSets.G),   
                                            :y => Dict(g => JuMP.value(model[:y_copy][g]) - x_interior[:y][g] for g in indexSets.G)),                                                                                                                                                          ## obj gradient
                                    Dict(1 => negative_∇F ));
    elseif cutSelection == "LC"
        # objective function
        @objective(model, Min,  sum(model[:h][g] +
                                        paramOPF.C_start[g] * model[:v][g] + 
                                            paramOPF.C_down[g] * model[:w][g] for g in indexSets.G) + 
                                                sum(paramDemand.w[d] * (1 - model[:x][d]) for d in indexSets.D) + sum(model[:θ]) -
                                                    sum(sum(x₀[:λ][g][i] * model[:λ_copy][g, i] for i in 1:κ[g]) + x₀[:y][g] * model[:y_copy][g] for g in indexSets.G) 
                    );
        ## ==================================================== solve the model and display the result ==================================================== ##
        optimize!(model);
        F  = JuMP.objective_value(model);
        currentInfo = CurrentInfo( x₀,                                                                                                                                                  ## current point
                                    - F - sum(sum(x₀[:λ][g][i] * stageDecision[:λ][g][i] for i in 1:κ[g]) + x₀[:y][g] * stageDecision[:y][g] for g in indexSets.G),                       ## obj function value
                                    Dict(1 => 0.0 ),                                                                                                                                    ## constraint value
                                    Dict( :y => Dict(g => JuMP.value(model[:y_copy][g]) - stageDecision[:y][g] for g in indexSets.G),
                                            :λ => Dict(g => Dict(i => JuMP.value(model[:λ_copy][g, i]) - stageDecision[:λ][g][i] for i in 1:κ[g]) for g in indexSets.G)
                                            ),                                                                                                                                        ## obj gradient
                                    Dict(1 => Dict( :y => Dict(g => 0.0 for g in indexSets.G), :λ => Dict(g => Dict(i => 0.0 for i in 1:κ[g]) for g in indexSets.G)))                     ## constraint gradient
                                        );
    elseif cutSelection == "SMC"
        # objective function
        @objective(model, Min,  sum(model[:h][g] +
                                        paramOPF.C_start[g] * model[:v][g] + 
                                            paramOPF.C_down[g] * model[:w][g] for g in indexSets.G) + 
                                                sum(paramDemand.w[d] * (1 - model[:x][d]) for d in indexSets.D) + sum(model[:θ]) +
                                                    sum(sum(x₀[:λ][g][i] * (stageDecision[:λ][g][i] - model[:λ_copy][g, i]) for i in 1:κ[g]) + x₀[:y][g] * (stageDecision[:y][g] - model[:y_copy][g]) for g in indexSets.G) 
                    );
        ## ==================================================== solve the model and display the result ==================================================== ##
        optimize!(model);
        F  = JuMP.objective_value(model);
        negative_∇F = Dict( :y => Dict(g => JuMP.value(model[:y_copy][g]) - stageDecision[:y][g] for g in indexSets.G),
                            :λ => Dict(g => Dict(i => JuMP.value(model[:λ_copy][g, i]) - stageDecision[:λ][g][i] for i in 1:κ[g]) for g in indexSets.G)
                            );

        currentInfo = CurrentInfo( x₀,                                                                                                ## current point
                                    1/2 * sum(x₀[:y][g] * x₀[:y][g] for g in indexSets.G) + 1/2 * sum(sum(x₀[:λ][g][i] * x₀[:λ][g][i] for i in 1:κ[g]) for g in indexSets.G),                                     ## obj function value
                                    Dict(1 => f_star_value - F - δ),                                                            ## constraint value
                                    x₀,                                                                                               ## obj gradient
                                    Dict(1 => negative_∇F )                                                                           ## constraint gradient
                                    );
    end
    return (currentInfo = currentInfo, currentInfo_f = F)
end    


"""
LevelSetMethod_optimization!(; stageDecision::Dict{Symbol, Dict{Int64, Float64}} = stageDecision,
                                    f_star_value::Float64 = f_star_value,
                                        x_interior::Union{Dict{Symbol, Dict{Int64, Float64}}, Nothing} = nothing,
                                            x₀::Dict{Symbol, Dict{Int64, Float64}} = x₀
                            )

# Arguments

    1. `stageDecision::Dict{Symbol, Dict{Int64, Float64}}` : the decision of the last stage
    2. `f_star_value::Float64` : the optimal value of the current approximate value function
    3. `x_interior::Union{Dict{Symbol, Dict{Int64, Float64}}, Nothing}` : an interior point
    4. `x₀::Dict{Symbol, Dict{Int64, Float64}}` : the initial point of the lagrangian dual variables
    5. `model::Model` : backward model
  
# Returns
    1. `cutInfo::Array{Any,1}` : the cut information
  
"""

function LevelSetMethod_optimization!(; model::Model = model, 
                                        levelSetMethodParam::LevelSetMethodParam = levelSetMethodParam, 
                                        cutSelection::String = cutSelection,## "ELC", "LC", "ShrinkageLC" 
                                        stageDecision::Dict = stageDecision,
                                        x_interior::Union{Dict, Nothing} = nothing, 
                                        x₀::Dict = x₀, κ::Dict = κ, mipGap::Float64 = 1e-3, timelimit::Int64 = 3,
                                        indexSets::IndexSets = indexSets, paramDemand::ParamDemand = paramDemand, paramOPF::ParamOPF = paramOPF, ϵ::Float64 = 1e-4, δ::Float64 = 5.
                                        )

    ## ==================================================== auxiliary function for function information ==================================================== ##
    # μ larger is better
    (μ, λ, threshold, nxt_bound, max_iter, Output, Output_Gap, f_star_value) = (levelSetMethodParam.μ, levelSetMethodParam.λ, levelSetMethodParam.threshold, levelSetMethodParam.nxt_bound, levelSetMethodParam.max_iter, levelSetMethodParam.Output,levelSetMethodParam.Output_Gap, levelSetMethodParam.valValFun);
    (D, G, L, B) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B);
    
    ## ==================================================== Levelset Method ============================================== ##    
    iter = 1;
    α = 1/2;

    # trajectory
    currentInfo, currentInfo_f = function_info(x₀ = x₀, model = model, f_star_value = f_star_value, stageDecision = stageDecision, cutSelection = cutSelection, x_interior = x_interior, κ = κ,
                                                    paramDemand = paramDemand, paramOPF = paramOPF, indexSets = indexSets, ϵ = ϵ, δ = δ);

    functionHistory = FunctionHistory(  Dict(1 => currentInfo.f), 
                                        Dict(1 => maximum(currentInfo.G[k] for k in keys(currentInfo.G)) )
                                        );

    # model for oracle
    oracleModel = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                                    "OutputFlag" => Output, 
                                                    "Threads" => 0)); 
    MOI.set(oracleModel, MOI.Silent(), true);
    set_optimizer_attribute(oracleModel, "MIPGap", mipGap);
    set_optimizer_attribute(oracleModel, "TimeLimit", timelimit);

    ## ==================================================== Levelset Method ============================================== ##
    para_oracle_bound = abs(currentInfo.f);
    z_rhs = 5 * 10^(ceil(log10(para_oracle_bound)));
    @variable(oracleModel, z ≥ - z_rhs);
    @variable(oracleModel, xλ_oracle[g in G, i in 1:κ[g]]);
    @variable(oracleModel, xy_oracle[G]);
    @variable(oracleModel, y ≤ 0);

    @objective(oracleModel, Min, z);
    oracleInfo = ModelInfo(oracleModel, xλ_oracle, xy_oracle, y, z);


    nxtModel = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
    "OutputFlag" => Output, 
    "Threads" => 0)); 
    MOI.set(nxtModel, MOI.Silent(), true);
    set_optimizer_attribute(nxtModel, "MIPGap", mipGap);
    set_optimizer_attribute(nxtModel, "TimeLimit", timelimit);

    @variable(nxtModel, xλ[g in G, i in 1:κ[g]]);
    @variable(nxtModel, xy[G]);
    @variable(nxtModel, z1);
    @variable(nxtModel, y1);
    nxtInfo = ModelInfo(nxtModel, xλ, xy, y1, z1);

    Δ = Inf; τₖ = 1; τₘ = .5; μₖ = 1;

    if cutSelection == "ELC"
        cutInfo =  [ - currentInfo.f - sum(sum(currentInfo.x[:λ][g][i] * x_interior[:λ][g][i] for i in 1:κ[g]) + currentInfo.x[:y][g] * x_interior[:y][g] for g in G), currentInfo.x] 
    elseif cutSelection == "LC"
        # cutInfo =  [ - currentInfo.f - sum(currentInfo.x[:s][g] * stageDecision[:s][g] + currentInfo.x[:y][g] * stageDecision[:y][g] for g in G), currentInfo.x] 
        cutInfo =  [ currentInfo_f, currentInfo.x] 
    elseif cutSelection == "SMC"
        cutInfo =  [ currentInfo_f - sum(sum(currentInfo.x[:λ][g][i] * stageDecision[:λ][g][i] for i in 1:κ[g]) + currentInfo.x[:y][g] * stageDecision[:y][g] for g in G), currentInfo.x] 
    end 

    while true
        add_constraint(currentInfo, oracleInfo; κ = κ);
        optimize!(oracleModel);
        st = termination_status(oracleModel);
        if st == MOI.OPTIMAL
            f_star = JuMP.objective_value(oracleModel);
        else 
            # @info "Oracle Model is $st"
            return (cutInfo = cutInfo, iter = iter)
        end

        # formulate alpha model
        result = Δ_model_formulation(functionHistory, f_star, iter, Output = Output);
        previousΔ = copy.(Δ);
        Δ, a_min, a_max = result[1], result[2], result[3];

        if Output_Gap # && (iter % 30 == 0)
            if iter == 1
                println("------------------------------------ Iteration Info --------------------------------------")
                println("Iter |   Gap                              Objective                             Constraint")
            end
            @printf("%3d  |   %5.3g                         %5.3g                              %5.3g\n", iter, Δ, - currentInfo.f, currentInfo.G[1])
        end

        # push!(gap_list, Δ);
        x₀ = currentInfo.x;
        if (round(previousΔ) > round(Δ)) || ((currentInfo.G[1] ≤ 0.0))
            # x₀ = currentInfo.x;
            τₖ = μₖ * τₖ;
            if cutSelection == "ELC"
                cutInfo =  [ - currentInfo.f - sum(sum(currentInfo.x[:λ][g][i] * x_interior[:λ][g][i] for i in 1:κ[g]) + currentInfo.x[:y][g] * x_interior[:y][g] for g in G), currentInfo.x] 
            elseif cutSelection == "LC"
                # cutInfo =  [ - currentInfo.f - sum(currentInfo.x[:s][g] * stageDecision[:s][g] + currentInfo.x[:y][g] * stageDecision[:y][g] for g in G), currentInfo.x] 
                cutInfo =  [ currentInfo_f, currentInfo.x] 
            elseif cutSelection == "SMC"
                cutInfo =  [ currentInfo_f - sum(sum(currentInfo.x[:λ][g][i] * stageDecision[:λ][g][i] for i in 1:κ[g]) + currentInfo.x[:y][g] * stageDecision[:y][g] for g in G), currentInfo.x] 
            end 
            τₖ = (τₖ + τₘ) / 2;
        end

        # update α
        if μ/2 ≤ (α-a_min)/(a_max-a_min) .≤ 1-μ/2
            α = α;
        else
            α = (a_min+a_max)/2;
        end

        # update level
        w = α * f_star;
        W = minimum( α * functionHistory.f_his[j] + (1-α) * functionHistory.G_max_his[j] for j in 1:iter);

        # λ = iter ≤ 10 ? 0.1 : 0.2;
        # λ = iter ≥ 20 ? 0.3 : λ;
        # λ = iter ≥ 40 ? 0.5 : λ;
        # λ = iter ≥ 60 ? 0.7 : λ;
        
        level = w + λ * (W - w)
        
        ## ==================================================== next iteration point ============================================== ##
        # obtain the next iteration point
        if iter == 1
            @constraint(nxtModel, levelConstraint, α * z1 + (1 - α) * y1 ≤ level);
        else 
            delete(nxtModel, nxtModel[:levelConstraint]);
            unregister(nxtModel, :levelConstraint);
            @constraint(nxtModel, levelConstraint, α * z1 + (1 - α) * y1 ≤ level);
        end
        add_constraint(currentInfo, nxtInfo; κ = κ);
        @objective(nxtModel, Min, sum(sum((xλ[g, i] - x₀[:λ][g][i]) * (xλ[g, i] - x₀[:λ][g][i]) for i in 1:κ[g]) for g in G) +
                                            sum((xy[g] - x₀[:y][g]) * (xy[g] - x₀[:y][g]) for g in G)  #+ 2 * (α * z1 + (1 - α) * y1) * τₖ 
                                    );
        optimize!(nxtModel);
        st = termination_status(nxtModel)
        if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
            x_nxt = Dict(:λ => Dict(g => Dict(i => JuMP.value(xλ[g, i]) for i in 1:κ[g]) for g in indexSets.G), 
                            :y => Dict( g => JuMP.value(xy[g]) for g in indexSets.G));
            λₖ = abs(dual(levelConstraint)); μₖ = λₖ + 1; 
        elseif st == MOI.NUMERICAL_ERROR ## need to figure out why this case happened and fix it
            # @info "Numerical Error occures! -- Build a new nxtModel"
            nxtModel = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 0)); 
            MOI.set(nxtModel, MOI.Silent(), true);
            set_optimizer_attribute(nxtModel, "MIPGap", mipGap);
            set_optimizer_attribute(nxtModel, "TimeLimit", timelimit);

            @variable(nxtModel, xλ[g in G, i in 1:κ[g]]);
            @variable(nxtModel, xy[G]);
            @variable(nxtModel, z1);
            @variable(nxtModel, y1);
            nxtInfo = ModelInfo(nxtModel, xλ, xy, y1, z1);
            @constraint(nxtModel, levelConstraint, α * z1 + (1 - α) * y1 ≤ level);
            add_constraint(currentInfo, nxtInfo; κ = κ);
            @objective(nxtModel, Min, sum(sum((xλ[g, i] - x₀[:λ][g][i]) * (xλ[g, i] - x₀[:λ][g][i]) for i in 1:κ[g]) for g in G) +
                                            sum((xy[g] - x₀[:y][g]) * (xy[g] - x₀[:y][g]) for g in G)  #+ 2 * (α * z1 + (1 - α) * y1) * τₖ 
                                    );
            optimize!(nxtModel);
            st = termination_status(nxtModel);
            if st != MOI.OPTIMAL 
                return (cutInfo = cutInfo, iter = iter)
            end
            x_nxt = Dict(:λ => Dict(g => Dict(i => JuMP.value(xλ[g, i]) for i in 1:κ[g]) for g in indexSets.G), 
                            :y => Dict( g => JuMP.value(xy[g]) for g in indexSets.G));
            λₖ = abs(dual(levelConstraint)); μₖ = λₖ + 1; 
        else
            # @info "Re-compute Next Iteration Point -- change to a safe level!"
            set_normalized_rhs( levelConstraint, w + .99 * (W - w));
            optimize!(nxtModel);
            st = termination_status(nxtModel);
            if st != MOI.OPTIMAL 
                return (cutInfo = cutInfo, iter = iter)
            end
            x_nxt = Dict(:λ => Dict(g => Dict(i => JuMP.value(xλ[g, i]) for i in 1:κ[g]) for g in indexSets.G), 
                            :y => Dict( g => JuMP.value(xy[g]) for g in indexSets.G));
        end

        ## stop rule: gap ≤ .07 * function-value && constraint ≤ 0.05 * LagrangianFunction
        if Δ ≤ threshold * f_star_value || iter > max_iter
            return (cutInfo = cutInfo, iter = iter)
        end
        
        ## ==================================================== end ============================================== ##
        ## save the trajectory
        currentInfo, currentInfo_f = function_info(x₀ = x_nxt, model = model, f_star_value = f_star_value, x_interior = x_interior, stageDecision = stageDecision, cutSelection = cutSelection, paramDemand = paramDemand, paramOPF = paramOPF, indexSets = indexSets, ϵ = ϵ, δ = δ, κ = κ)
        iter = iter + 1;
        functionHistory.f_his[iter] = currentInfo.f;
        functionHistory.G_max_his[iter] = maximum(currentInfo.G[k] for k in keys(currentInfo.G));
    end

end