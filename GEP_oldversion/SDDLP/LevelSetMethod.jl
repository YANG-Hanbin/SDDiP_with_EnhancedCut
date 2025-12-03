#############################################################################################
###########################  auxiliary functions for level set method #######################
#############################################################################################

"""
This function is to constraint the model for solving gap and alpha
"""
function Δ_model_formulation(
    functionHistory::FunctionHistory, 
    f_star::Float64, 
    iter::Int64; 
    Output::Int64 = 0
)::Dict{Int64, Float64}
    
    alphaModel = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 0
        )
    ); 
    MOI.set(alphaModel, MOI.Silent(), true);
    set_optimizer_attribute(alphaModel, "MIPGap", 1e-4);
    set_optimizer_attribute(alphaModel, "TimeLimit", 5);

    @variable(alphaModel, z)
    @variable(alphaModel, 0 ≤ α ≤ 1)
    @constraint(
        alphaModel, 
        con[j = 1:iter], 
        z ≤  α * ( functionHistory.f_his[j] - f_star) + (1 - α) * functionHistory.G_max_his[j] 
    );
    
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

    return Dict(
        1 => Δ, 
        2 => a_min, 
        3 => a_max
    )
end


"""
    This function is to add constraints for the model f_star and nxt pt.
"""
function add_constraint(
    currentInfo::CurrentInfo, 
    modelInfo::ModelInfo
)::Nothing
    m = length(currentInfo.G);
    xⱼ = currentInfo.x
    # add constraints     
    @constraint(
        modelInfo.model, 
        modelInfo.z .≥ currentInfo.f + 
        currentInfo.df[:St]' * (modelInfo.x .- xⱼ[:St]) + 
        sum(
            sum(
                currentInfo.df[:sur][g][k] * (modelInfo.x_sur[g, k] .- xⱼ[:sur][g][k]) for k in keys(currentInfo.df[:sur][g])
                ) for g in 1:length(modelInfo.x)
        )                                                              
    );

    @constraint(
        modelInfo.model, 
        [k = 1:m], 
        modelInfo.y .≥ currentInfo.G[k] + 
        sum(currentInfo.dG[k][:St] .* (modelInfo.x .- xⱼ[:St])) + 
        sum(
            sum(
                currentInfo.dG[k][:sur][g][km] * (modelInfo.x_sur[g, km] .- xⱼ[:sur][g][km]) for km in keys(currentInfo.dG[k][:sur][g])
            ) for g in 1:length(modelInfo.x)
        )                                                                          
    ); 
    return 
end


"""
    This function is to collect the information from the objecive f, and constraints G
"""
function FuncInfo_LevelSetMethod(
    x₀::Dict{Symbol, Any}; 
    backwardInfo::BackwardModelInfo = backwardInfo,
    cutSelection::String = cutSelection, 
    f_star_value::Union{Float64, Nothing} = f_star_value, 
    stageData::StageData = stageData, 
    Ŝ::Dict{Symbol, Any} =  Ŝ, 
    S̃::Union{Dict{Symbol, Any}, Nothing} =  S̃, 
    binaryInfo::BinaryInfo = binaryInfo,
    param::NamedTuple = param
)::CurrentInfo

    if cutSelection == "ELC"
        @objective(
            backwardInfo.model, 
            Min, 
            stageData.c1' * backwardInfo.x + 
            stageData.c2' * backwardInfo.y + 
            backwardInfo.θ + 
            stageData.penalty * backwardInfo.slack + 
            x₀[:St]' * ( Ŝ[:St] .- backwardInfo.Sc) + 
            sum(
                sum(
                    x₀[:sur][g][k] * (Ŝ[:sur][g][k] - backwardInfo.model[:sur_copy][g, k]) for k in keys(Ŝ[:sur][g])
                ) for g in 1:binaryInfo.d
            ) 
        );

        optimize!(backwardInfo.model);
        if param.tightness
            F_solution = ( 
                F = JuMP.objective_value(backwardInfo.model), 
                negative_∇F = Dict(
                    :St => round.(JuMP.value.(backwardInfo.Sc), digits = 2) - Ŝ[:St], 
                    :sur => Dict(
                        g => Dict(k => round(JuMP.value(backwardInfo.model[:sur_copy][g, k]), digits = 2) - Ŝ[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d)
                    )
            );    
        else
            F_solution = ( 
                F = JuMP.objective_value(backwardInfo.model), 
                negative_∇F = Dict(
                    :St => JuMP.value.(backwardInfo.Sc) - Ŝ[:St], 
                    :sur => Dict(
                        g => Dict(k => JuMP.value(backwardInfo.model[:sur_copy][g, k]) - Ŝ[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d)
                    )
            );   
        end
                        
        currentInfo  = CurrentInfo(
            x₀, 
            - F_solution.F - x₀[:St]' * (S̃[:St] .-  Ŝ[:St]) - sum(sum(x₀[:sur][g][k] * (S̃[:sur][g][k] - Ŝ[:sur][g][k]) for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d),
            Dict(1 => f_star_value - F_solution.F - param.ε),
            Dict( 
                :St => JuMP.value.(backwardInfo.Sc) .- S̃[:St], 
                :sur => Dict(g => Dict(k => JuMP.value(backwardInfo.model[:sur_copy][g, k]) - S̃[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d)
            ),
            Dict(1 => F_solution.negative_∇F),
            F_solution.F 
        );

    elseif cutSelection ∈ ["LC", "SBC"]
        @objective(
            backwardInfo.model, 
            Min, 
            stageData.c1' * backwardInfo.x + 
            stageData.c2' * backwardInfo.y + 
            backwardInfo.θ + 
            stageData.penalty * backwardInfo.slack - 
            x₀[:St]' * backwardInfo.Sc - sum(sum(x₀[:sur][g][k] * backwardInfo.model[:sur_copy][g, k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d) 
        );
        optimize!(backwardInfo.model);
        if param.tightness
            F_solution = ( 
                F = JuMP.objective_value(backwardInfo.model), 
                negative_∇F = Dict( 
                    :St => round.(JuMP.value.(backwardInfo.Sc), digits = 2), 
                    :sur => Dict(g => Dict(k => round(JuMP.value(backwardInfo.model[:sur_copy][g, k]), digits = 2) for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d)
                ),
                zero_∇F = Dict( 
                    :St => round.(JuMP.value.(backwardInfo.Sc), digits = 2) .* 0, 
                    :sur => Dict(g => Dict(k => 0.0 for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d)
                )
            );           
        else
            F_solution = ( 
                F = JuMP.objective_value(backwardInfo.model), 
                negative_∇F = Dict( 
                    :St => JuMP.value.(backwardInfo.Sc), 
                    :sur => Dict(g => Dict(k => JuMP.value(backwardInfo.model[:sur_copy][g, k]) for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d)
                ),
                zero_∇F = Dict( 
                    :St => JuMP.value.(backwardInfo.Sc) .* 0, 
                    :sur => Dict(g => Dict(k => 0.0 for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d)
                )
            );           
        end        
                        
        currentInfo  = CurrentInfo(
            x₀, 
            - F_solution.F - x₀[:St]' * Ŝ[:St] - sum(sum(x₀[:sur][g][k] * Ŝ[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d),                  
            Dict(1 => 0.0),
            Dict( 
                :St => JuMP.value.(backwardInfo.Sc) .- Ŝ[:St],                  
                :sur => Dict(g => Dict(k => JuMP.value(backwardInfo.model[:sur_copy][g, k]) - Ŝ[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d)),
            Dict(1 => F_solution.zero_∇F), 
            F_solution.F
        );
    elseif cutSelection == "ShrinkageLC"
        @objective(
            backwardInfo.model, 
            Min, 
            stageData.c1' * backwardInfo.x + 
            stageData.c2' * backwardInfo.y + 
            backwardInfo.θ + 
            stageData.penalty * backwardInfo.slack + 
            x₀[:St]' * ( Ŝ[:St] .- backwardInfo.Sc) + 
            sum(sum(x₀[:sur][g][k] * (Ŝ[:sur][g][k] - backwardInfo.model[:sur_copy][g, k]) for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d) 
        );

        optimize!(backwardInfo.model);
        if param.tightness
            F_solution = ( 
                F = JuMP.objective_value(backwardInfo.model), 
                negative_∇F = Dict( 
                    :St => round.(JuMP.value.(backwardInfo.Sc), digits = 2) - Ŝ[:St], 
                    :sur => Dict(g => Dict(k => round(JuMP.value(backwardInfo.model[:sur_copy][g, k]), digits = 2) - Ŝ[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d)
                )
            );           
        else
            F_solution = ( 
                F = JuMP.objective_value(backwardInfo.model), 
                negative_∇F = Dict( 
                    :St => JuMP.value.(backwardInfo.Sc) - Ŝ[:St], 
                    :sur => Dict(g => Dict(k => JuMP.value(backwardInfo.model[:sur_copy][g, k]) - Ŝ[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d))
            );   
        end        
                        
        currentInfo  = CurrentInfo(
            x₀, 
            1/2 * sum(x₀[:St] .* x₀[:St]) + 1/2 * sum(sum(x₀[:sur][g][k] * x₀[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d),                
            Dict(1 => f_star_value - F_solution.F - param.ε),
            Dict( 
                :St => x₀[:St], 
                :sur => Dict(g => Dict(k => x₀[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d)),    
            Dict(1 => F_solution.negative_∇F), 
            F_solution.F
        );

    end
    return currentInfo
end

"""
    This function is to collect the necessary parameters for the level set method.
"""
function setupLevelsetPara(
    forwardInfo::ForwardModelInfo, 
    backwardInfo::BackwardModelInfo,
    stageData::StageData, 
    demand::Vector{Float64}, 
    stateInfo::Any,
    param::NamedTuple;                                            
    binaryInfo::BinaryInfo = binaryInfo,   
    Output_Gap::Bool = false, 
    λ::Union{Float64, Nothing} = .3
)::NamedTuple
    L̂ = stateInfo.stageSolution; sur = stateInfo.stageSur; 
    Ŝ = Dict( :St => L̂, :sur => Dict(g => Dict(k => sur[g][k] for k in keys(sur[g])) for g in 1:binaryInfo.d));
    if param.cutSelection == "ELC"
        forward_modify_constraints!(
            forwardInfo,                   
            stageData,                
            demand,                  
            L̂,                     
            binaryInfo = binaryInfo                         
        );

        optimize!(forwardInfo.model); 
        f_star_value = JuMP.objective_value(forwardInfo.model);
        L̃ = Dict( 
            :St => L̂ .* 0 .+ .5,    
            :sur => Dict(
                g => Dict(k => .5 for k in keys(sur[g])) for g in 1:binaryInfo.d
            )
        );

        Output = 0; threshold = 1.0; 
        
        levelSetMethodParam = LevelSetMethodParam(
            0.95, 
            λ, 
            threshold,                            
            param.nxt_bound, 
            param.MaxIter, 
            Output, Output_Gap,                               
            Ŝ, 
            param.cutSelection, 
            L̃, 
            f_star_value
        );
    elseif param.cutSelection == "ShrinkageLC" 
        forward_modify_constraints!(
            forwardInfo,                 
            stageData,                
            demand,                  
            L̂,                     
            binaryInfo = binaryInfo                         
        );
        optimize!(forwardInfo.model); f_star_value = JuMP.objective_value(forwardInfo.model);
        Output = 0; threshold = 1.0; 
        levelSetMethodParam = LevelSetMethodParam(
            0.95, 
            λ, 
            threshold, 
            param.nxt_bound, 
            param.MaxIter, 
            Output, 
            Output_Gap,
            Ŝ,  
            param.cutSelection, 
            nothing, 
            f_star_value
        );
    elseif param.cutSelection == "LC" 
        f_star_value = 0.0;
        Output = 0; threshold = 1.0; 
        levelSetMethodParam = LevelSetMethodParam(
            0.95, 
            λ, 
            threshold, 
            param.nxt_bound, 
            param.MaxIter, 
            Output, 
            Output_Gap, 
            Ŝ,  
            param.cutSelection, 
            nothing, 
            f_star_value
        );
    elseif param.cutSelection == "SBC" 
        f_star_value = 0.0;
        Output = 0; threshold = 1.0; 
        levelSetMethodParam = LevelSetMethodParam(
            0.95, 
            λ, 
            threshold, 
            param.nxt_bound, 
            param.MaxIter, 
            Output, 
            Output_Gap, 
            Ŝ,  
            param.cutSelection, 
            nothing, 
            f_star_value
        );
    end

    if param.cutSelection == "SBC"
        x₀ = get_Benders_coefficient!(
            backwardInfo, 
            stateInfo
        );
    else
        x₀ = Dict( 
            :St => L̂ .* 0.0,                 
            :sur => Dict(
                g => Dict(
                    k => 0.0 for k in keys(sur[g])
                ) for g in 1:binaryInfo.d
            )
        );
    end

    return (
        levelSetMethodParam = levelSetMethodParam, 
        x₀ = x₀
    )
end

######################################################################################################################
## -------------------------------------------------- Main Function -------------------------------------------- ##
######################################################################################################################
function LevelSetMethod_optimization!( 
    backwardInfo::BackwardModelInfo, 
    x₀::Dict{Symbol, Any}; 
    levelSetMethodParam::LevelSetMethodParam = levelSetMethodParam,   
    stageData::StageData = stageData,                         
    binaryInfo::BinaryInfo = binaryInfo,
    param::NamedTuple = param
)::Any 
    
    ######################################################################################################################
    ###############################   auxiliary function for function information   ######################################
    ######################################################################################################################
    ##  μ larger is better
    (μ, λ, threshold, nxt_bound, max_iter, Output, Output_Gap) = (levelSetMethodParam.μ, levelSetMethodParam.λ, levelSetMethodParam.threshold, levelSetMethodParam.nxt_bound, levelSetMethodParam.max_iter, levelSetMethodParam.Output,levelSetMethodParam.Output_Gap);
    (Ŝ, cutSelection, S̃, f_star_value) = (levelSetMethodParam.Ŝ, levelSetMethodParam.cutSelection, levelSetMethodParam.S̃, levelSetMethodParam.f_star_value)
    (A, n, d) = (binaryInfo.A, binaryInfo.n, binaryInfo.d);
    
    ## ==================================================== Levelset Method ============================================== ##
    iter = 1;
    α = 1/2;

    ## trajectory
    currentInfo = FuncInfo_LevelSetMethod(
        x₀, 
        cutSelection = cutSelection, 
        backwardInfo = backwardInfo, 
        f_star_value = f_star_value, 
        stageData = stageData, 
        Ŝ = Ŝ, 
        S̃ = S̃, 
        binaryInfo = binaryInfo,
        param = param
    );


    if cutSelection == "ELC" 
        cutInfo =  [ 
            - currentInfo.f - 
            currentInfo.x[:St]' *  S̃[:St] - 
            sum(sum(currentInfo.x[:sur][g][k] * S̃[:sur][g][k] for k in keys(S̃[:sur][g])) for g in 1:binaryInfo.d),                                                              
            currentInfo.x
        ];
    elseif cutSelection == "LC"
        cutInfo = [ 
            - currentInfo.f - 
            currentInfo.x[:St]' *  Ŝ[:St] -  
            sum(sum(currentInfo.x[:sur][g][k] * Ŝ[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d),                                                              
            currentInfo.x
        ];
    elseif cutSelection == "ShrinkageLC"
        cutInfo = [ 
            currentInfo.S_at_x̂ - 
            currentInfo.x[:St]' *  Ŝ[:St] -  
            sum(sum(currentInfo.x[:sur][g][k] * Ŝ[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d),  
            currentInfo.x
        ];
    elseif cutSelection == "SBC"
        cutInfo = [ 
            - currentInfo.f - 
            currentInfo.x[:St]' *  Ŝ[:St] -  
            sum(sum(currentInfo.x[:sur][g][k] * Ŝ[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d),                                                              
            currentInfo.x
        ];

        return cutInfo
    end

    functionHistory = FunctionHistory(  
        Dict(1 => currentInfo.f), 
        Dict(1 => maximum(currentInfo.G[k] for k in keys(currentInfo.G)) )                              
    );

    ## model for oracle
    oracleModel = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                        "OutputFlag" => Output, 
                                        "Threads" => 0)); 
    MOI.set(oracleModel, MOI.Silent(), true);
    set_optimizer_attribute(oracleModel, "MIPGap", 1e-4);
    set_optimizer_attribute(oracleModel, "TimeLimit", 5);

    para_oracle_bound = abs(currentInfo.f);
    z_rhs = 20 * 10^(ceil(log10(para_oracle_bound)));
    @variable(oracleModel, z  ≥  - z_rhs);
    @variable(oracleModel, x[i = 1:d]);
    @variable(oracleModel, x_sur[g in 1:d, k in keys(x₀[:sur][g])]);
    @variable(oracleModel, y ≤ 0);

    @objective(oracleModel, Min, z);
    oracleInfo = ModelInfo(oracleModel, x, x_sur, y, z);

    nxtModel = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 0
        )
    ); 
    MOI.set(nxtModel, MOI.Silent(), true);
    set_optimizer_attribute(nxtModel, "MIPGap", 1e-3);
    set_optimizer_attribute(nxtModel, "TimeLimit", 5);

    @variable(nxtModel, x1[i = 1:d]);
    @variable(nxtModel, x_sur1[g in 1:d, k in keys(x₀[:sur][g])]);
    @variable(nxtModel, z1 );
    @variable(nxtModel, y1 );
    nxtInfo = ModelInfo(nxtModel, x1, x_sur1, y1, z1);

    Δ = Inf; τₖ = 1; τₘ = .5; μₖ = 1;

    while true
        add_constraint(
            currentInfo, 
            oracleInfo
        );
        optimize!(oracleModel);

        st = termination_status(oracleModel)
        # @info "oracle, $st, grad = $(currentInfo.dG), G = $(currentInfo.G)"
        
        if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
            f_star = JuMP.objective_value(oracleModel)
        else
            return cutInfo
        end
        # formulate alpha model
        result = Δ_model_formulation(functionHistory, f_star, iter, Output = Output);
        previousΔ = Δ;
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
        if round(previousΔ) > round(Δ)
            x₀ = currentInfo.x; τₖ = μₖ * τₖ;
            if cutSelection == "ELC" 
                cutInfo =  [ 
                    - currentInfo.f - 
                    currentInfo.x[:St]' *  S̃[:St] - 
                    sum(sum(currentInfo.x[:sur][g][k] * S̃[:sur][g][k] for k in keys(S̃[:sur][g])) for g in 1:binaryInfo.d),                                                              
                    currentInfo.x
                ];
            elseif cutSelection == "LC"
                cutInfo = [ 
                    - currentInfo.f - 
                    currentInfo.x[:St]' *  Ŝ[:St] -  
                    sum(sum(currentInfo.x[:sur][g][k] * Ŝ[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d),                                                              
                    currentInfo.x
                ];
            elseif cutSelection == "ShrinkageLC"
                cutInfo = [ 
                    currentInfo.S_at_x̂ - 
                    currentInfo.x[:St]' *  Ŝ[:St] -  
                    sum(sum(currentInfo.x[:sur][g][k] * Ŝ[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d),                                                            
                    currentInfo.x
                ]; 
            end 
        else
            τₖ = (τₖ + τₘ) / 2;
        end
        
        ## update α
        if μ/2 ≤ (α-a_min)/(a_max-a_min) .≤ 1-μ/2
            α = α
        else
            α = (a_min+a_max)/2
        end

        ## update level
        w = α * f_star;
        W = minimum( α * functionHistory.f_his[j] + (1-α) * functionHistory.G_max_his[j] for j in 1:iter);

        λ = iter ≤ 10 ? 0.05 : 0.1;
        λ = iter ≥ 20 ? 0.2 : λ;
        λ = iter ≥ 40 ? 0.3 : λ;
        λ = iter ≥ 50 ? 0.4 : λ;
        λ = iter ≥ 60 ? 0.5 : λ;
        λ = iter ≥ 70 ? 0.6 : λ;
        λ = iter ≥ 80 ? 0.7 : λ;
        λ = iter ≥ 90 ? 0.8 : λ;
        
        level = w + λ * (W - w);
        

        ## ==================================================== next iteration point ============================================== ##
        # obtain the next iteration point
        if iter == 1
            @constraint(nxtModel, level_constraint, α * z1 + (1 - α) * y1 ≤ level);
        else 
            delete(nxtModel, nxtModel[:level_constraint]);
            unregister(nxtModel, :level_constraint);
            @constraint(nxtModel, level_constraint, α * z1 + (1 - α) * y1 ≤ level);
        end

        add_constraint(currentInfo, nxtInfo);
        @objective(nxtModel, Min, sum((x1 .- x₀[:St]) .* (x1 .- x₀[:St])) + sum(sum(x_sur1[g, k] * x_sur1[g, k] for k in keys(x₀[:sur][g])) for g in 1:d) );
        optimize!(nxtModel)
        st = termination_status(nxtModel)
        # @info "$st"
        if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
            x_nxt = Dict( :St => JuMP.value.(x1), 
                                            :sur => Dict(g => Dict(k => JuMP.value.(x_sur1[g, k]) for k in keys(x₀[:sur][g])) for g in 1:d)
                                            );
            λₖ = abs(dual(level_constraint)); μₖ = λₖ + 1; 
        else
            return cutInfo
        end

        ## stop rule: gap ≤ .07 * function-value && constraint ≤ 0.05 * LagrangianFunction
        if Δ ≤ abs(f_star_value) * 1e-5 || iter > max_iter
            # @info "yes"
            return cutInfo
        end

        ## ==================================================== end ============================================== ##
        ## save the trajectory
        currentInfo = FuncInfo_LevelSetMethod(
            x_nxt, 
            cutSelection = cutSelection, 
            backwardInfo = backwardInfo, 
            f_star_value = f_star_value, 
            stageData = stageData, 
            Ŝ = Ŝ, 
            S̃ = S̃, 
            binaryInfo = binaryInfo,
            param = param);
        iter = iter + 1;
        
        functionHistory.f_his[iter] = currentInfo.f;
        functionHistory.G_max_his[iter] = maximum(currentInfo.G[k] for k in keys(currentInfo.G));

    end

end
