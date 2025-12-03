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
    currentInfo::NormalizedCurrentInfo, 
    modelInfo::ModelInfo
)::Nothing
    m = length(currentInfo.G);
    xⱼ = currentInfo.x;
    # add constraints     
    @constraint(
        modelInfo.model, 
        modelInfo.z .≥ currentInfo.f + 
        currentInfo.df[:St]' * (modelInfo.x .- xⱼ[1][:St]) + 
        sum(
            sum(
                currentInfo.df[:sur][g][k] * (modelInfo.x_sur[g, k] .- xⱼ[1][:sur][g][k]) for k in keys(currentInfo.df[:sur][g])
                ) for g in 1:length(modelInfo.x)
        ) + 
        currentInfo.df[:obj] * (modelInfo.model[:xθ] - xⱼ[2])                                                        
    );

    @constraint(
        modelInfo.model, 
        [k = 1:1], 
        modelInfo.y .≥ currentInfo.G[k] + 
        sum(currentInfo.dG[k][1][:St] .* (modelInfo.x .- xⱼ[1][:St])) + 
        sum(
            sum(
                currentInfo.dG[k][1][:sur][g][km] * (modelInfo.x_sur[g, km] .- xⱼ[1][:sur][g][km]) for km in keys(currentInfo.dG[k][1][:sur][g])
            ) for g in 1:length(modelInfo.x)
        )  + 
        currentInfo.dG[k][2] * (modelInfo.model[:xθ] - currentInfo.x[2])                                                                     
    ); 

    @constraint(
        modelInfo.model, 
        [k = 2:2], 
        modelInfo.y .≥ currentInfo.G[k] +  
        currentInfo.dG[k][2] * (modelInfo.model[:xθ] - currentInfo.x[2])                                                                     
    ); 
    return 
end


"""
    This function is to collect the information from the objecive f, and constraints G
"""
function FuncInfo_LevelSetMethod(
    πₙ₀::Float64, 
    πₙ::Dict{Symbol, Any}; 
    backwardInfo::BackwardModelInfo = backwardInfo,
    f_star_value::Union{Float64, Nothing} = f_star_value, 
    Ŝ::Dict{Symbol, Any} =  Ŝ, 
    S̃::Dict{Symbol, Any} =  S̃, 
    binaryInfo::BinaryInfo = binaryInfo
)::NormalizedCurrentInfo

    @objective(
        backwardInfo.model, 
        Min, 
        πₙ₀ * (f_star_value - backwardInfo.model[:primal_objective_expression]) + 
        πₙ[:St]' * (Ŝ[:St] .- backwardInfo.Sc) + 
        sum(
            sum(
                πₙ[:sur][g][k] * (Ŝ[:sur][g][k] - backwardInfo.model[:sur_copy][g, k]) for k in keys(πₙ[:sur][g])
            ) for g in 1:binaryInfo.d
        )
    );
    optimize!(backwardInfo.model);
    uₙ = Dict(
        :St => S̃[:St], 
        :obj => 1.0
    );
    normalization_function = πₙ₀ * uₙ[:obj] + πₙ[:St]' * uₙ[:St];

    currentInfo  = NormalizedCurrentInfo(
        (πₙ, πₙ₀), 
        - JuMP.objective_value(backwardInfo.model),
        Dict(
            1 => normalization_function - 1, 
            2 => πₙ₀
        ),
        Dict( 
            :obj => JuMP.value(backwardInfo.model[:primal_objective_expression]) - f_star_value,
            :St => JuMP.value.(backwardInfo.Sc) .- Ŝ[:St], 
            :sur => Dict(
                g => Dict(
                    k => JuMP.value(backwardInfo.model[:sur_copy][g, k]) - Ŝ[:sur][g][k] for k in keys(Ŝ[:sur][g])
                ) for g in 1:binaryInfo.d
            )
        ),
        Dict(
            1 => (
                Dict(
                :St => uₙ[:St], 
                :sur => Dict(
                    g => Dict(
                        k => 0.0 for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d
                    )
                ), uₙ[:obj]
            ), 
            2 => (nothing, 1.0)
        ),
        JuMP.objective_value(backwardInfo.model) 
    );
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
    L̂ = stateInfo.stageSolution; 
    sur = stateInfo.stageSur; 
    Ŝ = Dict(
        :St => L̂, 
        :sur => Dict(
            g => Dict(k => sur[g][k] for k in keys(sur[g])
        ) for g in 1:binaryInfo.d)
    );
    forward_modify_constraints!(
        forwardInfo,                   
        stageData,                
        demand,                  
        L̂,                     
        binaryInfo = binaryInfo                         
    );
    optimize!(forwardInfo.model); 
    f_star_value = JuMP.objective_value(forwardInfo.model);
    S̃ = Dict( 
        :St => L̂ .* 0 .+ .5,    
        :sur => Dict(
            g => Dict(
                k => .5 for k in keys(sur[g])
            ) for g in 1:binaryInfo.d
        )
    );
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
        S̃, 
        f_star_value
    );

    πₙ = Dict(
        :St => L̂ .* 1.0,                 
        :sur => Dict(
            g => Dict(
                k => 0.0 for k in keys(sur[g])
            ) for g in 1:binaryInfo.d
        )
    );
    return (
        levelSetMethodParam = levelSetMethodParam, 
        πₙ = πₙ
    )
end

######################################################################################################################
## -------------------------------------------------- Main Function -------------------------------------------- ##
######################################################################################################################
function LevelSetMethod_optimization!( 
    backwardInfo::BackwardModelInfo, 
    πₙ::Dict{Symbol, Any}; 
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
        - 10.0, 
        πₙ; 
        backwardInfo = backwardInfo,
        f_star_value = f_star_value, 
        Ŝ = Ŝ, 
        S̃ = S̃, 
        binaryInfo = binaryInfo
    );

    cutInfo = [ 
        - currentInfo.f - 
        currentInfo.x[2] * f_star_value - 
        currentInfo.x[1][:St]' *  Ŝ[:St] -  
        sum(sum(currentInfo.x[1][:sur][g][k] * Ŝ[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d),                                                              
        currentInfo.x[1],
         - currentInfo.x[2]
    ];

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
    @variable(oracleModel, xθ);
    @variable(oracleModel, x_sur[g in 1:d, k in keys(πₙ[:sur][g])]);
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
    @variable(nxtModel, xθ);
    @variable(nxtModel, x_sur1[g in 1:d, k in keys(πₙ[:sur][g])]);
    @variable(nxtModel, z1);
    @variable(nxtModel, y1);
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
        πₙ = currentInfo.x[1]; πₙ₀ = currentInfo.x[2];
        if round(previousΔ) > round(Δ)
            πₙ = currentInfo.x[1]; πₙ₀ = currentInfo.x[2]; τₖ = μₖ * τₖ;
            cutInfo = [ 
                - currentInfo.f - 
                currentInfo.x[2] * f_star_value - 
                currentInfo.x[1][:St]' *  Ŝ[:St] -  
                sum(sum(currentInfo.x[1][:sur][g][k] * Ŝ[:sur][g][k] for k in keys(Ŝ[:sur][g])) for g in 1:binaryInfo.d),                                                              
                currentInfo.x[1],
                - currentInfo.x[2]
            ]; 
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
        @objective(
            nxtModel, 
            Min, 
            (nxtModel[:xθ] - πₙ₀)^2 +
            sum(
                (x1 .- πₙ[:St]) .* (x1 .- πₙ[:St])) + 
                sum(sum((x_sur1[g, k] - πₙ[:sur][g][k]) * (x_sur1[g, k] - πₙ[:sur][g][k]) for k in keys(πₙ[:sur][g])) for g in 1:d
            ) 
        );
        optimize!(nxtModel)
        st = termination_status(nxtModel)
        # @info "$st"
        if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
            πₙ = Dict(
                :St => JuMP.value.(x1), 
                :sur => Dict(g => Dict(k => JuMP.value.(x_sur1[g, k]) for k in keys(πₙ[:sur][g])) for g in 1:d)            
            );
            πₙ₀ = JuMP.value(xθ) ≤ 0 ? JuMP.value(xθ) : 0.0;
            # πₙ₀ = JuMP.value(xθ);
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
            πₙ₀, 
            πₙ; 
            backwardInfo = backwardInfo,
            f_star_value = f_star_value, 
            Ŝ = Ŝ, 
            S̃ = S̃, 
            binaryInfo = binaryInfo
        );
        iter = iter + 1;
        
        functionHistory.f_his[iter] = currentInfo.f;
        functionHistory.G_max_his[iter] = maximum(currentInfo.G[k] for k in keys(currentInfo.G));

    end

end
