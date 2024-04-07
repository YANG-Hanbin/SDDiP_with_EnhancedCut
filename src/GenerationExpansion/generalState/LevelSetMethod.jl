#############################################################################################
###########################  auxiliary functions for level set method #######################
#############################################################################################

"""
This function is to constraint the model for solving gap and alpha
"""

function Δ_model_formulation(functionHistory::FunctionHistory, f_star::Float64, iter::Int64; Output::Int64 = 0)
    
    alphaModel = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV),
            "OutputFlag" => Output, 
            "Threads" => 0)
            )

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
function add_constraint(currentInfo::CurrentInfo, modelInfo::ModelInfo)
    m = length(currentInfo.G)

    xⱼ = currentInfo.x
    # add constraints     
    @constraint(modelInfo.model, modelInfo.z .≥ currentInfo.f + currentInfo.df' * (modelInfo.x .- xⱼ) 
                                                                    )

    @constraint(modelInfo.model, [k = 1:m], modelInfo.y .≥ currentInfo.G[k] + sum(currentInfo.dG[k] .* (modelInfo.x .- xⱼ))
                                                                                 ) 
end


"""
    This function is to collect the information from the objecive f, and constraints G
"""
function FuncInfo_LevelSetMethod(x₀::Vector{Float64}; 
                                        backwardInfo::BackwardModelInfo = backwardInfo,
                                            cutSelection::String = cutSelection, 
                                                f_star_value::Union{Float64, Nothing} = f_star_value, 
                                                    stageData::StageData = stageData, 
                                                        Ŝ::Vector{Float64} =  Ŝ, S̃::Union{Vector{Float64}, Nothing} =  S̃, 
                                                            ϵ::Float64 = ϵ )

    if cutSelection == "ELC"
        @objective(backwardInfo.model, Min, stageData.c1' * backwardInfo.x + stageData.c2' * backwardInfo.y + backwardInfo.θ + stageData.penalty * backwardInfo.slack + 
                                                            x₀' * ( Ŝ .- backwardInfo.Sc) );
        optimize!(backwardInfo.model);
        F_solution = ( F = JuMP.objective_value(backwardInfo.model), 
                            ∇F = Ŝ .- JuMP.value.(backwardInfo.Sc) );

        currentInfo  = CurrentInfo(x₀, 
                                    - F_solution.F - x₀' * ( S̃ .-  Ŝ),
                                    Dict(1 => (1 - ϵ) * f_star_value - F_solution.F),
                                    - F_solution.∇F - ( S̃ .-  Ŝ),
                                    Dict(1 => - F_solution.∇F), 
                                    F_solution.F
                                    )                                        
    elseif cutSelection == "LC"
        @objective(backwardInfo.model, Min, stageData.c1' * backwardInfo.x + stageData.c2' * backwardInfo.y + backwardInfo.θ + stageData.penalty * backwardInfo.slack - 
                                                            x₀' * backwardInfo.Sc );
        optimize!(backwardInfo.model);
        F_solution = (F = JuMP.objective_value(backwardInfo.model), 
                                ∇F = - JuMP.value.(backwardInfo.Sc) );

        currentInfo  = CurrentInfo(x₀, 
                                    - F_solution.F - x₀' *  Ŝ, 
                                    Dict(1 => 0.0),
                                    - F_solution.∇F -  Ŝ,
                                    Dict(1 => - F_solution.∇F * 0), 
                                    F_solution.F
                                    );
    elseif cutSelection == "ELCwithoutConstraint"
        @objective(backwardInfo.model, Min, stageData.c1' * backwardInfo.x + stageData.c2' * backwardInfo.y + backwardInfo.θ + stageData.penalty * backwardInfo.slack + 
                                                            x₀' * ( Ŝ .- backwardInfo.Sc) );
        optimize!(backwardInfo.model);
        F_solution = ( F = JuMP.objective_value(backwardInfo.model), 
                        ∇F = Ŝ .- JuMP.value.(backwardInfo.Sc) );

        currentInfo  = CurrentInfo(x₀, 
                                    - F_solution.F - x₀' * ( S̃ .-  Ŝ),
                                    Dict(1 => 0.0 ),
                                    - F_solution.∇F - ( S̃ .-  Ŝ),
                                    Dict(1 => - F_solution.∇F * 0), 
                                    F_solution.F
                                    ) 
    elseif cutSelection == "ShrinkageLC"
        @objective(backwardInfo.model, Min, stageData.c1' * backwardInfo.x + stageData.c2' * backwardInfo.y + backwardInfo.θ + stageData.penalty * backwardInfo.slack + 
                                                            x₀' * ( Ŝ .- backwardInfo.Sc) );
        optimize!(backwardInfo.model);
        F_solution = ( F = JuMP.objective_value(backwardInfo.model), 
                        ∇F = Ŝ .- JuMP.value.(backwardInfo.Sc) );

        currentInfo  = CurrentInfo(x₀, 
                                    1/2 * sum(x₀ .* x₀),
                                    Dict(1 => (1 - ϵ) * f_star_value - F_solution.F),
                                    x₀,
                                    Dict(1 => - F_solution.∇F), 
                                    F_solution.F
                                    ) 

    end

    return currentInfo
end

######################################################################################################################
## -------------------------------------------------- Main Function -------------------------------------------- ##
######################################################################################################################

function LevelSetMethod_optimization!( backwardInfo::BackwardModelInfo, x₀::Vector{Float64}; 
                                        levelSetMethodParam::LevelSetMethodParam = levelSetMethodParam, 
                                        stageData::StageData = stageData, 
                                        ϵ::Float64 = 1e-4, 
                                        binaryInfo::BinaryInfo = binaryInfo) 
    
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
    currentInfo = FuncInfo_LevelSetMethod(x₀, cutSelection = cutSelection, backwardInfo = backwardInfo, f_star_value = f_star_value, stageData = stageData, Ŝ = Ŝ, S̃ = S̃, ϵ = ϵ);

    functionHistory = FunctionHistory(  Dict(1 => currentInfo.f), 
                                            Dict(1 => maximum(currentInfo.G[k] for k in keys(currentInfo.G)) )
                                            );

    ## model for oracle
    oracleModel = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 0)
            );

    para_oracle_bound = abs(currentInfo.f);
    z_rhs = 10 * 10^(ceil(log10(para_oracle_bound)));
    @variable(oracleModel, z  ≥  - z_rhs);
    @variable(oracleModel, x[i = 1:d]);
    @variable(oracleModel, y ≤ 0);

    @objective(oracleModel, Min, z);
    oracleInfo = ModelInfo(oracleModel, x, y, z);

    nxtModel = Model(
        optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
        "OutputFlag" => Output, 
        "Threads" => 0)
        );

    @variable(nxtModel, x1[i = 1:d]);
    @variable(nxtModel, z1 );
    @variable(nxtModel, y1 );
    nxtInfo = ModelInfo(nxtModel, x1, y1, z1);


    Δ = Inf; τₖ = 1; τₘ = .5; μₖ = 1;

    if cutSelection == "ELC" || cutSelection == "ELCwithoutConstraint"
        cutInfo =  [ - currentInfo.f - currentInfo.x' *  S̃,  
                                                                    currentInfo.x] 
    elseif cutSelection == "LC"
        cutInfo = [ - currentInfo.f - currentInfo.x' *  Ŝ,  
                                                                    currentInfo.x] 
    elseif cutSelection == "ShrinkageLC"
        cutInfo = [ currentInfo.S_at_x̂ - currentInfo.x' *  Ŝ,  
                                                                    currentInfo.x] 
    end 

    while true
        add_constraint(currentInfo, oracleInfo);
        optimize!(oracleModel);

        st = termination_status(oracleModel)
        # @info "oracle, $st, grad = $(currentInfo.dG), G = $(currentInfo.G)"
        
        if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
            f_star = JuMP.objective_value(oracleModel)
        else
            @info "oracle is infeasible, new start!"
            oracleModel = Model(
            optimizer_with_attributes(
                ()->Gurobi.Optimizer(GRB_ENV), 
                "OutputFlag" => Output, 
                "Threads" => 0)
                );

            para_oracle_bound = abs(currentInfo.f);
            z_rhs = 10 * 10^(ceil(log10(para_oracle_bound)));
            @variable(oracleModel, z  ≥  - z_rhs);
            @variable(oracleModel, x[i = 1:d]);
            @variable(oracleModel, y ≤ 0);

            @objective(oracleModel, Min, z);
            oracleInfo = ModelInfo(oracleModel, x, y, z);
            add_constraint(currentInfo, oracleInfo);
            optimize!(oracleModel);
            f_star = JuMP.objective_value(oracleModel)
        end

        # f_star = JuMP.objective_value(oracleModel)

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
            if cutSelection == "ELC" || cutSelection == "ELCwithoutConstraint"
                cutInfo =  [ - currentInfo.f - currentInfo.x' *  S̃,  
                                                                            currentInfo.x];
            elseif cutSelection == "LC"
                cutInfo = [ - currentInfo.f - currentInfo.x' *  Ŝ,  
                                                                            currentInfo.x];
            elseif cutSelection == "ShrinkageLC"
                cutInfo = [ currentInfo.S_at_x̂ - currentInfo.x' *  Ŝ,  
                                                                            currentInfo.x];
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
        @objective(nxtModel, Min, sum((x1 .- x₀) .* (x1 .- x₀)));
        optimize!(nxtModel)
        st = termination_status(nxtModel)
        # @info "$st"
        if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
            x_nxt = JuMP.value.(x1);
            λₖ = abs(dual(level_constraint)); μₖ = λₖ + 1; 
        elseif st == MOI.NUMERICAL_ERROR || st == MOI.INFEASIBLE_OR_UNBOUNDED
            # @info "Numerical Error occures! -- Build a new nxtModel"
            nxtModel = Model(
                optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                "OutputFlag" => Output, 
                "Threads" => 0)
                );
        
            @variable(nxtModel, x1[i = 1:d]);
            @variable(nxtModel, z1 );
            @variable(nxtModel, y1 );
            nxtInfo = ModelInfo(nxtModel, x1, y1, z1);
            @constraint(nxtModel, level_constraint, α * z1 + (1 - α) * y1 ≤ level);
            add_constraint(currentInfo, nxtInfo);
            @objective(nxtModel, Min, (x1 - x₀)' * (x1 - x₀));
            optimize!(nxtModel);
            x_nxt = JuMP.value.(x1);
            λₖ = abs(dual(level_constraint)); μₖ = λₖ + 1; 
        else
            set_normalized_rhs( level_constraint, w + .1 * (W - w));
            optimize!(nxtModel);
            x_nxt = JuMP.value.(x1);
        end

        ## stop rule: gap ≤ .07 * function-value && constraint ≤ 0.05 * LagrangianFunction
        if Δ ≤ (abs(f_star_value) * 1e-5 + 5) || iter > max_iter
            # @info "yes"
            return cutInfo
        end

        ## ==================================================== end ============================================== ##
        ## save the trajectory
        currentInfo = FuncInfo_LevelSetMethod(x_nxt, cutSelection = cutSelection, backwardInfo = backwardInfo, f_star_value = f_star_value, stageData = stageData, Ŝ = Ŝ, S̃ = S̃, ϵ = ϵ);
        iter = iter + 1;
        
        functionHistory.f_his[iter] = currentInfo.f;
        functionHistory.G_max_his[iter] = maximum(currentInfo.G[k] for k in keys(currentInfo.G));

    end

end
