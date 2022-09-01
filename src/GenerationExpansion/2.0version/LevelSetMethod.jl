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



######################################################################################################################
## -------------------------------------------------- Main Function -------------------------------------------- ##
######################################################################################################################

function LevelSetMethod_optimization!( backwardInfo::BackwardModelInfo, f_star_value::Float64; 
                                        levelSetMethodParam::LevelSetMethodParam = levelSetMethodParam, 
                                        stageData::StageData = stageData, 
                                        sum_generator::Vector{Float64} = sum_generator, 
                                        ϵ::Float64 = 1e-4, 
                                        Enhanced_Cut = true, 
                                        binaryInfo::BinaryInfo = binaryInfo) 
    
    ######################################################################################################################
    ###############################   auxiliary function for function information   ######################################
    ######################################################################################################################
    ##  μ larger is better
    (μ, λ, threshold, nxt_bound, max_iter, Output, Output_Gap) = (levelSetMethodParam.μ, levelSetMethodParam.λ, levelSetMethodParam.threshold, levelSetMethodParam.nxt_bound, levelSetMethodParam.max_iter, levelSetMethodParam.Output,levelSetMethodParam.Output_Gap);
    (A, n, d) = (binaryInfo.A, binaryInfo.n, binaryInfo.d);
    # l_interior = [.8 for i in 1:n] - .5 * sum_generator;
    l_interior = sum_generator .* .8 .+ .1;

    ## collect the information from the objecive f, and constraints G
    function compute_f_G(x₀::Vector{Float64}; 
                                            Enhanced_Cut::Bool = true, 
                                            f_star_value::Float64 = f_star_value, 
                                            stageData::StageData = stageData, 
                                            sum_generator::Vector{Float64} = sum_generator, 
                                            ϵ::Float64 = ϵ )

        if Enhanced_Cut 
            @objective(backwardInfo.model, Min, stageData.c1' * backwardInfo.x + stageData.c2' * backwardInfo.y + backwardInfo.θ + stageData.penalty * backwardInfo.slack + 
                                                                x₀' * (sum_generator .- backwardInfo.Lc) );
            optimize!(backwardInfo.model);
            F_solution = [ JuMP.objective_value(backwardInfo.model), sum_generator .- round.(JuMP.value.(backwardInfo.Lc)) ];

            currentInfo  = CurrentInfo(x₀, 
                                        - F_solution[1] - x₀' * (l_interior .- sum_generator),
                                        Dict(1 => (1 - ϵ) * f_star_value - F_solution[1]),
                                        round.(- F_solution[2] - (l_interior .- sum_generator), digits = 2),
                                        Dict(1 => - F_solution[2])

                                        )                                        
        else
            @objective(backwardInfo.model, Min, stageData.c1' * backwardInfo.x + stageData.c2' * backwardInfo.y + backwardInfo.θ + stageData.penalty * backwardInfo.slack - 
                                                                x₀' * backwardInfo.Lc );
            optimize!(backwardInfo.model);
            F_solution = [ JuMP.objective_value(backwardInfo.model), - round.(JuMP.value.(backwardInfo.Lc)) ];

            currentInfo  = CurrentInfo(x₀, 
                                        - F_solution[1] - x₀' * sum_generator, 
                                        Dict(1 => 0.0),
                                        - F_solution[2] - sum_generator,
                                        Dict(1 => - F_solution[2] * 0)
                                        );
        end

        return currentInfo
    end


    ## ==================================================== Levelset Method ============================================== ##
    if Enhanced_Cut
        x₀ = sum_generator .* f_star_value .- .5 .* f_star_value;
    else
        x₀ = zeros(n);
    end 
    # x₀ = zeros(n);

    iter = 1;
    α = 1/2;

    ## trajectory
    currentInfo = compute_f_G(x₀, Enhanced_Cut = Enhanced_Cut);
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
    z_rhs = 1 * 10^(ceil(log10(para_oracle_bound)));
    @variable(oracleModel, z  ≥  - z_rhs);
    @variable(oracleModel, x[i = 1:n]);
    @variable(oracleModel, y ≤ 0);

    @objective(oracleModel, Min, z);
    oracleInfo = ModelInfo(oracleModel, x, y, z);

    nxtModel = Model(
        optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
        "OutputFlag" => Output, 
        "Threads" => 0)
        );

    @variable(nxtModel, x1[i = 1:n]);
    @variable(nxtModel, z1 );
    @variable(nxtModel, y1 );
    nxtInfo = ModelInfo(nxtModel, x1, y1, z1);


    Δ = Inf; τₖ = 1; τₘ = .5; μₖ = 1;

    if Enhanced_Cut
        cutInfo =  [ - currentInfo.f - currentInfo.x' * l_interior,  
                                                                    currentInfo.x] 
    else
        cutInfo = [ - currentInfo.f - currentInfo.x[:zb]' * sum_generator,  
                                                                    currentInfo.x] 
    end 

    while true
        add_constraint(currentInfo, oracleInfo);
        optimize!(oracleModel);

        st = termination_status(oracleModel)
        if st != MOI.OPTIMAL
            @info "oracle is infeasible"
            break
        end

        f_star = JuMP.objective_value(oracleModel)

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
        if round(previousΔ) > round(Δ)
            x₀ = currentInfo.x; τₖ = μₖ * τₖ;
            if Enhanced_Cut
                cutInfo =  [ - currentInfo.f - currentInfo.x' * l_interior,  
                                                                            currentInfo.x] 
            else
                cutInfo = [ - currentInfo.f - currentInfo.x[:zb]' * sum_generator,  
                                                                            currentInfo.x] 
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

        # λ = iter ≤ 8 ? 0.05 : 0.15;
        # λ = iter ≥ 15 ? 0.25 : λ;
        # λ = iter ≥ 25 ? 0.4 : λ;
        # λ = iter ≥ 35 ? 0.6 : λ;
        # λ = iter ≥ 40 ? 0.7 : λ;
        # λ = iter ≥ 45 ? 0.8 : λ;
        
        level = w + λ * (W - w);
        

        ## ==================================================== next iteration point ============================================== ##
        # obtain the next iteration point
        if iter == 1
            @constraint(nxtModel, level_constraint, α * z1 + (1 - α) * y1 ≤ level);
        else 
            delete(nxtModel, nxtModel[:level_constraint])
            unregister(nxtModel, :level_constraint)
            @constraint(nxtModel, level_constraint, α * z1 + (1 - α) * y1 ≤ level);
        end

        add_constraint(currentInfo, nxtInfo);
        @objective(nxtModel, Min, sum((x1 .- x₀) .* (x1 .- x₀)))
        optimize!(nxtModel)
        st = termination_status(nxtModel)
        if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
            x_nxt = JuMP.value.(x1);
            λₖ = abs(dual(level_constraint)); μₖ = λₖ + 1; 
        elseif st == MOI.NUMERICAL_ERROR ## need to figure out why this case happened and fix it
            @info "Numerical Error occures! -- Build a new nxtModel"

            nxtModel = Model(
                optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                "OutputFlag" => Output, 
                "Threads" => 0)
                );
        
            @variable(nxtModel, x1[i = 1:n]);
            @variable(nxtModel, z1 );
            @variable(nxtModel, y1 );
            nxtInfo = ModelInfo(nxtModel, x1, y1, z1);

            @constraint(nxtModel, level_constraint, α * z1 + (1 - α) * y1 ≤ level);
            add_constraint(currentInfo, nxtInfo);
            @objective(nxtModel, Min, (x1 - x₀)' * (x1 - x₀))
            optimize!(nxtModel)
            x_nxt = JuMP.value.(x1)
            λₖ = abs(dual(level_constraint)); μₖ = λₖ + 1; 
        else
            set_normalized_rhs( level_constraint, w + .99 * (W - w))
            optimize!(nxtModel)
            x_nxt = JuMP.value.(x1)
        end

        ## stop rule: gap ≤ .07 * function-value && constraint ≤ 0.05 * LagrangianFunction
        if ( Δ ≤ 100 && currentInfo.G[1] ≤ threshold ) || iter > max_iter
            return cutInfo
        end

        ## ==================================================== end ============================================== ##
        ## save the trajectory
        currentInfo = compute_f_G(x_nxt, Enhanced_Cut = Enhanced_Cut)
        iter = iter + 1
        
        functionHistory.f_his[iter] = currentInfo.f;
        functionHistory.G_max_his[iter] = maximum(currentInfo.G[k] for k in keys(currentInfo.G));

    end

end
