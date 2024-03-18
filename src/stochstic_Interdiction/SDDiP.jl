function setupLevelsetPara(f_star_value::Real, x̂::Vector{Float64};
                                    cutSelection::String = "ELC", 
                                    Output_Gap::Bool = false,
                                    λ::Union{Float64, Nothing} = .3, ℓ1::Real = 1.0, ℓ2::Real = 0.8)
    if cutSelection == "ELC"
        x̃ = x̂ .* ℓ2 .+ (1 - ℓ2)/2;

        Output = 0; threshold = 1.0; 
        levelSetMethodParam = LevelSetMethodParam(0.95, λ, threshold, 
                                                            1e14, 2e2, Output, Output_Gap,
                                                            x̂, cutSelection, x̃, f_star_value)


    elseif cutSelection == "ShrinkageLC" 
        Output = 0; threshold = 1.0; 
        levelSetMethodParam = LevelSetMethodParam(0.95, λ, threshold, 
                                                            1e14, 1e2, Output, Output_Gap,
                                                            x̂,  cutSelection, nothing, f_star_value)

    elseif cutSelection == "ELCwithoutConstraint" 
        L̃ = L̂ .* ℓ2 .+ (1 - ℓ2)/2;
        Output = 0; threshold = 1.0; f_star_value = 0.0;
        levelSetMethodParam = LevelSetMethodParam(0.95, λ, threshold, 
                                                            1e14, 1e2, Output, Output_Gap,
                                                            x̂,  cutSelection, x̃, f_star_value)

    elseif cutSelection == "LC" 
        L̃ = nothing; f_star_value = 0.0;
        Output = 0; threshold = 1.0; 
        levelSetMethodParam = LevelSetMethodParam(0.95, λ, threshold, 
                                                            1e14, 1e2, Output, Output_Gap,
                                                            x̂,  cutSelection, x̃, f_star_value)
    end

    x₀ =  x̂ .* f_star_value .* ℓ1 .- f_star_value * (ℓ1 ./ 2);

    return (levelSetMethodParam = levelSetMethodParam, x₀ = x₀)
end


function SDDiP_algorithm( ;stageData::StageData = stageData, 
                            indexSets::IndexSets = indexSets,
                            prob::Prob = prob,
                            Ω::Dict{Int64, ScenarioData} = Ω, 
                            ϵ::Float64 = 0.001, max_iter::Int64 = 200, cutSelection::String = "ELC"
                                )
    ## d: dimension of x
    ## M: num of scenarios when doing one iteration
    initial = now(); iter_time = 0; total_Time = 0; t0 = 0.0;
    T = length(keys(Ω));
    i = 1; LB = - Inf; UB = Inf; incumbent = [0. for l in indexSets.D]; x̂ = [0. for l in indexSets.D]; cutInfo = Dict();

    col_names = [:iter, :LB, :OPT, :UB, :gap, :time, :Time]; # needs to be a vector Symbols
    col_types = [Int64, Float64, Float64, Float64, String, Float64, Float64];
    named_tuple = (; zip(col_names, type[] for type in col_types )...);
    sddipResult = DataFrame(named_tuple); # 0×7 DataFrame
    gapList = [];
    @time gurobiResult = gurobiOptimize!();
    OPT = gurobiResult.OPT 

    forward_stage1_Info = forward_stage1_Model!(); forward_stage2_Info = forward_stage2_Model!(Ω[1], incumbent);
    backward_stage2_Info = backwardModel!(Ω[1], incumbent)

    println("---------------- print out iteration information -------------------")
    while true
        t0 = now();
      
        ## Forward Step
        if true             # the first-stage problem
            optimize!(forward_stage1_Info.model); 
            LB = JuMP.objective_value(forward_stage1_Info.model);
            x̂ = round.(JuMP.value.(forward_stage1_Info.x), digits = 3);
            stage1_Value = JuMP.objective_value(forward_stage1_Info.model) - sum(JuMP.value.(forward_stage1_Info.θ))
        end

        if true             # the second-stage problem
            stage2_Value = Dict();
            for ω in indexSets.Ω
                forward_modify_constraints!(forward_stage2_Info, Ω[ω], x̂); 
                optimize!(forward_stage2_Info.model); 
                stage2_Value[ω] = JuMP.objective_value(forward_stage2_Info.model)

                # obtain dual variables 
                vω = JuMP.objective_value(forward_stage2_Info.model)
                πω = JuMP.dual.(forward_stage2_Info.model[:nonantipativity])
                cutInfo[ω] = LagrangianCut(vω, πω, x̂)
            end
        end

        ## compute the upper bound
        if stage1_Value + sum(prob.p[ω] * stage2_Value[ω] for ω in indexSets.Ω) < UB 
            incumbent = x̂
        end
        UB = minimum([stage1_Value + sum(prob.p[ω] * stage2_Value[ω] for ω in indexSets.Ω), UB])

        gap = round((UB-LB)/UB * 100 ,digits = 2);
        gapString = string(gap,"%");
        push!(sddipResult, [i, LB, OPT, UB, gapString, iter_time, total_Time]); push!(gapList, gap);
        if i == 1
            println("---------------------------------- Iteration Info ------------------------------------")
            println("Iter |   LB                              UB                             gap")
        end
        @printf("%3d  |   %5.3g                         %5.3g                              %1.3f%s\n", i, LB, UB, gap, "%")
        if OPT-LB ≤ 1e-2 * OPT || i > max_iter
            return Dict(:solHistory => sddipResult, 
                            :solution => incumbent, 
                            :gapHistory => gapList) 
        end


        for ω in indexSets.Ω 
            backward_modify_constraints!(backward_stage2_Info, Ω[ω], x̂)
            (levelSetMethodParam, x₀) = setupLevelsetPara(stage2_Value[ω], x̂; cutSelection = cutSelection,  # "ShrinkageLC", "ELCwithoutConstraint", "LC", "ELC"
                                                                        Output_Gap = false, λ = .3, ℓ1 = 0.0, ℓ2 = 0.0);
            (vω, πω) = LevelSetMethod_optimization!(backward_stage2_Info, x₀; ϵ = ϵ, levelSetMethodParam = levelSetMethodParam) 
            cutInfo[ω] = LagrangianCut(vω, πω, nothing)
        end


        ## Add cuts 
        for ω in indexSets.Ω 
            if cutInfo[ω] === nothing
                @constraint(forward_stage1_Info.model, forward_stage1_Info.θ[ω] ≥ cutInfo[ω].v + sum(cutInfo[ω].π .* forward_stage1_Info.x) )
            else
                @constraint(forward_stage1_Info.model, forward_stage1_Info.θ[ω] ≥ cutInfo[ω].v + sum(cutInfo[ω].π .* (forward_stage1_Info.x .- cutInfo[ω].x)))
            end
        end


        t1 = now();
        iter_time = (t1 - t0).value/1000;
        total_Time = (t1 - initial).value/1000;
        
        i = i + 1;

    end

end

# using DataFrames
# using Latexify
# df = DataFrame(A = 'x':'z', B = ["M", "F", "F"])
# latexify(df; env=:table, latex=false)

# save("src/GenerationExpansion/3.0version/testData_5/sddipResult_LCinterval.jld2", "sddipResult", sddipResult)
# save("src/GenerationExpansion/3.0version/testData_5/sddipResult_LCbinary.jld2", "sddipResult", sddipResult)
# save("src/GenerationExpansion/3.0version/testData_5/sddipResult_ELCinterval0.jld2", "sddipResult", sddipResult)
# save("src/GenerationExpansion/3.0version/testData_5/sddipResult_ELCbinary0.jld2", "sddipResult", sddipResult)



## ------------------------- review --------------------- ##
# 1. if we set the interior point is close the center, when need lambda to be close to 1
# 2. ...