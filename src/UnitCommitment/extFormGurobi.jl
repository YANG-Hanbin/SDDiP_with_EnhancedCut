################################################################################################################################################
############################################################     Gurobi function   #############################################################
################################################################################################################################################
# https://www.gurobi.com/documentation/9.5/refman/logging.html
function gurobiOptimize!(indexSets::IndexSets, 
                            paramDemand::ParamDemand, 
                            paramOPF::ParamOPF, Ξ::Dict; timelimit::Int64 = 6000, 
                            mipGap::Float64 = 5e-3, β::Float64 = .5)  

    (D, G, L, B, T) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T);
    (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L);

    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => 1, 
                                          "Threads" =>0, 
                                          "MIPGap" => mipGap, 
                                          "TimeLimit" => timelimit) 
                                          );
    ## the first stage variables
    @variable(model, θ_angle[B, 1:T]);      ## phase angle of the bus i
    @variable(model, P[L, 1:T]);            ## real power flow on line l; elements in L is Tuple (i, j)
    @variable(model, s[G, 1:T] ≥ 0);        ## real power generation at generator g
    @variable(model, 0 ≤ x[D, 1:T] ≤ 1);    ## load shedding


    @variable(model, zg[G, 1:T], Bin);      ## binary status indicator of generator g
    @variable(model, zb[B, 1:T], Bin);      ## binary status indicator of bus i
    @variable(model, zl[L, 1:T], Bin);      ## binary status indicator of line l

    @variable(model, rg[G, 1:T], Bin);      ## binary restoration indicator of generator g
    @variable(model, rb[B, 1:T], Bin);      ## binary restoration indicator of bus i
    @variable(model, rl[L, 1:T], Bin);      ## binary restoration indicator of line l

    @variable(model, S[D] ≥ 0);             ## cumulative load-shedding
    @variable(model, maxS ≥ 0)      
    @variable(model, minS ≥ 0)  


    ## the subproblems' variables
    @variable(model, θω[B, 1:T, keys(Ξ)]);        ## phase angle of the bus i
    @variable(model, Pω[L, 1:T, keys(Ξ)]);        ## real power flow on line l; elements in L is Tuple (i, j)  
    @variable(model, sω[G, 1:T, keys(Ξ)]);        ## real power generation at generator g
    @variable(model, 0 ≤ xω[D, 1:T, keys(Ξ)] ≤ 1);## load shedding

    @variable(model, zbω[B, 1:T, keys(Ξ)], Bin);
    @variable(model, zgω[G, 1:T, keys(Ξ)], Bin);
    @variable(model, zlω[L, 1:T, keys(Ξ)], Bin);

    @variable(model, νb[B, keys(Ξ)], Bin);
    @variable(model, νg[G, keys(Ξ)], Bin);
    @variable(model, νl[L, keys(Ξ)], Bin);

    # @variable(model, copyb[B, keys(Ξ)], Bin)
    # @variable(model, copyg[G, keys(Ξ)], Bin)
    # @variable(model, copyl[L, keys(Ξ)], Bin) 

    # constraint 1 a b
    for l in L
      i = l[1];
      j = l[2];
      @constraint(model, [t in 1:T], P[l, t] ≤ - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmax * (1 - zl[l, t] ) ) );
      @constraint(model, [t in 1:T], P[l, t] ≥ - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmin * (1 - zl[l, t] ) ) );
    end

    ## constraint 1c
    @constraint(model, [l in L, t in 1:T], P[l, t] ≥ - paramOPF.W[l] * zl[l, t] );
    @constraint(model, [l in L, t in 1:T], P[l, t] ≤ paramOPF.W[l] * zl[l, t] );

    ## constraint 1d
    @constraint(model, [i in B, t in 1:T], sum(s[g, t] for g in Gᵢ[i]) + 
                                              sum(P[(i, j), t] for j in out_L[i] ) - 
                                                sum(P[(j, i), t] for j in in_L[i] ) 
                                                  .== sum(paramDemand.demand[t][d] * x[d, t] for d in Dᵢ[i]) );
    
    ## constraint 1e
    @constraint(model, [g in G, t in 1:T], s[g, t] ≥ paramOPF.smin[g] * zg[g, t] );
    @constraint(model, [g in G, t in 1:T], s[g, t] ≤ paramOPF.smax[g] * zg[g, t] );

    ## constraint 1f g h i
    @constraint(model, [i in B, t in 1:T, d in Dᵢ[i]], zb[i, t] ≥ x[d, t] );
    @constraint(model, [i in B, t in 1:T, g in Gᵢ[i]], zb[i, t] ≥ zg[g, t]);
    @constraint(model, [i in B, t in 1:T, j in out_L[i]], zb[i, t] ≥ zl[(i, j), t] );
    @constraint(model, [i in B, t in 1:T, j in in_L[i]], zb[i, t] ≥ zl[(j, i), t] );

    ## constraint 1j k
    @constraint(model, [i in B, t in 1:T-1], rb[i, t] - rb[i, t+1] ≥ 0);
    @constraint(model, [g in G, t in 1:T-1], rg[g, t] - rg[g, t+1] ≥ 0);
    @constraint(model, [l in L, t in 1:T-1], rl[l, t] - rl[l, t+1] ≥ 0);
    @constraint(model, [i in B, t in 1:T-1], rb[i, t] - rb[i, t+1] ≥ zb[i, t+1] - zb[i, t]);
    @constraint(model, [g in G, t in 1:T-1], rg[g, t] - rg[g, t+1] ≥ zg[g, t+1] - zg[g, t]);
    @constraint(model, [l in L, t in 1:T-1], rl[l, t] - rl[l, t+1] ≥ zl[l, t+1] - zl[l, t]);

    ## constraint 1l m: fairness
    @constraint(model, fairness[d1 in D, d2 in D], sum((x[d1, t] - x[d2, t]) for t in 1:T) ≤ β * T)

    # @constraint(model, [d in D], S[d] == sum((1 - x[d, t]) for t in 1:T));
    # @constraint(model, [d in D], maxS ≥ S[d]);
    # @constraint(model, [d in D], minS ≤ S[d]);
    # @constraint(model, fairness, β * maxS ≥ maxS - minS);

    ex = @expression(model, 0);
    ## second stage constraints
    for ω in keys(Ξ) 
        ξ = Ξ[ω]; rv = ξ[maximum([1, length(ξ)-1])].rv; 
        p = ξ[length(ξ)].p; 
        ex = @expression(model, ex + p * (sum( sum(paramDemand.w[d] * (1 - xω[d, t, ω]) for d in D ) for t in 1:T) +
                                                    sum(paramDemand.cb[i] * νb[i, ω] for i in B) + 
                                                        sum(paramDemand.cg[g] * νg[g, ω] for g in G) + 
                                                            sum(paramDemand.cl[l] * νl[l, ω] for l in L)
                                                        ) );

        ## constraints 2 a b
        @constraint(model, [l in L, t in minimum([rv.τ, T]):T], Pω[l, t, ω] ≤ - paramOPF.b[l] * (θω[l[1], t, ω] - θω[l[2], t, ω] + paramOPF.θmax * (1 - zlω[l, t, ω]) ));
        @constraint(model, [l in L, t in minimum([rv.τ, T]):T], Pω[l, t, ω] ≥ - paramOPF.b[l] * (θω[l[1], t, ω] - θω[l[2], t, ω] + paramOPF.θmin * (1 - zlω[l, t, ω]) ));

        ## constraint 2c
        @constraint(model, [l in L, t in minimum([rv.τ, T]):T], Pω[l, t, ω] ≥ - paramOPF.W[l] * zlω[l, t, ω] );
        @constraint(model, [l in L, t in minimum([rv.τ, T]):T], Pω[l, t, ω] ≤ paramOPF.W[l] * zlω[l, t, ω] );

        ## constraint 2d
        @constraint(model, [i in B, t in minimum([rv.τ, T]):T], sum(sω[g, t, ω] for g in Gᵢ[i]) + 
                                                            sum(Pω[(i, j), t, ω] for j in out_L[i] ) - 
                                                                sum(Pω[(j, i), t, ω] for j in in_L[i] ) 
                                                                    .== sum(paramDemand.demand[t][d] * xω[d, t, ω] for d in Dᵢ[i]) );

        ## constraint 2e
        @constraint(model, [g in G, t in minimum([rv.τ, T]):T], sω[g, t, ω] ≥ paramOPF.smin[g] * zgω[g, t, ω]);
        @constraint(model, [g in G, t in minimum([rv.τ, T]):T], sω[g, t, ω] ≤ paramOPF.smax[g] * zgω[g, t, ω]);

        ## constraint 2f g h i 
        @constraint(model, [i in B, t in minimum([rv.τ, T]):T, d in Dᵢ[i]], zbω[i, t, ω] ≥ xω[d, t, ω]);
        @constraint(model, [i in B, t in minimum([rv.τ, T]):T, g in Gᵢ[i]], zbω[i, t, ω] ≥ zgω[g, t, ω]);
        @constraint(model, [i in B, t in minimum([rv.τ, T]):T, j in out_L[i]], zbω[i, t, ω] ≥ zlω[(i, j), t, ω]);
        @constraint(model, [i in B, t in minimum([rv.τ, T]):T, j in in_L[i]], zbω[i, t, ω] ≥ zlω[(j, i), t, ω]);

        ## constriant 2j 
        @constraint(model, [i in B, t in minimum([rv.τ, T])-1:T-1], zbω[i, t, ω] ≥ zbω[i, t+1, ω]);
        @constraint(model, [g in G, t in minimum([rv.τ, T])-1:T-1], zgω[g, t, ω] ≥ zgω[g, t+1, ω]);
        @constraint(model, [l in L, t in minimum([rv.τ, T])-1:T-1], zlω[l, t, ω] ≥ zlω[l, t+1, ω]);


        # constraint 2k
        @constraint(model, [i in B, t in minimum([rv.τ, T]):T], zbω[i, t, ω] ≤ 1- νb[i, ω] );
        @constraint(model, [g in G, t in minimum([rv.τ, T]):T], zgω[g, t, ω] ≤ 1- νg[g, ω] );
        @constraint(model, [l in L, t in minimum([rv.τ, T]):T], zlω[l, t, ω] ≤ 1- νl[l, ω] );

        # constraint 2l
        @constraint(model, [i in B], νb[i, ω] ≥ rv.vb[i] );
        @constraint(model, [g in G], νg[g, ω] ≥ rv.vg[g] );
        @constraint(model, [l in L], νl[l, ω] ≥ rv.vl[l] );

        ## constraint 2m
        @constraint(model, [i in B, j in unique(rv.Ibb[i])], νb[j, ω] ≥ rv.ub[i] * zbω[i, rv.τ - 1, ω] );
        @constraint(model, [i in B, j in unique(rv.Ibg[i])], νg[j, ω] ≥ rv.ub[i] * zbω[i, rv.τ - 1, ω] );
        @constraint(model, [i in B, j in unique(rv.Ibl[i])], νl[j, ω] ≥ rv.ub[i] * zbω[i, rv.τ - 1, ω] );

        @constraint(model, [i in G, j in unique(rv.Igb[i])], νb[j, ω] ≥ rv.ug[i] * zgω[i, rv.τ - 1, ω] );
        @constraint(model, [i in G, j in unique(rv.Igg[i])], νg[j, ω] ≥ rv.ug[i] * zgω[i, rv.τ - 1, ω] );
        @constraint(model, [i in G, j in unique(rv.Igl[i])], νl[j, ω] ≥ rv.ug[i] * zgω[i, rv.τ - 1, ω] );

        @constraint(model, [i in L, j in unique(rv.Ilb[i])], νb[j, ω] ≥ rv.ul[i] * zlω[i, rv.τ - 1, ω] );
        @constraint(model, [i in L, j in unique(rv.Ilg[i])], νg[j, ω] ≥ rv.ul[i] * zlω[i, rv.τ - 1, ω] );
        @constraint(model, [i in L, j in unique(rv.Ill[i])], νl[j, ω] ≥ rv.ul[i] * zlω[i, rv.τ - 1, ω] );
    end 
    @objective(model, Min, ex);


    # scenarioPath[ω1][ω2] is a list that gives the common realization in scenarios ω1 and ω2
    scenarioPath = Dict(); 
    for ω in keys(Ξ)
        ξ = Ξ[ω]; l1 = length(ξ); scenarioPath[ω] = Dict();
        for j in keys(Ξ) 
            if j != ω
                scenarioPath[ω][j] = []; _ξ = Ξ[j]; l2 = length(_ξ);
                for i in 1:minimum([l1, l2])
                    if ξ[i] == _ξ[i]
                        push!(scenarioPath[ω][j], i)
                    end
                end
            end
        end
    end


    # non-anticipativity
    for ω in keys(Ξ)
        ξ = Ξ[ω]; τ0 = ξ[1].rv.τ; l1 = length(ξ);
        @constraint(model, [i in B, t in 1:τ0-1], θω[i, t, ω] == θ_angle[i, t]);
        @constraint(model, [l in L, t in 1:τ0-1], Pω[l, t, ω] == P[l, t]);
        @constraint(model, [g in G, t in 1:τ0-1], sω[g, t, ω] == s[g, t]);
        @constraint(model, [d in D, t in 1:τ0-1], xω[d, t, ω] == x[d, t]);
        @constraint(model, [i in B], zbω[i, τ0-1, ω] == zb[i, τ0-1]);
        @constraint(model, [l in L], zlω[l, τ0-1, ω] == zl[l, τ0-1]);
        @constraint(model, [g in G], zgω[g, τ0-1, ω] == zg[g, τ0-1]);

        commonset = scenarioPath[ω];
        for j in keys(commonset)
            l = length(commonset[j]);
            _ξ = Ξ[j]; l2 = length(_ξ); 
            if l > 1
                for i in 1:l 
                    τ1 = Ξ[j][i].rv.τ; τ2 = minimum([Ξ[ω][i+1].rv.τ, Ξ[j][i+1].rv.τ]); 
                    @constraint(model, [i in B, t in τ1:τ2-1], θω[i, t, ω] == θω[i, t, j]);
                    @constraint(model, [l in L, t in τ1:τ2-1], Pω[l, t, ω] == Pω[l, t, j]);
                    @constraint(model, [g in G, t in τ1:τ2-1], sω[g, t, ω] == sω[g, t, j]);
                    @constraint(model, [d in D, t in τ1:τ2-1], xω[d, t, ω] == xω[d, t, j]);
                    @constraint(model, [i in B, t in τ1:τ2-1], zbω[i, t, ω] == zbω[i, t, j]);
                    @constraint(model, [l in L, t in τ1:τ2-1], zlω[l, t, ω] == zlω[l, t, j]);
                    @constraint(model, [g in G, t in τ1:τ2-1], zgω[g, t, ω] == zgω[g, t, j]);
                end
            end
        end
    end


    for ω in keys(Ξ)  
        ξ = Ξ[ω]; l1 = length(ξ); commonset = scenarioPath[ω];
        for j in keys(commonset)
            _ξ = Ξ[j]; l2 = length(_ξ); 
            if (l1 < l2 && (commonset[j] != []))
                @constraint(model, [i in B], νb[i, ω] ≤ νb[i, j]);
                @constraint(model, [g in G], νg[g, ω] ≤ νg[g, j]);
                @constraint(model, [l in L], νl[l, ω] ≤ νl[l, j]);
            end
        end
    end

                                                                                     
    ####################################################### solve the model and display the result ###########################################################
    optimize!(model);
    first_stage_variable = Dict(:zg => round.(JuMP.value.(zg)), 
                                    :zb => round.(JuMP.value.(zb)), 
                                    :zl => round.(JuMP.value.(zl))
                                    );
    subproblem_stage_variable = Dict(:zg => round.(JuMP.value.(zgω)), 
                                    :zb => round.(JuMP.value.(zbω)), 
                                    :zl => round.(JuMP.value.(zlω))
                                    );
    first_stage_load_shed = JuMP.value.(x);
    subproblem_load_shed = JuMP.value.(xω);
                                    
    return (OPT = JuMP.objective_value(model), 
                first_stage_variable = first_stage_variable, 
                    subproblem_stage_variable = subproblem_stage_variable,
                        first_stage_load_shed = first_stage_load_shed, 
                            subproblem_load_shed = subproblem_load_shed
                        )
end
