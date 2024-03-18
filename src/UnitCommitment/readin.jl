"""
prepareIndexSets(T::Int64 = 10, network_data::Dict{String, Any} = network_data, branchInfo::DataFrame = branchInfo)

Readin the data and prepare the index sets for the power network
"""
function prepareIndexSets(  ; T::Int64 = 10, 
                                network_data::Dict{String, Any} = network_data,
                                    branchInfo::DataFrame = branchInfo
                                                    )
    D = Vector{Int64}()
    G = Vector{Int64}()
    B = Vector{Int64}()
    L = Vector{Tuple{Int64, Int64}}() 

    Dᵢ =    Dict{Int64,Vector{Int64}}()
    Gᵢ =    Dict{Int64,Vector{Int64}}()
    out_L = Dict{Int64,Vector{Int64}}()
    in_L =  Dict{Int64,Vector{Int64}}()


    _b = Dict{Tuple{Int64, Int64}, Float64}()                   ## total line charging susceptance
    θmax = 100
    θmin = - 100
    W = Dict{Tuple{Int64, Int64}, Float64}()
    smax = Dict{Int64, Float64}()
    smin = Dict{Int64, Float64}()
    M = Dict{Int64, Float64}()
    slope = Dict{Int64, Float64}()
    intercept = Dict{Int64, Float64}()
    C_start = Dict{Int64, Float64}()
    C_down = Dict{Int64, Float64}()

    demand = Dict{Int64, Float64}()

    w = Dict{Int64, Float64}()                                  ## priority level of load D
    cb = Dict{Int64, Float64}()                                 ## set of fire damage cost cᵢ at :b ∈ B
    cg = Dict{Int64, Float64}()
    cl = Dict{Tuple{Int64, Int64}, Float64}()

    for i in keys(network_data["bus"])
        b = network_data["bus"][i]["bus_i"]
        push!(B, b)
        Dᵢ[b]    = Vector{Int64}()
        Gᵢ[b]    = Vector{Int64}()
        out_L[b] = Vector{Int64}()
        in_L[b]  = Vector{Int64}()
        cb[b] = 5
    end

    for i in keys(network_data["load"])
        d = network_data["load"][i]["index"]
        b = network_data["load"][i]["load_bus"]
        w[d] = wsample([50, 100, 150, 200, 250, 300, 500, 600, 700, 1000], [8, 8, 10, 8, 2, 3, 1, 1, 1, .5], 1)[1];                                ## priority level of load d

        push!(Dᵢ[b], d)
        push!(D, d)
        demand[d] = network_data["load"][i]["pd"];
    end

    for i in keys(network_data["gen"])
        g = network_data["gen"][i]["index"]
        b = network_data["gen"][i]["gen_bus"]

        push!(G, g)
        push!(Gᵢ[b], g)

        smax[g] = network_data["gen"][i]["pmax"]
        smin[g] = network_data["gen"][i]["pmin"]
        M[g] = network_data["gen"][i]["ramp_30"]
        cg[g] = wsample([50, 1000, 2500], [0.2, .75, 0.05], 1)[1]  
        slope[g] = network_data["gen"][i]["cost"][2]
        intercept[g] = network_data["gen"][i]["cost"][1]
        C_start[g] = network_data["gen"][i]["startup"]
        C_down[g] = network_data["gen"][i]["shutdown"]

    end


    for i in keys(network_data["branch"])
        l = (network_data["branch"][i]["f_bus"], network_data["branch"][i]["t_bus"])

        if l ∉ L 
            push!(L, l) 
            push!(out_L[l[1]], l[2])
            push!(in_L[l[2]], l[1])

            _b[l] = - 1/network_data["branch"][i]["br_x"]                          ## total line charging susceptance
            W[l] = network_data["branch"][i]["rate_a"]         
            cl[l] = 0.285 *  branchInfo[parse(Int64,i), :Length] 
        else
            _b[l] = _b[l] - 1/network_data["branch"][i]["br_x"] 
            W[l] = W[l] + network_data["branch"][i]["rate_a"] 
            cl[l] = cl[l] + 0.285 *  branchInfo[parse(Int64,i), :Length] 
        end

    end

    paramOPF = ParamOPF(_b, θmax, θmin, W, smax, smin, M, slope, intercept, C_start, C_down)
    indexSets = IndexSets(D, G, L, B ,T, Dᵢ, Gᵢ, out_L, in_L)
    paramDemand = ParamDemand(demand, w, cb, cg, cl, 1e4)
 
     return (indexSets = indexSets, 
             paramOPF = paramOPF, 
             paramDemand = paramDemand)
end

"""
scenario_tree_generation(T::Int64 = 10, numRealization::Int64 = 3, indexSets::IndexSets = indexSets)

Simulation

# Arguments

  1. `T`: Time horizon.
  2. `numRealization`: Each stage has `numRealization` realizations.
# Returns
  1. `scenarioTree`: A scenario tree with `T` stages and `numRealization` realizations at each stage.
"""
function scenario_tree_generation(; T::Int64 = 10, numRealization::Int64 = 3, indexSets::IndexSets = indexSets)
    tree = Dict{Int64, StageRealization}()
    nodes = Dict{Int64, RandomVariables}()
    for n in 1:1
        deviation = Dict{Int64, Float64}()
        for d in indexSets.D
            deviation[d] = 1.
        end
        nodes[n] = RandomVariables(deviation)
    end
    stageRealization = StageRealization(nodes, Dict{Int64, Float64}(n => 1/numRealization for n in 1:numRealization))
    tree[1] = stageRealization
    for t in 2:T
        nodes = Dict{Int64, RandomVariables}()
        for n in 1:numRealization
            deviation = Dict{Int64, Float64}()
            for d in indexSets.D
                deviation[d] = rand(Uniform(0.9, 1.2))
            end
            nodes[n] = RandomVariables(deviation)
        end
        stageRealization = StageRealization(nodes, Dict{Int64, Float64}(n => 1/numRealization for n in 1:numRealization))
        tree[t] = stageRealization
    end
    scenarioTree = ScenarioTree(tree)
    return scenarioTree
end

"""
sample_scenarios(numRealization::Int64 = 3, scenarioTree::ScenarioTree = scenarioTree)

Simulation

# Arguments

  1. `numScenarios`: The number of scenarios will be sampled

# Returns
  1. `Ξ`: A subset of scenarios.
"""
function sample_scenarios(; numScenarios::Int64 = 10, scenarioTree::ScenarioTree = scenarioTree)
    Ξ = Dict{Int64, Dict{Int64, RandomVariables}}()
    for ω in 1:numScenarios
        ξ = Dict{Int64, RandomVariables}()
        ξ[1] = scenarioTree.tree[1].nodes[1]
        n = wsample(collect(keys(scenarioTree.tree[1].prob)), collect(values(scenarioTree.tree[1].prob)), 1)[1]
        for t in 2:length(keys(scenarioTree.tree))
            ξ[t] = scenarioTree.tree[t].nodes[n]
            n = wsample(collect(keys(scenarioTree.tree[t].prob)), collect(values(scenarioTree.tree[t].prob)), 1)[1]
        end
        Ξ[ω] = ξ
    end
    return Ξ
end
