using LightGraphs
using StatsBase
using SparseArrays



function Scenarios(nv::Int64 = 5, ne::Int64 = 10, S::Int64 = 10)
    network = SimpleGraph(nv, ne)
    collect(edges(network))

    N = Vector{Int64}(1:nv)

    A = Vector{Tuple{Int64, Int64}}()
    D = Vector{Tuple{Int64, Int64}}()
    Dᶜ = Vector{Tuple{Int64, Int64}}()
    c = Dict{Tuple{Int64, Int64}, Real}()
    q = Dict{Tuple{Int64, Int64}, Real}()
    r = Dict{Tuple{Int64, Int64}, Real}()
    for e in collect(edges(network))
        push!(A, (e.src, e.dst))
        q[(e.src, e.dst)] = rand()
        if rand() < 0.5
            push!(D , (e.src, e.dst))
            c[(e.src, e.dst)] = floor(rand() * 100)
            r[(e.src, e.dst)] = q[(e.src, e.dst)] * .5
        else 
            push!(Dᶜ, (e.src, e.dst))
        end
    end


    distmx = fill(Inf, nv, nv)
    for key in keys(q)
        distmx[key...] = q[key]
    end

    b = sum(values(c)) * .5

    Ω = Dict{Int64, ScenarioData}()


    for ω in 1:S
        sω =  sample(N, Weights([1:nv...]))
        tω = sample(N, Weights([1:nv...]))
        ϕω = Dict{Int64, Real}()
        for (i,j) in D 
            dist = dijkstra_shortest_paths(network, j, distmx).dists
            ϕω[j] = dist[j]
        end
        Ω[ω] = ScenarioData(sω, tω, ϕω)
    end

    stageData = StageData(c, b)
    indexSets = IndexSets(N, A, D, Dᶜ, Ω)
    return (stageData = stageData, indexSets = indexSets, Ω = Ω)
end






