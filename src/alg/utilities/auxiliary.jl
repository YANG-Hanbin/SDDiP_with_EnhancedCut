""" 
    function setup_core_point(
        stateInfo::StateInfo;
        indexSets::IndexSets = indexSets,
        paramOPF::ParamOPF = paramOPF, 
        param_cut::NamedTuple = param_cut   
    )::StateInfo
"""
function setup_core_point(
    stateInfo::StateInfo;
    indexSets::IndexSets = indexSets,
    paramOPF::ParamOPF = paramOPF, 
    param::NamedTuple = param,
    param_cut::NamedTuple = param_cut   
)::StateInfo
    if param_cut.core_point_strategy == "Mid"
        BinVar = Dict{Any, Dict{Any, Any}}(:y => Dict{Any, Any}(
            g => 0.5 for g in indexSets.G)
        );
        ContVar = Dict{Any, Dict{Any, Any}}(:s => Dict{Any, Any}(
            g => (paramOPF.smin[g] + paramOPF.smax[g])/2 for g in indexSets.G)
        );
        if stateInfo.ContAugState == nothing
            ContAugState = nothing
        else
            ContAugState = Dict{Any, Dict{Any, Dict{Any, Any}}}(
                :s => Dict{Any, Dict{Any, Any}}(
                    g => Dict{Any, Any}(
                        k => .5 for k in keys(stateInfo.ContAugState[:s][g])
                    ) for g in indexSets.G
                )
            );
        end

        if stateInfo.ContStateBin == nothing
            ContStateBin = nothing
        else
            ContStateBin = Dict{Any, Dict{Any, Dict{Any, Any}}}(
                :s => Dict{Any, Dict{Any, Any}}(
                    g => Dict{Any, Any}(
                        i => .5 for i in 1:param.κ[g]
                    ) for g in indexSets.G
                )
            );
        end
    elseif param_cut.core_point_strategy == "Eps"
        BinVar = Dict{Any, Dict{Any, Any}}(:y => Dict{Any, Any}(
            g => stateInfo.BinVar[:y][g] * param_cut.ℓ + (1 - param_cut.ℓ)/2 for g in indexSets.G)
        );
        ContVar = Dict{Any, Dict{Any, Any}}(:s => Dict{Any, Any}(
            g => stateInfo.ContVar[:s][g] * param_cut.ℓ + (1 - param_cut.ℓ) * (paramOPF.smin[g] + paramOPF.smax[g])/2 for g in indexSets.G)
        );
        if stateInfo.ContAugState == nothing
            ContAugState = nothing
        else
            ContAugState = Dict{Any, Dict{Any, Dict{Any, Any}}}(
                :s => Dict{Any, Dict{Any, Any}}(
                    g => Dict{Any, Any}(
                        k => stateInfo.ContAugState[:s][g][k] * param_cut.ℓ + (1 - param_cut.ℓ)/2 for k in keys(stateInfo.ContAugState[:s][g])
                    ) for g in indexSets.G
                )
            );
        end

        if stateInfo.ContStateBin == nothing
            ContStateBin = nothing
        else
            ContStateBin = Dict{Any, Dict{Any, Dict{Any, Any}}}(
                :s => Dict{Any, Dict{Any, Any}}(
                    g => Dict{Any, Any}(
                        i => stateInfo.ContStateBin[:s][g][i] * param_cut.ℓ + (1 - param_cut.ℓ)/2 for i in 1:param.κ[g]
                    ) for g in indexSets.G
                )
            );
        end
    end

    return StateInfo(
        BinVar, 
        nothing, 
        ContVar, 
        nothing, 
        nothing, 
        nothing, 
        nothing, 
        nothing, 
        ContAugState,
        nothing,
        ContStateBin
    );
end

""" 
    function binarize_continuous_variable(
        state::Float64, 
        smax::Float64, 
        param::NamedTuple  
    )::StateInfo
"""
function binarize_continuous_variable(
    state::Float64, 
    smax::Float64, 
    param::NamedTuple  
)::Vector
    num_binary_vars = floor(Int, log2(smax / param.ε)) + 1  # the number of binary variables needed
    max_integer = round(Int, state / param.ε)  # 将连续变量映射到整数
    binary_representation = digits(max_integer, base=2, pad=num_binary_vars)  # 获取二进制表示（低位在前）
    return binary_representation
end