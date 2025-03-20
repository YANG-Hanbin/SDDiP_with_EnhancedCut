using Pkg;
Pkg.activate(".");
using JuMP, Gurobi, ParallelDataTransfer;
using Distributions, Statistics, StatsBase, Distributed, Random;
using Test, Dates, Printf;
using CSV, DataFrames;
using JLD2, FileIO;
using PrettyTables;


const GRB_ENV = Gurobi.Env();

project_root = @__DIR__;

include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "def.jl"));
include(joinpath(project_root, "src", "GenerationExpansion", "SDDP", "setting.jl"));

T = 10; # 10, 15
num = 5; # 5, 10
algorithm = "SDDLP";
cutSelection = "ELC";
tightness = true;

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------- TO GENERATE TABLES ---------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
# 初始化DataFrame
algorithm = "SDDiP";
tightness = false;
result_df = DataFrame(cut=String[], T=Int[], num=Int[], best_LB=Float64[], 
                      final_gap=Float64[], total_iter=Int[], avg_iter_time=String[], 
                      best_LB_time=Float64[], best_LB_iter=Int[])

for cut in ["LC", "ELC", "ShrinkageLC"]
    for T in [10, 15]
        for num in [5, 10]
            try
                # 加载数据
                file_path = "/Users/aaron/SDDiP_with_EnhancedCut/src/GenerationExpansion/logger/Periods$T-Real$num/$algorithm-$cut-$tightness.jld2"
                solHistory = load(file_path)["sddpResults"][:solHistory]

                # 计算所需的统计数据
                best_LB, best_LB_idx = findmax(solHistory.LB)  # 最优LB及其索引
                final_gap = parse(Float64, replace(solHistory.gap[end], "%" => ""))  # 最终gap
                total_iter = solHistory.iter[end]  # 总迭代数
                iter_times = diff(solHistory.Time)  # 计算每次迭代的时间间隔
                avg_time = mean(iter_times)  # 计算平均迭代时间
                std_time = std(iter_times)   # 计算标准差
                avg_iter_time = @sprintf "%.2f (%.2f)" avg_time std_time  # 格式化字符串
                best_LB_time = solHistory.Time[best_LB_idx]  # 到达best LB的时间
                best_LB_iter = solHistory.iter[best_LB_idx]  # 到达best LB的迭代数

                # 添加到DataFrame
                push!(result_df, (cut, T, num, round(best_LB, digits = 1), round(final_gap, digits = 1), total_iter, avg_iter_time, round(best_LB_time, digits = 0), best_LB_iter))
            catch e
                @warn "Error processing file: $file_path" exception=(e, catch_backtrace())
            end
        end
    end
end

# 定义格式化函数，保留一位小数
column_formatter = function(x, i, j)
    if x isa Float64
        return @sprintf("%.1f", x)  # 保留一位小数
    elseif x isa Tuple  # 处理 iter_range 之类的元组数据
        return "$(x[1])--$(x[2])"
    else
        return string(x)  # 其他数据类型转换为字符串
    end
end

# 生成 LaTeX 表格
latex_table = pretty_table(
    String, 
    result_df, 
    backend=Val(:latex),
    formatters=(column_formatter,)
)

# 输出 LaTeX 代码
println(latex_table)
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------- TO GENERATE TABLES ---------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##