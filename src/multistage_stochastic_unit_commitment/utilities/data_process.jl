using Pkg
Pkg.activate(".")
using JuMP, Gurobi, PowerModels
using Statistics, StatsBase, Random, Dates, Distributions
using Distributed, ParallelDataTransfer
using CSV, DataFrames, Printf
using JLD2, FileIO
using StatsPlots, PlotThemes
using PrettyTables;
using VegaLite, VegaDatasets;

project_root = @__DIR__;
include(joinpath(project_root, "src", "multistage_stochastic_unit_commitment", "utilities", "structs.jl"))
theme(:default)

case = "case30"; 
tightness = true; 
cutSelection = :SMC; 
num = 10; T = 12; 
algorithm = :SDDPL;

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------- TO GENERATE TABLES ---------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
# 初始化DataFrame
med_method = :IntervalMed # :ExactPoint, :IntervalMed
result_df = DataFrame(
    cut=Symbol[], 
    T=Int[], 
    num=Int[], 
    best_LB=Float64[],         
    final_gap=Float64[], 
    total_iter=Int[], 
    avg_iter_time=String[],         
    # best_LB_time=Float64[], 
    # best_LB_iter=Int[],
    gap_under_1_time=Union{Missing, Float64}[],
    gap_under_1_iter=Union{Missing, Int}[]
);

for cut in [:LC, :PLC, :SMC]
    for T in [6, 8, 12]
        for num in [5, 10]
            try
                if cut == :PLC 
                    file_path = "/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDPL-$cut-$med_method-0.0.jld2"
                else
                    file_path = "/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDPL-$cut-$med_method.jld2"
                end
                solHistory = load(file_path)["sddpResults"][:solHistory]

                # file_path = "/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDiP-$cut-256.jld2"
                # solHistory = load(file_path)["sddpResults"][:solHistory]

                # file_path = "/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDP-$cut.jld2"
                # solHistory = load(file_path)["sddpResults"][:solHistory]

                # 计算所需的统计数据
                best_LB, best_LB_idx = findmax(solHistory.LB)  # 最优LB及其索引
                final_gap = parse(Float64, replace(solHistory.gap[end], "%" => ""))  # 最终gap
                total_iter = solHistory.Iter[end]  # 总迭代数
                iter_times = diff(solHistory.Time)  # 计算每次迭代的时间间隔
                avg_time = mean(iter_times)  # 计算平均迭代时间
                std_time = std(iter_times)   # 计算标准差
                avg_iter_time = @sprintf "%.1f (%.1f)" avg_time std_time  # 格式化字符串
                best_LB_time = solHistory.Time[best_LB_idx]  # 到达best LB的时间
                best_LB_iter = solHistory.Iter[best_LB_idx]  # 到达best LB的迭代数

                # 将 gap 列（字符串）转换为 Float64 含义的百分数
                gap_vals = parse.(Float64, replace.(solHistory.gap, "%" => ""))

                # 找到 gap 第一次小于 1.0 的位置
                below1_idx = findfirst(<(1.0), gap_vals)

                # 初始化默认值
                gap_under_1_iter = missing
                gap_under_1_time = missing

                if below1_idx !== nothing
                    gap_under_1_iter = solHistory.Iter[below1_idx]
                    gap_under_1_time = solHistory.Time[below1_idx]
                end

                # 添加到DataFrame
                push!(result_df, (
                    cut, T, num, best_LB, final_gap, total_iter, 
                    avg_iter_time, 
                    # best_LB_time, best_LB_iter,
                    gap_under_1_time, gap_under_1_iter
                    )
                );
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
## ----------------------------------------------------------------------- the same instance with different cuts -------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ## 
med_method = :IntervalMed # :ExactPoint, :IntervalMed
algorithm = :SDDPL
for T in [6, 8, 12]
    for num in [5, 10]
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]  # 蓝色、橙色、绿色
        # 读取数据
        # sddlpResultLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-LC.jld2")["sddpResults"][:solHistory]
        # sddlpResultPLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-PLC.jld2")["sddpResults"][:solHistory]
        # sddlpResultSMC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-SMC.jld2")["sddpResults"][:solHistory]

        sddlpResultLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-LC-$med_method.jld2")["sddpResults"][:solHistory]
        sddlpResultPLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-PLC-$med_method-0.0.jld2")["sddpResults"][:solHistory]
        sddlpResultSMC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-SMC-$med_method.jld2")["sddpResults"][:solHistory]

        # 处理 gap 数据
        sddlpResultLC.gap_float = parse.(Float64, replace.(sddlpResultLC.gap, "%" => "")) 
        sddlpResultPLC.gap_float = parse.(Float64, replace.(sddlpResultPLC.gap, "%" => ""))
        sddlpResultSMC.gap_float = parse.(Float64, replace.(sddlpResultSMC.gap, "%" => ""))

        # 统一数据格式
        df_LC = DataFrame(Iter=sddlpResultLC.Iter, Time=sddlpResultLC.Time, LB=sddlpResultLC.LB ./ 10^3, Cut="LC")
        df_PLC = DataFrame(Iter=sddlpResultPLC.Iter, Time=sddlpResultPLC.Time, LB=sddlpResultPLC.LB ./ 10^3, Cut="PLC")
        df_SMC = DataFrame(Iter=sddlpResultSMC.Iter, Time=sddlpResultSMC.Time, LB=sddlpResultSMC.LB ./ 10^3, Cut="SMC")

        # 合并数据
        df = vcat(df_LC, df_PLC, df_SMC)

        df |> @vlplot(
            :line,
            x={:Time, axis={title="Time (sec.)", titleFontSize=25, labelFontSize=25,}},
            y={:LB, axis={title="Lower bounds (× 10³)", titleFontSize=25, labelFontSize=25}},
            color={
                :Cut, 
                legend={title=nothing, orient="top", columns=3}, 
                scale={domain=["LC", "PLC", "SMC"],  # 这里显式定义颜色顺序
                    range=colors}  # 绑定对应颜色
            },  
            strokeDash={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], 
                    range=[[5, 3], [10, 2], [10, 5, 2, 5]]}  # 虚线样式
            },  
            shape={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], 
                    range=["circle", "diamond", "cross"]}  # 形状
            },  
            strokeWidth={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], 
                    range=[1, 1, 1]}  # LC 粗，PLC 中等，SMC 细
            },  
            width=500,
            height=350,
            config={ 
                axis={
                    labelFont="Times New Roman", 
                    titleFont="Times New Roman"
                    }, 
                legend={
                    labelFont="Times New Roman", 
                    titleFont="Times New Roman",
                    labelFontSize=25,   # 调整 legend 标签字体大小
                    symbolSize=150,      # 增大 legend 符号大小
                    symbolStrokeWidth=3  # 增加 legend 线条粗细
                }, 
                title={font="Times New Roman"} 
            }
        ) |> save("$(homedir())/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/lower_bound_Time_Period$T-Real$num-$med_method.pdf")


        sddlpResultLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDPL-LC-$med_method.jld2")["sddpResults"][:solHistory]
        sddlpResultPLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDPL-PLC-$med_method-0.0.jld2")["sddpResults"][:solHistory]
        sddlpResultSMC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDPL-SMC-$med_method.jld2")["sddpResults"][:solHistory]

        # 处理 gap 数据
        sddlpResultLC.gap_float = parse.(Float64, replace.(sddlpResultLC.gap, "%" => "")) 
        sddlpResultPLC.gap_float = parse.(Float64, replace.(sddlpResultPLC.gap, "%" => ""))
        sddlpResultSMC.gap_float = parse.(Float64, replace.(sddlpResultSMC.gap, "%" => ""))

        # 统一数据格式
        df_LC = DataFrame(Iter=sddlpResultLC.Iter, Time=sddlpResultLC.Time, LB=sddlpResultLC.LB ./ 10^3, Cut="LC")
        df_PLC = DataFrame(Iter=sddlpResultPLC.Iter, Time=sddlpResultPLC.Time, LB=sddlpResultPLC.LB ./ 10^3, Cut="PLC")
        df_SMC = DataFrame(Iter=sddlpResultSMC.Iter, Time=sddlpResultSMC.Time, LB=sddlpResultSMC.LB ./ 10^3, Cut="SMC")

        df_LC = filter(row -> row.Iter < 35, df_LC)
        df_PLC = filter(row -> row.Iter < 35, df_PLC)
        df_SMC = filter(row -> row.Iter < 35, df_SMC)
        

        df = vcat(df_LC, df_PLC, df_SMC)

        df |> @vlplot(
            :line,
            x={:Iter, 
                axis={
                    title="Iteration", 
                    titleFontSize=25, 
                    labelFontSize=25, 
                    values=collect(0:5:35)  # 只显示 0, 5, 10, 15, ..., 35
                }
            },
            y={:LB, axis={title="Lower bounds (× 10³)", titleFontSize=25, labelFontSize=25}},
            color={
                :Cut, 
                legend={title=nothing, orient="top", columns=3}, 
                scale={domain=["LC", "PLC", "SMC"], range=colors}
            },  
            strokeDash={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], range=[[5, 3], [10, 2], [10, 5, 2, 5]]}
            },  
            shape={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], range=["circle", "diamond", "cross"]}
            },  
            strokeWidth={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], range=[1, 1, 1]}
            },  
            width=500,
            height=350,
            config={ 
                axis={
                    labelFont="Times New Roman", 
                    titleFont="Times New Roman"
                }, 
                legend={
                    labelFont="Times New Roman", 
                    titleFont="Times New Roman",
                    labelFontSize=25,   
                    symbolSize=150,     
                    symbolStrokeWidth=3  
                }, 
                title={font="Times New Roman"} 
            }
        ) |> save("$(homedir())/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/lower_bound_Iter_Period$T-Real$num-$med_method.pdf")
    end
end

for T in [6, 8, 12]
    for num in [5, 10]
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]  # 蓝色、橙色、绿色
        # 读取数据
        # sddlpResultLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-LC.jld2")["sddpResults"][:solHistory]
        # sddlpResultPLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-PLC.jld2")["sddpResults"][:solHistory]
        # sddlpResultSMC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-SMC.jld2")["sddpResults"][:solHistory]

        sddlpResultLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDiP-LC-256.jld2")["sddpResults"][:solHistory]
        sddlpResultPLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDiP-PLC-256.jld2")["sddpResults"][:solHistory]
        sddlpResultSMC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDiP-SMC-256.jld2")["sddpResults"][:solHistory]

        # 处理 gap 数据
        sddlpResultLC.gap_float = parse.(Float64, replace.(sddlpResultLC.gap, "%" => "")) 
        sddlpResultPLC.gap_float = parse.(Float64, replace.(sddlpResultPLC.gap, "%" => ""))
        sddlpResultSMC.gap_float = parse.(Float64, replace.(sddlpResultSMC.gap, "%" => ""))

        # 统一数据格式
        df_LC = DataFrame(Iter=sddlpResultLC.Iter, Time=sddlpResultLC.Time, LB=sddlpResultLC.LB ./ 10^3, Cut="LC")
        df_PLC = DataFrame(Iter=sddlpResultPLC.Iter, Time=sddlpResultPLC.Time, LB=sddlpResultPLC.LB ./ 10^3, Cut="PLC")
        df_SMC = DataFrame(Iter=sddlpResultSMC.Iter, Time=sddlpResultSMC.Time, LB=sddlpResultSMC.LB ./ 10^3, Cut="SMC")

        # 合并数据
        df = vcat(df_LC, df_PLC, df_SMC)

        df |> @vlplot(
            :line,
            x={:Time, axis={title="Time (sec.)", titleFontSize=25, labelFontSize=25,}},
            y={:LB, axis={title="Lower bounds (× 10³)", titleFontSize=25, labelFontSize=25}},
            color={
                :Cut, 
                legend={title=nothing, orient="top", columns=3}, 
                scale={domain=["LC", "PLC", "SMC"],  # 这里显式定义颜色顺序
                    range=colors}  # 绑定对应颜色
            },  
            strokeDash={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], 
                    range=[[5, 3], [10, 2], [10, 5, 2, 5]]}  # 虚线样式
            },  
            shape={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], 
                    range=["circle", "diamond", "cross"]}  # 形状
            },  
            strokeWidth={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], 
                    range=[1, 1, 1]}  # LC 粗，PLC 中等，SMC 细
            },  
            width=500,
            height=350,
            config={ 
                axis={
                    labelFont="Times New Roman", 
                    titleFont="Times New Roman"
                    }, 
                legend={
                    labelFont="Times New Roman", 
                    titleFont="Times New Roman",
                    labelFontSize=25,   # 调整 legend 标签字体大小
                    symbolSize=150,      # 增大 legend 符号大小
                    symbolStrokeWidth=3  # 增加 legend 线条粗细
                }, 
                title={font="Times New Roman"} 
            }
        ) |> save("$(homedir())/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/sddip-lower_bound_Time_Period$T-Real$num-$med_method.pdf")


        sddlpResultLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDiP-LC-256.jld2")["sddpResults"][:solHistory]
        sddlpResultPLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDiP-PLC-256.jld2")["sddpResults"][:solHistory]
        sddlpResultSMC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDiP-SMC-256.jld2")["sddpResults"][:solHistory]

        # 处理 gap 数据
        sddlpResultLC.gap_float = parse.(Float64, replace.(sddlpResultLC.gap, "%" => "")) 
        sddlpResultPLC.gap_float = parse.(Float64, replace.(sddlpResultPLC.gap, "%" => ""))
        sddlpResultSMC.gap_float = parse.(Float64, replace.(sddlpResultSMC.gap, "%" => ""))

        # 统一数据格式
        df_LC = DataFrame(Iter=sddlpResultLC.Iter, Time=sddlpResultLC.Time, LB=sddlpResultLC.LB ./ 10^3, Cut="LC")
        df_PLC = DataFrame(Iter=sddlpResultPLC.Iter, Time=sddlpResultPLC.Time, LB=sddlpResultPLC.LB ./ 10^3, Cut="PLC")
        df_SMC = DataFrame(Iter=sddlpResultSMC.Iter, Time=sddlpResultSMC.Time, LB=sddlpResultSMC.LB ./ 10^3, Cut="SMC")

        df_LC = filter(row -> row.Iter < 35, df_LC)
        df_PLC = filter(row -> row.Iter < 35, df_PLC)
        df_SMC = filter(row -> row.Iter < 35, df_SMC)
        

        df = vcat(df_LC, df_PLC, df_SMC)

        df |> @vlplot(
            :line,
            x={:Iter, 
                axis={
                    title="Iteration", 
                    titleFontSize=25, 
                    labelFontSize=25, 
                    values=collect(0:5:35)  # 只显示 0, 5, 10, 15, ..., 35
                }
            },
            y={:LB, axis={title="Lower bounds (× 10³)", titleFontSize=25, labelFontSize=25}},
            color={
                :Cut, 
                legend={title=nothing, orient="top", columns=3}, 
                scale={domain=["LC", "PLC", "SMC"], range=colors}
            },  
            strokeDash={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], range=[[5, 3], [10, 2], [10, 5, 2, 5]]}
            },  
            shape={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], range=["circle", "diamond", "cross"]}
            },  
            strokeWidth={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], range=[1, 1, 1]}
            },  
            width=500,
            height=350,
            config={ 
                axis={
                    labelFont="Times New Roman", 
                    titleFont="Times New Roman"
                }, 
                legend={
                    labelFont="Times New Roman", 
                    titleFont="Times New Roman",
                    labelFontSize=25,   
                    symbolSize=150,     
                    symbolStrokeWidth=3  
                }, 
                title={font="Times New Roman"} 
            }
        ) |> save("$(homedir())/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/sddip-lower_bound_Iter_Period$T-Real$num-$med_method.pdf")
    end
end

for T in [6, 8, 12]
    for num in [5, 10]
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]  # 蓝色、橙色、绿色
        # 读取数据
        # sddlpResultLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-LC.jld2")["sddpResults"][:solHistory]
        # sddlpResultPLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-PLC.jld2")["sddpResults"][:solHistory]
        # sddlpResultSMC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-SMC.jld2")["sddpResults"][:solHistory]

        sddlpResultLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDP-LC.jld2")["sddpResults"][:solHistory]
        sddlpResultPLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDP-PLC.jld2")["sddpResults"][:solHistory]
        sddlpResultSMC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDP-SMC.jld2")["sddpResults"][:solHistory]

        # 处理 gap 数据
        sddlpResultLC.gap_float = parse.(Float64, replace.(sddlpResultLC.gap, "%" => "")) 
        sddlpResultPLC.gap_float = parse.(Float64, replace.(sddlpResultPLC.gap, "%" => ""))
        sddlpResultSMC.gap_float = parse.(Float64, replace.(sddlpResultSMC.gap, "%" => ""))

        # 统一数据格式
        df_LC = DataFrame(Iter=sddlpResultLC.Iter, Time=sddlpResultLC.Time, LB=sddlpResultLC.LB ./ 10^3, Cut="LC")
        df_PLC = DataFrame(Iter=sddlpResultPLC.Iter, Time=sddlpResultPLC.Time, LB=sddlpResultPLC.LB ./ 10^3, Cut="PLC")
        df_SMC = DataFrame(Iter=sddlpResultSMC.Iter, Time=sddlpResultSMC.Time, LB=sddlpResultSMC.LB ./ 10^3, Cut="SMC")

        # 合并数据
        df = vcat(df_LC, df_PLC, df_SMC)

        df |> @vlplot(
            :line,
            x={:Time, axis={title="Time (sec.)", titleFontSize=25, labelFontSize=25,}},
            y={:LB, axis={title="Lower bounds (× 10³)", titleFontSize=25, labelFontSize=25}},
            color={
                :Cut, 
                legend={title=nothing, orient="top", columns=3}, 
                scale={domain=["LC", "PLC", "SMC"],  # 这里显式定义颜色顺序
                    range=colors}  # 绑定对应颜色
            },  
            strokeDash={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], 
                    range=[[5, 3], [10, 2], [10, 5, 2, 5]]}  # 虚线样式
            },  
            shape={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], 
                    range=["circle", "diamond", "cross"]}  # 形状
            },  
            strokeWidth={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], 
                    range=[1, 1, 1]}  # LC 粗，PLC 中等，SMC 细
            },  
            width=500,
            height=350,
            config={ 
                axis={
                    labelFont="Times New Roman", 
                    titleFont="Times New Roman"
                    }, 
                legend={
                    labelFont="Times New Roman", 
                    titleFont="Times New Roman",
                    labelFontSize=25,   # 调整 legend 标签字体大小
                    symbolSize=150,      # 增大 legend 符号大小
                    symbolStrokeWidth=3  # 增加 legend 线条粗细
                }, 
                title={font="Times New Roman"} 
            }
        ) |> save("$(homedir())/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/sddp-lower_bound_Time_Period$T-Real$num.pdf")


        sddlpResultLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDP-LC.jld2")["sddpResults"][:solHistory]
        sddlpResultPLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDP-PLC.jld2")["sddpResults"][:solHistory]
        sddlpResultSMC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDP-SMC.jld2")["sddpResults"][:solHistory]

        # 处理 gap 数据
        sddlpResultLC.gap_float = parse.(Float64, replace.(sddlpResultLC.gap, "%" => "")) 
        sddlpResultPLC.gap_float = parse.(Float64, replace.(sddlpResultPLC.gap, "%" => ""))
        sddlpResultSMC.gap_float = parse.(Float64, replace.(sddlpResultSMC.gap, "%" => ""))

        # 统一数据格式
        df_LC = DataFrame(Iter=sddlpResultLC.Iter, Time=sddlpResultLC.Time, LB=sddlpResultLC.LB ./ 10^3, Cut="LC")
        df_PLC = DataFrame(Iter=sddlpResultPLC.Iter, Time=sddlpResultPLC.Time, LB=sddlpResultPLC.LB ./ 10^3, Cut="PLC")
        df_SMC = DataFrame(Iter=sddlpResultSMC.Iter, Time=sddlpResultSMC.Time, LB=sddlpResultSMC.LB ./ 10^3, Cut="SMC")

        df_LC = filter(row -> row.Iter < 35, df_LC)
        df_PLC = filter(row -> row.Iter < 35, df_PLC)
        df_SMC = filter(row -> row.Iter < 35, df_SMC)
        

        df = vcat(df_LC, df_PLC, df_SMC)

        df |> @vlplot(
            :line,
            x={:Iter, 
                axis={
                    title="Iteration", 
                    titleFontSize=25, 
                    labelFontSize=25, 
                    values=collect(0:5:35)  # 只显示 0, 5, 10, 15, ..., 35
                }
            },
            y={:LB, axis={title="Lower bounds (× 10³)", titleFontSize=25, labelFontSize=25}},
            color={
                :Cut, 
                legend={title=nothing, orient="top", columns=3}, 
                scale={domain=["LC", "PLC", "SMC"], range=colors}
            },  
            strokeDash={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], range=[[5, 3], [10, 2], [10, 5, 2, 5]]}
            },  
            shape={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], range=["circle", "diamond", "cross"]}
            },  
            strokeWidth={
                :Cut, 
                scale={domain=["LC", "PLC", "SMC"], range=[1, 1, 1]}
            },  
            width=500,
            height=350,
            config={ 
                axis={
                    labelFont="Times New Roman", 
                    titleFont="Times New Roman"
                }, 
                legend={
                    labelFont="Times New Roman", 
                    titleFont="Times New Roman",
                    labelFontSize=25,   
                    symbolSize=150,     
                    symbolStrokeWidth=3  
                }, 
                title={font="Times New Roman"} 
            }
        ) |> save("$(homedir())/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/sddp-lower_bound_Iter_Period$T-Real$num.pdf")
    end
end



## ---------------------------------------------------------------------------------------------------------------------------------------- ##
## ------------------------------------------------------------  Bar Chart  --------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------------------------------------------------------------- ##
results = DataFrame(T = Int[], num = Int[], method = String[], avg_time = Float64[], std_time = Float64[], avg_LM_iter = Float64[], std_LM_iter = Float64[])

# 遍历不同的 (T, num) 组合
for T in [6, 8, 12]
    for num in [5, 10]
        for method in ["LC", "PLC", "SMC"]
            # 读取数据
            if method == "PLC" 
                data = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDPL-$method-$med_method-0.0.jld2")["sddpResults"][:solHistory]
            else
                data = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDPL-$method-$med_method.jld2")["sddpResults"][:solHistory]
            end
            

            df = DataFrame(Iter=data.Iter, time=data.time, LM_iters=data.LM_iter, Cut=method)
            
            # 计算平均值和标准差
            avg_time = mean(df.time)
            std_time = std(df.time)
            avg_LM_iter = mean(df.LM_iters)
            std_LM_iter = std(df.LM_iters)

            # 存入 DataFrame
            push!(results, (T, num, method, avg_time, std_time, avg_LM_iter, std_LM_iter))
        end
    end
end


# 颜色方案
# colors = ["#1E90FF", "#DC143C", "#006400"]  
colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]  # 蓝色、橙色、绿色
# 画第一个 Grouped Bar Chart（平均迭代时间）
results |>
@vlplot(
    :bar,
    x={"T:n", title="T", axis={labelFont="Times New Roman", labelFontSize=25, titleFontSize=25, labelAngle=0}},
    xOffset={"method:n", title="Cut"},
    y={"avg_time:q", title="Average Iteration Time", axis={labelFontSize=25, titleFontSize=25}},
    color={"method:n", scale={range=colors}, title=nothing, labelFontSize=25, titleFontSize=25},
    column={"num:n", title="R", 
            header={labelFont="Times New Roman", titleFont="Times New Roman", 
                    labelFontSize=25, titleFontSize=25}},  # 确保 R 也是 Times New Roman
    tooltip=[{ "T:n"}, {"num:n"}, {"method:n"}, {"avg_time:q"}, {"std_time:q"}],
    width=300, height=250,
    config={ 
        axis={labelFont="Times New Roman", titleFont="Times New Roman"}, 
        legend={
            labelFont="Times New Roman", titleFont="Times New Roman",
            labelFontSize=25, symbolSize=150, symbolStrokeWidth=3
        }, 
        title={font="Times New Roman"},
        bar={width=20}
    }
) + 
@vlplot(
    :errorbar,
    x={"T:n"},
    y={"avg_time:q", scale={zero=false}},
    yError={"std_time:q"}
) |> save("$(homedir())/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/$algorithm-AverageIterTime-$med_method.pdf")

# 画第二个 Grouped Bar Chart（平均 LM_iter 次数）
results |>
@vlplot(
    :bar,
    x={"T:n", title="T", axis={labelFont="Times New Roman", labelFontSize=25, titleFontSize=25, labelAngle=0}},
    xOffset={"method:n", title="Cut"},
    y={"avg_LM_iter:q", title="Average Iteration Counts", axis={labelFontSize=25, titleFontSize=25}},
    color={"method:n", scale={range=colors}, title=nothing, labelFontSize=25, titleFontSize=25},
    column={"num:n", title="R", 
            header={labelFont="Times New Roman", titleFont="Times New Roman", 
                    labelFontSize=25, titleFontSize=25}},  # 确保 R 也是 Times New Roman
    tooltip=[{ "T:n"}, {"num:n"}, {"method:n"}, {"avg_LM_iter:q"}, {"std_LM_iter:q"}],
    width=300, height=250,
    config={ 
        axis={labelFont="Times New Roman", titleFont="Times New Roman"}, 
        legend={
            labelFont="Times New Roman", titleFont="Times New Roman",
            labelFontSize=25, symbolSize=150, symbolStrokeWidth=3
        }, 
        title={font="Times New Roman"},
        bar={width=20}
    }
) + 
@vlplot(
    :errorbar,
    x={"T:n"},
    y={"avg_LM_iter:q", scale={zero=false}},
    yError={"std_LM_iter:q"}
) |> save("$(homedir())/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/$algorithm-LMiter-$med_method.pdf")





T = 12; num = 10;
results = DataFrame(alg = String[], method = String[], avg_time = Float64[], std_time = Float64[], avg_LM_iter = Float64[], std_LM_iter = Float64[])

# 遍历不同的 (T, num) 组合
for method in ["LC", "PLC", "SMC"]
    # 读取数据
    if method == "PLC" 
        data = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDPL-$method-$med_method-0.0.jld2")["sddpResults"][:solHistory]
    else
        data = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDPL-$method-$med_method.jld2")["sddpResults"][:solHistory]
    end
    

    df = DataFrame(Iter=data.Iter, time=data.time, LM_iters=data.LM_iter, Cut=method)
    
    # 计算平均值和标准差
    avg_time = mean(df.time)
    std_time = std(df.time)
    avg_LM_iter = mean(df.LM_iters)
    std_LM_iter = std(df.LM_iters)

    # 存入 DataFrame
    push!(results, ("SDDP-L", method, avg_time, std_time, avg_LM_iter, std_LM_iter))
end

for method in ["LC", "PLC", "SMC"]
    # 读取数据
    data = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDiP-$method-256.jld2")["sddpResults"][:solHistory]
    

    df = DataFrame(Iter=data.Iter, time=data.time, LM_iters=data.LM_iter, Cut=method)
    
    # 计算平均值和标准差
    avg_time = mean(df.time)
    std_time = std(df.time)
    avg_LM_iter = mean(df.LM_iters)
    std_LM_iter = std(df.LM_iters)

    # 存入 DataFrame
    push!(results, ("SDDiP", method, avg_time, std_time, avg_LM_iter, std_LM_iter))
end

for method in ["LC", "PLC", "SMC"]
    # 读取数据
    data = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/SDDP-$method.jld2")["sddpResults"][:solHistory]

    df = DataFrame(Iter=data.Iter, time=data.time, LM_iters=data.LM_iter, Cut=method)
    
    # 计算平均值和标准差
    avg_time = mean(df.time)
    std_time = std(df.time)
    avg_LM_iter = mean(df.LM_iters)
    std_LM_iter = std(df.LM_iters)

    # 存入 DataFrame
    push!(results, ("SDDP", method, avg_time, std_time, avg_LM_iter, std_LM_iter))
end

# 颜色方案
# colors = ["#1E90FF", "#DC143C", "#006400"]  
colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]  # 蓝色、橙色、绿色
# 画第一个 Grouped Bar Chart（平均迭代时间）
results |>
@vlplot(
    :bar,
    x={"alg:n", title = nothing, axis={labelFont="Times New Roman", labelFontSize=25, titleFontSize=25, labelAngle=0}},
    xOffset={"method:n", title="Cut"},
    y={"avg_time:q", title="Average Iteration Time", axis={labelFontSize=25, titleFontSize=25}},
    color={"method:n", scale={range=colors}, orient="top", title=nothing, labelFontSize=25, titleFontSize=25},
    tooltip=[{ "alg:n"},{"method:n"}, {"avg_time:q"}, {"std_time:q"}],
    width=300, height=250,
    config={ 
        axis={labelFont="Times New Roman", titleFont="Times New Roman"}, 
        legend={
            labelFont="Times New Roman", titleFont="Times New Roman",
            labelFontSize=25, symbolSize=150, orient="top", symbolStrokeWidth=3
        }, 
        title={font="Times New Roman"},
        bar={width=20}
    }
) |> save("$(homedir())/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/AverageIterTime-$med_method.pdf")

# 画第二个 Grouped Bar Chart（平均 LM_iter 次数）
results |>
@vlplot(
    :bar,
    x={"alg:n", title=nothing, axis={labelFont="Times New Roman", labelFontSize=25, titleFontSize=25, labelAngle=0}},
    xOffset={"method:n", title="Cut"},
    y={"avg_LM_iter:q", title="Average Iteration Counts", axis={labelFontSize=25, titleFontSize=25}},
    color={"method:n", scale={range=colors}, orient="top", title=nothing, labelFontSize=25, titleFontSize=25},
    tooltip=[{ "T:n"}, {"num:n"}, {"method:n"}, {"avg_LM_iter:q"}, {"std_LM_iter:q"}],
    width=300, height=250,
    config={ 
        axis={labelFont="Times New Roman", titleFont="Times New Roman"}, 
        legend={
            labelFont="Times New Roman", titleFont="Times New Roman",
            labelFontSize=25, symbolSize=150, orient="top", symbolStrokeWidth=3
        }, 
        title={font="Times New Roman"},
        bar={width=20}
    }
) |> save("$(homedir())/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-$case/Periods$T-Real$num/LMiter-$med_method.pdf")


## ---------------------------------------------------------------------------------------------------------------------------------------- ##
## --------------------------------------------------------  Deal with SDDiP  ------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------------------------------------------------------------- ##
# 初始化DataFrame
result_df = DataFrame(cut=Symbol[], T=Int[], num=Int[], best_LB=Float64[], 
                      final_gap=Float64[], total_iter=Int[], avg_iter_time=String[], 
                      best_LB_time=Int64[], best_LB_iter=Int[])

for cut in [:LC, :PLC, :SMC]
    for T in [6, 8, 12]
        for num in [5, 10]
            try
                # 加载数据
                
                file_path = "/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-case30pwl/Periods$T-Real$num/SDDiP-$cut-64.jld2"
                solHistory = load(file_path)["sddpResults"][:solHistory]

                # 计算所需的统计数据
                best_LB, best_LB_idx = findmax(solHistory.LB)  # 最优LB及其索引
                final_gap = parse(Float64, replace(solHistory.gap[end], "%" => ""))  # 最终gap
                total_iter = solHistory.Iter[end]  # 总迭代数
                iter_times = diff(solHistory.Time)  # 计算每次迭代的时间间隔
                avg_time = mean(iter_times)  # 计算平均迭代时间
                std_time = std(iter_times)   # 计算标准差
                avg_iter_time = @sprintf "%.1f (%.1f)" avg_time std_time  # 格式化字符串
                best_LB_time = solHistory.Time[best_LB_idx]  # 到达best LB的时间
                best_LB_iter = solHistory.Iter[best_LB_idx]  # 到达best LB的迭代数

                # 添加到DataFrame
                push!(result_df, (cut, T, num, round(best_LB, digits = 1), round(final_gap, digits = 1), total_iter, avg_iter_time, Int(round(best_LB_time)), best_LB_iter))
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


## ---------------------------------------------------------------------------------------------------------------------------------------- ##
## ------------------------------------------------------- Core Point Selection  ---------------------------------------------------------- ##
## ---------------------------------------------------------------------------------------------------------------------------------------- ##

colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]  

sddlpResult1 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-case30pwl/Periods6-Real5/SDDPL-PLC-IntervalMed-1.0.jld2")["sddpResults"][:solHistory]
sddlpResult2 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-case30pwl/Periods6-Real5/SDDPL-PLC-IntervalMed-0.2.jld2")["sddpResults"][:solHistory]
sddlpResult4 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-case30pwl/Periods6-Real5/SDDPL-PLC-IntervalMed-0.4.jld2")["sddpResults"][:solHistory]
sddlpResult6 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-case30pwl/Periods6-Real5/SDDPL-PLC-IntervalMed-0.6.jld2")["sddpResults"][:solHistory]
sddlpResult8 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-case30pwl/Periods6-Real5/SDDPL-PLC-IntervalMed-0.8.jld2")["sddpResults"][:solHistory]

# 处理 gap 数据
for res in [sddlpResult1, sddlpResult2, sddlpResult4, sddlpResult6, sddlpResult8]
    res.gap_float = parse.(Float64, replace.(res.gap, "%" => ""))
end

# 统一数据格式
df_1 = DataFrame(Iter=sddlpResult1.Iter, Time=sddlpResult1.Time, LB=sddlpResult1.LB ./ 10^3, Param="0.0")
df_2 = DataFrame(Iter=sddlpResult2.Iter, Time=sddlpResult2.Time, LB=sddlpResult2.LB ./ 10^3, Param="0.2")    
df_4 = DataFrame(Iter=sddlpResult4.Iter, Time=sddlpResult4.Time, LB=sddlpResult4.LB ./ 10^3, Param="0.4")
df_6 = DataFrame(Iter=sddlpResult6.Iter, Time=sddlpResult6.Time, LB=sddlpResult6.LB ./ 10^3, Param="0.6")
df_8 = DataFrame(Iter=sddlpResult8.Iter, Time=sddlpResult8.Time, LB=sddlpResult8.LB ./ 10^3, Param="0.8")


# 合并数据
df = vcat(df_1, df_2, df_4, df_6, df_8)

df |> @vlplot(
    :line,
    transform=[{filter="datum.Time <= 100"}],
    x={:Time, axis={title="Time (sec.)", titleFontSize=25, labelFontSize=25}},
    y={:LB, 
    axis={title="Lower bounds (× 10³)", titleFontSize=25, labelFontSize=25},
    scale={domain=[15, 130]}},
    color={
        :Param, 
        legend={title=nothing, orient="top", columns=5}, 
        scale={domain=["0.0", "0.2", "0.4", "0.6", "0.8"], range=colors}  # 绑定颜色到参数
    },  
    strokeDash={
        :Param, 
        scale={domain=["0.0", "0.2", "0.4", "0.6", "0.8"], 
            range=[[5, 3], [10, 2], [10, 5, 2, 5], [8, 4, 2, 4], [12, 6]]}  # 自定义虚线样式
    },  
    shape={
        :Param, 
        scale={domain=["0.0", "0.2", "0.4", "0.6", "0.8"], 
            range=["circle", "diamond", "cross", "triangle", "square"]}  # 形状
    },  
    strokeWidth={
        :Param, 
        scale={domain=["0.0", "0.2", "0.4", "0.6", "0.8"], 
            range=[1, 1, 1, 1, 1]}  # 线条粗细
    },  
    width=500,
    height=450,
    config={ 
        axis={
            labelFont="Times New Roman", 
            titleFont="Times New Roman"
        }, 
        legend={
            labelFont="Times New Roman", 
            titleFont="Times New Roman",
            labelFontSize=25,  
            symbolSize=150,      
            symbolStrokeWidth=3  
        }, 
        title={font="Times New Roman"} 
    }
) |> save("$(homedir())/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-case30pwl/Periods6-Real5/core_point_time.pdf")


df |> @vlplot(
    :line,
    transform=[{filter="datum.Iter <= 15"}],
    x={:Iter, axis={title="Iteration", titleFontSize=25, labelFontSize=25}},
    y={:LB, axis={title="Lower bounds (× 10³)", titleFontSize=25, labelFontSize=25},
    scale={domain=[15, 130]}},
    color={
        :Param, 
        legend={title=nothing, orient="top", columns=5}, 
        scale={domain=["0.0", "0.2", "0.4", "0.6", "0.8"], range=colors}  # 绑定颜色到参数
    },  
    strokeDash={
        :Param, 
        scale={domain=["0.0", "0.2", "0.4", "0.6", "0.8"], 
            range=[[5, 3], [10, 2], [10, 5, 2, 5], [8, 4, 2, 4], [12, 6]]}  # 自定义虚线样式
    },  
    shape={
        :Param, 
        scale={domain=["0.0", "0.2", "0.4", "0.6", "0.8"], 
            range=["circle", "diamond", "cross", "triangle", "square"]}  # 形状
    },  
    strokeWidth={
        :Param, 
        scale={domain=["0.0", "0.2", "0.4", "0.6", "0.8"], 
            range=[1, 1, 1, 1, 1]}  # 线条粗细
    },  
    width=500,
    height=450,
    config={ 
        axis={
            labelFont="Times New Roman", 
            titleFont="Times New Roman"
        }, 
        legend={
            labelFont="Times New Roman", 
            titleFont="Times New Roman",
            labelFontSize=25,  
            symbolSize=150,      
            symbolStrokeWidth=3  
        }, 
        title={font="Times New Roman"} 
    }
) |> save("$(homedir())/SDDiP_with_EnhancedCut/src/multistage_stochastic_unit_commitment/logger/numericalResults-case30pwl/Periods6-Real5/core_point_Iter.pdf")