using Pkg
Pkg.activate(".")
using JuMP, Gurobi, PowerModels
using Statistics, StatsBase, Random, Dates, Distributions
using Distributed, ParallelDataTransfer
using CSV, DataFrames, Printf
using JLD2, FileIO
using StatsPlots, PlotThemes
using PrettyTables;

project_root = @__DIR__;
include(joinpath(project_root, "src", "alg", "utilities", "structs.jl"))
theme(:default)

case = "case30pwl"; 
tightness = true; 
cutSelection = :SMC; 
num = 10; T = 12; 
algorithm = :SDDPL;

for cut in [:LC, :PLC, :SMC]
    for T in [6, 8, 12] 
        for num in [3, 5, 10] 
            solHistory = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-$case/Periods$T-Real$num/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]
        end
    end
end


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------- TO GENERATE TABLES ---------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
# 初始化DataFrame
result_df = DataFrame(cut=Symbol[], T=Int[], num=Int[], best_LB=Float64[], 
                      final_gap=Float64[], total_iter=Int[], avg_iter_time=String[], 
                      best_LB_time=Int64[], best_LB_iter=Int[])

for cut in [:LC, :PLC, :SMC]
    for T in [6, 8, 12]
        for num in [3, 5, 10]
            try
                # 加载数据
                
                file_path = "/Users/aaron/SDDiP_with_EnhancedCut/src/alg/new_logger/numericalResults-case30pwl/Periods$T-Real$num/SDDP-$cut.jld2"
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




## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## ----------------------------------------------------------------------- the same cut with different instances -------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ## 
cutSelection = :SMC
sddlpResult63 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods6-Real3/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]
sddlpResult65 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods6-Real5/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]
sddlpResult610 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods6-Real10/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]

sddlpResult83 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods8-Real3/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]
sddlpResult85 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods8-Real5/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]
sddlpResult810 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods8-Real10/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]

sddlpResult123 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods12-Real3/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]
sddlpResult125 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods12-Real5/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]
sddlpResult1210 = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods12-Real10/$algorithm-$cutSelection.jld2")["sddpResults"][:solHistory]

sddlpResult63.new_gap = parse.(Float64, replace.(sddlpResult63.gap, "%" => ""))
sddlpResult65.new_gap = parse.(Float64, replace.(sddlpResult65.gap, "%" => ""))
sddlpResult610.new_gap = parse.(Float64, replace.(sddlpResult610.gap, "%" => ""))
sddlpResult83.new_gap = parse.(Float64, replace.(sddlpResult83.gap, "%" => ""))
sddlpResult85.new_gap = parse.(Float64, replace.(sddlpResult85.gap, "%" => ""))
sddlpResult810.new_gap = parse.(Float64, replace.(sddlpResult810.gap, "%" => ""))
sddlpResult123.new_gap = parse.(Float64, replace.(sddlpResult123.gap, "%" => ""))
sddlpResult125.new_gap = parse.(Float64, replace.(sddlpResult125.gap, "%" => ""))
sddlpResult1210.new_gap = parse.(Float64, replace.(sddlpResult1210.gap, "%" => ""))




timegap = @df sddlpResult63 plot(:Time, :new_gap, label="T = 6, #R = 3", 
                                                # title = "Gap vs. Time", 
                                                xlab = "Time", 
                                                ylab = "Gap (%)", 
                                                xlim = [0,500],
                                                titlefont = font(15,"Times New Roman"), 
                                                xguidefont=font(15,"Times New Roman"), 
                                                yguidefont=font(15,"Times New Roman"), 
                                                xtickfontsize=13, 
                                                ytickfontsize=13, 
                                                legendfontsize=11, 
                                                marker=(:xcross, 2, 1.), 
                                                color=:goldenrod,  # 使用蓝色
                                                yformatter=y->string(Int(y)),
                                                tickfont=font("Computer Modern"),
                                                legendfont=font("Times New Roman"), legend=:outerright,
                                                linestyle=:solid)  # 实线

@df sddlpResult65 plot!(:Time, :new_gap, marker=(:star, 2, 1.), label="T = 6, #R = 5", linestyle=:solid, color=:orange)
@df sddlpResult610 plot!(:Time, :new_gap, marker=(:hexagon, 2, 1.), label="T = 6, #R = 10", linestyle=:solid, color=:green)
@df sddlpResult83 plot!(:Time, :new_gap, marker=(:plus, 2, 1.), label="T = 8, #R = 3", linestyle=:solid, color=:Crimson)
@df sddlpResult85 plot!(:Time, :new_gap, marker=(:star, 2, 1.), label="T = 8, #R = 5", linestyle=:solid, color=:HotPink)
@df sddlpResult810 plot!(:Time, :new_gap, marker=(:hexagon, 2, 1.), label="T = 8, #R = 10", linestyle=:solid, color=:brown)
@df sddlpResult123 plot!(:Time, :new_gap, marker=(:xcross, 2, 1.), label="T = 12, #R = 3", linestyle=:solid, color=:lightslategray)
@df sddlpResult125 plot!(:Time, :new_gap, marker=(:square, 2, 1.), label="T = 12, #R = 5", linestyle=:solid, color=:cyan)
@df sddlpResult1210 plot!(:Time, :new_gap, marker=(:diamond, 2, 1.), label="T = 12, #R = 10", linestyle=:solid, color=:palevioletred)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## ----------------------------------------------------------------------- the same instance with different cuts -------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ## 
T = 6; 
num = 10;
for T in [6, 8, 12]
    for num in [3, 5, 10]
        sddlpResultLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods$T-Real$num/$algorithm-LC.jld2")["sddpResults"][:solHistory]
        sddlpResultPLC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods$T-Real$num/$algorithm-PLC.jld2")["sddpResults"][:solHistory]
        sddlpResultSMC = load("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods$T-Real$num/$algorithm-SMC.jld2")["sddpResults"][:solHistory]

        sddlpResultLC.gap_float = parse.(Float64, replace.(sddlpResultLC.gap, "%" => "")) 
        sddlpResultPLC.gap_float = parse.(Float64, replace.(sddlpResultPLC.gap, "%" => ""))
        sddlpResultSMC.gap_float = parse.(Float64, replace.(sddlpResultSMC.gap, "%" => ""))

        lbtime = @df sddlpResultLC plot(
            :Time, 
            :LB, 
            label="LC",
            xlab = "Time (sec.)",
            xlim = [0, 3600],
            # ylim = [0, 200000],
            ylab = "Lower bounds (× 10³)",
            titlefont = font(15,"Times New Roman"),
            xguidefont=font(15,"Times New Roman"), 
            yguidefont=font(15,"Times New Roman"), 
            xtickfontsize=13, 
            ytickfontsize=13, 
            # size=(400,300),
            yformatter = y -> y/10^3,
            tickfont=font("Computer Modern"),
            marker=(:plus, 2, 1.), 
            linewidth=1,
            linestyle=:dot, 
            legend=:outertop,  # legend 在顶部
            legendfontsize=11, 
            legendfont=font("Times New Roman"), 
            legend_column=3,  # legend 列数减少，使其松散
            legend_spacing=6,  # 控制 legend 之间的间距
            color=:DodgerBlue
        )  
        # @df sddlpResultLC plot!(:Time, :UB, label="", linewidth=1, marker=(:circle, 2, 1.), color=:DodgerBlue)  
        @df sddlpResultPLC plot!(:Time, :LB, label="PLC", linewidth=1, marker=(:vline, 2, 1.), linestyle=:dash, color=:red) 
        # @df sddlpResultPLC plot!(:Time, :UB, label="", linewidth=1, marker=(:star, 2, 1.), color=:red)  
        @df sddlpResultSMC plot!(:Time, :LB, label="SMC", linewidth=1, marker=(:cross, 2, 1.), linestyle=:dashdot, color=:black)
        # @df sddlpResultSMC plot!(:Time, :UB, linewidth=1, marker=(:hexagon, 2, 1.), label="", color=:orange)
        lbtime |> save("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods$T-Real$num/lower_bound_Time_Period$T-Real$num.pdf")



        lbiter = @df sddlpResultLC plot(
            :Iter, 
            :LB, 
            label="LC",
            xlab = "Iteration",
            xlim = [0, 40],
            # ylim = [0, 200000],
            ylab = "Lower bounds (× 10³)",
            titlefont = font(15,"Times New Roman"),
            xguidefont=font(15,"Times New Roman"), 
            yguidefont=font(15,"Times New Roman"), 
            xtickfontsize=13, 
            ytickfontsize=13, 
            # size=(400,300),
            yformatter = y -> y/10^3,
            tickfont=font("Computer Modern"),
            marker=(:plus, 2, 1.), 
            linewidth=1,
            linestyle=:dot, 
            legend=:outertop,  # legend 在顶部
            legendfontsize=11, 
            legendfont=font("Times New Roman"), 
            legend_column=3,  # legend 列数减少，使其松散
            legend_spacing=6,  # 控制 legend 之间的间距
            color=:DodgerBlue
        )  
        # @df sddlpResultLC plot!(:Time, :UB, label="", linewidth=1, marker=(:circle, 2, 1.), color=:DodgerBlue)  
        @df sddlpResultPLC plot!(:Iter, :LB, label="PLC", linewidth=1, marker=(:vline, 2, 1.), linestyle=:dash, color=:red) 
        # @df sddlpResultPLC plot!(:Time, :UB, label="", linewidth=1, marker=(:star, 2, 1.), color=:red)  
        @df sddlpResultSMC plot!(:Iter, :LB, label="SMC", linewidth=1, marker=(:cross, 2, 1.), linestyle=:dashdot, color=:black)
        # @df sddlpResultSMC plot!(:Time, :UB, linewidth=1, marker=(:hexagon, 2, 1.), label="", color=:orange)
        lbiter |> save("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods$T-Real$num/lower_bound_Iter_Period$T-Real$num.pdf")

        gaptime = @df sddlpResultLC plot(
            :Time, 
            :gap_float, 
            label="LC",
            xlab = "Time (sec.)",
            xlim = [0, 3600],
            ylim = [0, 100],
            ylab = "Gap (%)",
            titlefont = font(15,"Times New Roman"),
            xguidefont=font(15,"Times New Roman"), 
            yguidefont=font(15,"Times New Roman"), 
            xtickfontsize=13, 
            ytickfontsize=13, 
            # size=(400,300),
            yformatter = y -> y,
            tickfont=font("Computer Modern"),
            marker=(:plus, 2, 1.), 
            linewidth=1,
            linestyle=:dot, 
            legend=:outertop,  # legend 在顶部
            legendfontsize=11, 
            legendfont=font("Times New Roman"), 
            legend_column=3,  # legend 列数减少，使其松散
            legend_spacing=6,  # 控制 legend 之间的间距
            color=:DodgerBlue
        )  
        # @df sddlpResultLC plot!(:Time, :UB, label="", linewidth=1, marker=(:circle, 2, 1.), color=:DodgerBlue)  
        @df sddlpResultPLC plot!(:Time, :gap_float, label="PLC", linewidth=1, marker=(:vline, 2, 1.), linestyle=:dash, color=:red) 
        # @df sddlpResultPLC plot!(:Time, :UB, label="", linewidth=1, marker=(:star, 2, 1.), color=:red)  
        @df sddlpResultSMC plot!(:Time, :gap_float, label="SMC", linewidth=1, marker=(:cross, 2, 1.), linestyle=:dashdot, color=:black)
        # @df sddlpResultSMC plot!(:Time, :UB, linewidth=1, marker=(:hexagon, 2, 1.), label="", color=:orange)
        gaptime |> save("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods$T-Real$num/gpa_Time_Period$T-Real$num.pdf")

        gapiter = @df sddlpResultLC plot(
            :Iter, 
            :gap_float, 
            label="LC",
            xlab = "Iteration",
            xlim = [0, 40],
            # ylim = [0, 200000],
            ylab = "Gap (%)",
            titlefont = font(15,"Times New Roman"),
            xguidefont=font(15,"Times New Roman"), 
            yguidefont=font(15,"Times New Roman"), 
            xtickfontsize=13, 
            ytickfontsize=13, 
            # size=(400,300),
            yformatter = y -> y,
            tickfont=font("Computer Modern"),
            marker=(:plus, 2, 1.), 
            linewidth=1,
            linestyle=:dot, 
            legend=:outertop,  # legend 在顶部
            legendfontsize=11, 
            legendfont=font("Times New Roman"), 
            legend_column=3,  # legend 列数减少，使其松散
            legend_spacing=6,  # 控制 legend 之间的间距
            color=:DodgerBlue
        )  
        # @df sddlpResultLC plot!(:Time, :UB, label="", linewidth=1, marker=(:circle, 2, 1.), color=:DodgerBlue)  
        @df sddlpResultPLC plot!(:Iter, :gap_float, label="PLC", linewidth=1, marker=(:vline, 2, 1.), linestyle=:dash, color=:red) 
        # @df sddlpResultPLC plot!(:Time, :UB, label="", linewidth=1, marker=(:star, 2, 1.), color=:red)  
        @df sddlpResultSMC plot!(:Iter, :gap_float, label="SMC", linewidth=1, marker=(:cross, 2, 1.), linestyle=:dashdot, color=:black)
        # @df sddlpResultSMC plot!(:Time, :UB, linewidth=1, marker=(:hexagon, 2, 1.), label="", color=:orange)
        gapiter |> save("/Users/aaron/SDDiP_with_EnhancedCut/src/alg/logger/numericalResults-case30pwl/Periods$T-Real$num/gap_Iter_Period$T-Real$num.pdf")
    end
end


using Plots
gr()  # 使用 GR 后端

# 主图（左 y 轴）：Lower Bound vs Time
lbtime = @df sddlpResultLC plot(
    :Time, 
    :LB, 
    label="LC",
    xlab="Time (sec.)",
    xlim=[0, 3600],
    ylab="Lower bounds (× 10³)",
    titlefont=font(15, "Times New Roman"),
    xguidefont=font(15, "Times New Roman"), 
    yguidefont=font(15, "Times New Roman"), 
    xtickfontsize=13, 
    ytickfontsize=13, 
    yformatter=y -> y / 10^3,  # 让 LB 以千为单位显示
    tickfont=font("Computer Modern"),
    marker=(:plus, 2, 1.), 
    linewidth=1,
    linestyle=:dot, 
    legend=:outertop,  # legend 在顶部
    legendfontsize=11, 
    legendfont=font("Times New Roman"), 
    legend_column=3,  # legend 列数减少，使其松散
    legend_spacing=6,  # 控制 legend 之间的间距
    color=:DodgerBlue
)  

@df sddlpResultPLC plot!(:Time, :LB, label="PLC", linewidth=1, marker=(:vline, 2, 1.), linestyle=:dash, color=:red) 
@df sddlpResultSMC plot!(:Time, :LB, label="SMC", linewidth=1, marker=(:cross, 2, 1.), linestyle=:dashdot, color=:black)

# 右 y 轴（Gap vs Time）
twinx()  # 创建右侧 y 轴
@df sddlpResultLC plot!(:Time, :gap_float, label="Gap (%)", linewidth=2, linestyle=:solid, marker=(:circle, 3, 1.), color=:green, yaxis=:right)