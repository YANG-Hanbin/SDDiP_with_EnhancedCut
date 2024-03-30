T = 3; num = 5; cutSelection = "LC";
result = load("src/GenerationExpansion/data/testData_stage($T)_real($num)/generInt_Interval_result_stage($T)_real($num)_$cutSelection.jld2")["result"]
result[:solHistory]
describe(result[:solHistory])
(result[:solHistory][1,3] - result[:solHistory][2,2])/result[:solHistory][1,3]

for cutSelection in ["LC", "ELC", "ShrinkageLC"]
    for T in [3, 5]
        for num in [5, 10]
            result = load("src/GenerationExpansion/data/testData_stage($T)_real($num)/binary_Interval_result_stage($T)_real($num)_$cutSelection.jld2")["result"]
            result[:solHistory][1:500,: ]
            describe(result[:solHistory][1:500,: ])
        end
    end
end
