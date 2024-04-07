T = 3; num = 10; cutSelection = "ELC"; tightness = "";
result = load("src/GenerationExpansion/data/testData_stage($T)_real($num)/binary_result_stage($T)_real($num)_$cutSelection.jld2")["result"]
result[:solHistory]
describe(result[:solHistory])
(result[:solHistory][1,3] - result[:solHistory][2,2])/result[:solHistory][1,3]

