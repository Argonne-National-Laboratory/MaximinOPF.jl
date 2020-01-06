# "Case9 feas problem with no lines cut should have value 0.0"
# "Case9 feas problem with lines 1,4,7 cut should have value 3.451901668361513"
# "Case9 feas problem with lines 1--9 cut should have value 4.6"
# "Case30 feas problem with no lines cut should have value 0.0"
# "Case30 feas problem with lines 8,9,10,40 cut should have value 0.937"
# "Case30 feas problem with all lines cut should have value 2.57327375"
# "Case57 feas problem with no lines cut should have value 0.0"
# "Case57 feas problem with lines 41,80 cut should have value 0.7597369902683009"
# "Case57 feas problem with all lines cut should have value 6.2313456"

casesConic = [
	Dict("file" => "../data/case9.m", 
		"name" => "case9", 
		"attack_budget" => 3,
		"inactive_indices" => [],
		"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>3.769156, "SOC_Maxmin"=>3.4519, "SDP_FP"=>0.0, "SDP_Minmax"=>3.769156, "SDP_Maxmin"=>3.4519 )
	),
	Dict("file" => "../data/case9.m", 
		"name" => "case9", 
		"attack_budget" => 0,
		"inactive_indices" => [1,4,7],
		"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>3.4519, "SOC_Minmax"=>3.4519, "SOC_Maxmin"=>3.4519, "SDP_FP"=>3.4519, "SDP_Minmax"=>3.4519, "SDP_Maxmin"=>3.4519 )
	),
	Dict("file" => "../data/case9.m", 
	 	"name" => "case9", 
	 	"attack_budget" => 3,
	 	"inactive_indices" => [],
	 	"protected_indices" => [2,3,5,6,8,9],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>3.4519, "SOC_Maxmin"=>3.4519, "SDP_FP"=>0.0, "SDP_Minmax"=>3.4519, "SDP_Maxmin"=>3.4519 )
	 ),
	 Dict("file" => "../data/case30.m", 
	 	"name" => "case30", 
	 	"attack_budget" => 4,
	 	"inactive_indices" => [],
	 	"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>1.971, "SOC_Maxmin"=>0.937, "SDP_FP"=>0.0, "SDP_Minmax"=>2.00117, "SDP_Maxmin"=>0.937 )
	),
	Dict("file" => "../data/case30.m", 
	 	"name" => "case30", 
	 	"attack_budget" => 0,
	 	"inactive_indices" => [8,9,10,40],
	 	"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.937, "SOC_Minmax"=>0.937, "SOC_Maxmin"=>0.937, "SDP_FP"=>0.937, "SDP_Minmax"=>0.937, "SDP_Maxmin"=>0.937 )
	),
	Dict("file" => "../data/case30.m", 
	 	"name" => "case30", 
	 	"attack_budget" => 4,
	 	"inactive_indices" => [],
	 	"protected_indices" => filter(l->!(l in [8,9,10,40]), 1:41), 
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>0.937, "SOC_Maxmin"=>0.937, "SDP_FP"=>0.0, "SDP_Minmax"=>0.937, "SDP_Maxmin"=>0.937 )
	),
	Dict("file" => "../data/case57.m", 
	 	"name" => "case57", 
	 	"attack_budget" => 2,
	 	"inactive_indices" => [],
	 	"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>2.768, "SOC_Maxmin"=>0.759737, "SDP_FP"=>0.0, "SDP_Minmax"=>3.06229, "SDP_Maxmin"=>0.77726 )
	),
	Dict("file" => "../data/case57.m", 
		"name" => "case57", 
		"attack_budget" => 0,
		"inactive_indices" => [41, 80],
		"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.759737, "SOC_Minmax"=>0.759737, "SOC_Maxmin"=>0.759737, "SDP_FP"=>0.0, "SDP_Minmax"=>0.77726, "SDP_Maxmin"=>0.77726 )
	),
	Dict("file" => "../data/case57.m", 
	 	"name" => "case57", 
	 	"attack_budget" => 4,
	 	"inactive_indices" => [],
	 	"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>4.511, "SOC_Maxmin"=>1.5793, "SDP_FP"=>0.0, "SDP_Minmax"=>4.619, "SDP_Maxmin"=>1.6 )
	),
	Dict("file" => "../data/case118.m", 
	 	"name" => "case118", 
	 	"expectedvalue" => 2.645, 
	 	"attack_budget" => 2,
	 	"inactive_indices" => [],
	 	"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>2.645, "SOC_Maxmin"=>1.447, "SDP_FP"=>0.0, "SDP_Minmax"=>2.70888, "SDP_Maxmin"=>1.448 )
	),
	Dict("file" => "../data/case118.m", 
	 	"name" => "case118", 
	 	"expectedvalue" => 4.051, 
	 	"attack_budget" => 4,
	 	"inactive_indices" => [],
	 	"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>4.051, "SOC_Maxmin"=>2.567, "SDP_FP"=>0.0, "SDP_Minmax"=>4.296, "SDP_Maxmin"=>2.568 )
	),
]

#=
function getTestcasesConic()
	pm_datas = []

	for i in 1:length(casesConic)
		case = casesConic[i]["file"]
		pm_data = PowerModels.parse_file(case)
		pm_data["attacker_budget"]=nLineAttacked ###Adding another key and entry
		pm_data["inactive_branches"]=lineindexs ###Adding another key and entry
		pm_data["protected_branches"]=protected_indices ###Adding another key and entry
		casename = casesConic[i]["name"]
		expected_values = casesConic[i]["expected_values"]
		protected_branches = casesConic[i]["protected_indices"]
		inactive_branches = casesConic[i]["inactive_indices"]
		K = casesConic[i]["attack_budget"]
	    push!(pm_datas, 
	    	Dict("pm_data" => pm_data, 
	    		"K" => K, 
	    		"name" => casename,
	    		"protected_indices" => protected_branches,
	    		"inactive_indices" => inactive_branches,
	    		"expectedvalueFP" => expectFP,
	    		"expectedvalueMinmax" => expectMinmax,
	    		"expectedvalueMaxmin" => expectMaxmin
	    	))
	    push!(pm_datas,pm_data)
	end
	return pm_datas
end
=#
