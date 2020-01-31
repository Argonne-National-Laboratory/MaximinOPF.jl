# "Case9 feas problem with no lines cut should have value 0.0"
# "Case9 feas problem with lines 1,4,7 cut should have value 3.451901668361513"
# "Case9 feas problem with lines 1--9 cut should have value 4.6"
# "Case30 feas problem with no lines cut should have value 0.0"
# "Case30 feas problem with lines 8,9,10,40 cut should have value 0.937"
# "Case30 feas problem with all lines cut should have value 2.57327375"
# "Case57 feas problem with no lines cut should have value 0.0"
# "Case57 feas problem with lines 41,80 cut should have value 0.7597369902683009"
# "Case57 feas problem with all lines cut should have value 6.2313456"

# "test_cases[(FORM,INSTANCE,BUDGET,DESCRIPTION)]"
test_cases=Dict{Tuple{String,String,Int,Int},Dict{String,Any}}()

# "CASE 9"
# "SOC case 9"
test_cases[("SOC","case9",3,0)]=Dict{String,Any}("name"=>"case9K3SOC0", "attacker_budget"=>3)
test_cases[("SOC","case9",3,0)]["inactive_indices"]=[]
test_cases[("SOC","case9",3,0)]["protected_indices"]=[]
test_cases[("SOC","case9",3,0)]["expected_values"] = Dict{String,Float64}("FP"=>0.0, "Minmax"=>3.7544, "Maxmin"=>3.4519)

test_cases[("SOC","case9",3,1)]=Dict{String,Any}("name"=>"case9K3SOC1", "attacker_budget"=>0)
test_cases[("SOC","case9",3,1)]["inactive_indices"]=[1,4,7]
test_cases[("SOC","case9",3,1)]["protected_indices"]=[]
test_cases[("SOC","case9",3,1)]["expected_values"] = Dict{String,Float64}("FP"=>3.4519, "Minmax"=>3.4519, "Maxmin"=>3.4519)

test_cases[("SOC","case9",3,2)]=Dict{String,Any}("name"=>"case9K3SOC2", "attacker_budget"=>3)
test_cases[("SOC","case9",3,2)]["inactive_indices"]=[]
test_cases[("SOC","case9",3,2)]["protected_indices"]=[2,3,5,6,8,9]
test_cases[("SOC","case9",3,2)]["expected_values"] = Dict{String,Float64}("FP"=>0.0, "Minmax"=>3.4519, "Maxmin"=>3.4519)

# "SDP case 9"
test_cases[("SDP","case9",3,0)]=Dict{String,Any}("name"=>"case9K3SDP0", "attacker_budget"=>3)
test_cases[("SDP","case9",3,0)]["inactive_indices"]=[]
test_cases[("SDP","case9",3,0)]["protected_indices"]=[]
test_cases[("SDP","case9",3,0)]["expected_values"] = Dict{String,Float64}("FP"=>0.0, "Minmax"=>3.769156, "Maxmin"=>3.4519 )

test_cases[("SDP","case9",3,1)]=Dict{String,Any}("name"=>"case9K3SDP1", "attacker_budget"=>0)
test_cases[("SDP","case9",3,1)]["inactive_indices"]=[1,4,7]
test_cases[("SDP","case9",3,1)]["protected_indices"]=[]
test_cases[("SDP","case9",3,1)]["expected_values"] = Dict{String,Float64}("FP"=>3.4519, "Minmax"=>3.4519, "Maxmin"=>3.4519 )

test_cases[("SDP","case9",3,2)]=Dict{String,Any}("name"=>"case9K3SDP2", "attacker_budget"=>3)
test_cases[("SDP","case9",3,2)]["inactive_indices"]=[]
test_cases[("SDP","case9",3,2)]["protected_indices"]=[2,3,5,6,8,9]
test_cases[("SDP","case9",3,2)]["expected_values"] = Dict{String,Float64}("FP"=>0.0, "Minmax"=>3.4519, "Maxmin"=>3.4519)

# "ACR case 9"
test_cases[("ACR","case9",3,0)]=Dict{String,Any}("name"=>"case9K3ACR0", "attacker_budget"=>3)
test_cases[("ACR","case9",3,0)]["inactive_indices"]=[]
test_cases[("ACR","case9",3,0)]["protected_indices"]=[]
test_cases[("ACR","case9",3,0)]["expected_values"] = Dict{String,Float64}()

test_cases[("ACR","case9",3,1)]=Dict{String,Any}("name"=>"case9K3ACR1", "attacker_budget"=>0)
test_cases[("ACR","case9",3,1)]["inactive_indices"]=[1,4,7]
test_cases[("ACR","case9",3,1)]["protected_indices"]=[]
test_cases[("ACR","case9",3,1)]["expected_values"] = Dict{String,Float64}()

test_cases[("ACR","case9",3,2)]=Dict{String,Any}("name"=>"case9K3ACR2", "attacker_budget"=>3)
test_cases[("ACR","case9",3,2)]["inactive_indices"]=[]
test_cases[("ACR","case9",3,2)]["protected_indices"]=[2,3,5,6,8,9]
test_cases[("ACR","case9",3,2)]["expected_values"] = Dict{String,Float64}()

# "CASE 30"
# "SOC case 30"
test_cases[("SOC","case30",4,0)]=Dict{String,Any}("name"=>"case30K4SOC0", "attacker_budget"=>4)
test_cases[("SOC","case30",4,0)]["inactive_indices"]=[]
test_cases[("SOC","case30",4,0)]["protected_indices"]=[]
test_cases[("SOC","case30",4,0)]["expected_values"] = Dict{String,Float64}("FP"=>0.0, "Minmax"=>1.971, "Maxmin"=>0.937 )

test_cases[("SOC","case30",4,1)]=Dict{String,Any}("name"=>"case30K4SOC1", "attacker_budget"=>0)
test_cases[("SOC","case30",4,1)]["inactive_indices"]=[8,9,10,40]
test_cases[("SOC","case30",4,1)]["protected_indices"]=[]
test_cases[("SOC","case30",4,1)]["expected_values"] = Dict{String,Float64}("FP"=>0.937, "Minmax"=>0.937, "Maxmin"=>0.937)

test_cases[("SOC","case30",4,2)]=Dict{String,Any}("name"=>"case30K4SOC2", "attacker_budget"=>4)
test_cases[("SOC","case30",4,2)]["inactive_indices"]=[]
test_cases[("SOC","case30",4,2)]["protected_indices"]=filter(l->!(l in [8,9,10,40]), 1:41)
test_cases[("SOC","case30",4,2)]["expected_values"] = Dict{String,Float64}("FP"=>0.0, "Minmax"=>0.937, "Maxmin"=>0.937 )

# "SDP case 30"
test_cases[("SDP","case30",4,0)]=Dict{String,Any}("name"=>"case30K4SDP0", "attacker_budget"=>4)
test_cases[("SDP","case30",4,0)]["inactive_indices"]=[]
test_cases[("SDP","case30",4,0)]["protected_indices"]=[]
test_cases[("SDP","case30",4,0)]["expected_values"] = Dict{String,Float64}("FP"=>0.0, "Minmax"=>2.00117, "Maxmin"=>0.937 )

test_cases[("SDP","case30",4,1)]=Dict{String,Any}("name"=>"case30K4SDP1", "attacker_budget"=>0)
test_cases[("SDP","case30",4,1)]["inactive_indices"]=[8,9,10,40]
test_cases[("SDP","case30",4,1)]["protected_indices"]=[]
test_cases[("SDP","case30",4,1)]["expected_values"] = Dict{String,Float64}("FP"=>0.937, "Minmax"=>0.937, "Maxmin"=>0.937 )

test_cases[("SDP","case30",4,2)]=Dict{String,Any}("name"=>"case30K4SDP2", "attacker_budget"=>4)
test_cases[("SDP","case30",4,2)]["inactive_indices"]=[]
test_cases[("SDP","case30",4,2)]["protected_indices"]=filter(l->!(l in [8,9,10,40]), 1:41)
test_cases[("SDP","case30",4,2)]["expected_values"] = Dict{String,Float64}("FP"=>0.0, "Minmax"=>0.937, "Maxmin"=>0.937 )

# "ACR case 30"
test_cases[("ACR","case30",4,0)]=Dict{String,Any}("name"=>"case30K4ACR0", "attacker_budget"=>4)
test_cases[("ACR","case30",4,0)]["inactive_indices"]=[]
test_cases[("ACR","case30",4,0)]["protected_indices"]=[]
test_cases[("ACR","case30",4,0)]["expected_values"] = Dict{String,Float64}( )

test_cases[("ACR","case30",4,1)]=Dict{String,Any}("name"=>"case30K4ACR1", "attacker_budget"=>0)
test_cases[("ACR","case30",4,1)]["inactive_indices"]=[8,9,10,40]
test_cases[("ACR","case30",4,1)]["protected_indices"]=[]
test_cases[("ACR","case30",4,1)]["expected_values"] = Dict{String,Float64}( )

test_cases[("ACR","case30",4,2)]=Dict{String,Any}("name"=>"case30K4ACR2", "attacker_budget"=>4)
test_cases[("ACR","case30",4,2)]["inactive_indices"]=[]
test_cases[("ACR","case30",4,2)]["protected_indices"]=filter(l->!(l in [8,9,10,40]), 1:41)
test_cases[("ACR","case30",4,2)]["expected_values"] = Dict{String,Float64}( )

casesConic = [
	Dict("file" => "../data/case57.m", 
	 	"name" => "case57K2", 
	 	"attack_budget" => 2,
	 	"inactive_indices" => [],
	 	"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>2.768, "SOC_Maxmin"=>0.759737, "SDP_FP"=>0.0, "SDP_Minmax"=>3.06229, "SDP_Maxmin"=>0.77726 )
	),
	Dict("file" => "../data/case57.m", 
		"name" => "case57K0Inactive", 
		"attack_budget" => 0,
		"inactive_indices" => [41, 80],
		"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.759737, "SOC_Minmax"=>0.759737, "SOC_Maxmin"=>0.759737, "SDP_FP"=>0.0, "SDP_Minmax"=>0.77726, "SDP_Maxmin"=>0.77726 )
	),
	Dict("file" => "../data/case57.m", 
	 	"name" => "case57K4", 
	 	"attack_budget" => 4,
	 	"inactive_indices" => [],
	 	"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>4.511, "SOC_Maxmin"=>1.5793, "SDP_FP"=>0.0, "SDP_Minmax"=>4.619, "SDP_Maxmin"=>1.6 )
	),
	Dict("file" => "../data/case118.m", 
	 	"name" => "case118K2", 
	 	"expectedvalue" => 2.645, 
	 	"attack_budget" => 2,
	 	"inactive_indices" => [],
	 	"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>2.645, "SOC_Maxmin"=>1.447, "SDP_FP"=>0.0, "SDP_Minmax"=>2.70888, "SDP_Maxmin"=>1.448 )
	),
	Dict("file" => "../data/case118.m", 
	 	"name" => "case118K4", 
	 	"expectedvalue" => 4.051, 
	 	"attack_budget" => 4,
	 	"inactive_indices" => [],
	 	"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>4.051, "SOC_Maxmin"=>2.567, "SDP_FP"=>0.0, "SDP_Minmax"=>4.296, "SDP_Maxmin"=>2.568 )
	),
	Dict("file" => "../data/case300.m", 
	 	"name" => "case300K1", 
	 	"expectedvalue" => 0, 
	 	"attack_budget" => 1,
	 	"inactive_indices" => [],
	 	"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>0, "SOC_Maxmin"=>0, "SDP_FP"=>0.0, "SDP_Minmax"=>0, "SDP_Maxmin"=>0 )
	),
	Dict("file" => "../data/case300.m", 
	 	"name" => "case300K2", 
	 	"expectedvalue" => 0, 
	 	"attack_budget" => 2,
	 	"inactive_indices" => [],
	 	"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>0, "SOC_Maxmin"=>0, "SDP_FP"=>0.0, "SDP_Minmax"=>0, "SDP_Maxmin"=>0 )
	),
	Dict("file" => "../data/case300.m", 
	 	"name" => "case300K3", 
	 	"expectedvalue" => 0, 
	 	"attack_budget" => 3,
	 	"inactive_indices" => [],
	 	"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>0, "SOC_Maxmin"=>0, "SDP_FP"=>0.0, "SDP_Minmax"=>0, "SDP_Maxmin"=>0 )
	),
	Dict("file" => "../data/case300.m", 
	 	"name" => "case300K4", 
	 	"expectedvalue" => 0, 
	 	"attack_budget" => 4,
	 	"inactive_indices" => [],
	 	"protected_indices" => [],
		"expected_values" => Dict{String,Float64}("SOC_FP"=>0.0, "SOC_Minmax"=>0, "SOC_Maxmin"=>0, "SDP_FP"=>0.0, "SDP_Minmax"=>0, "SDP_Maxmin"=>0 )
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
