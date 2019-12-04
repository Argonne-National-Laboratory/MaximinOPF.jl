# "Case9 feas problem with no lines cut should have value 0.0"
# "Case9 feas problem with lines 1,4,7 cut should have value 3.451901668361513"
# "Case9 feas problem with lines 1--9 cut should have value 4.6"
# "Case30 feas problem with no lines cut should have value 0.0"
# "Case30 feas problem with lines 8,9,10,40 cut should have value 0.937"
# "Case30 feas problem with all lines cut should have value 2.57327375"
# "Case57 feas problem with no lines cut should have value 0.0"
# "Case57 feas problem with lines 41,80 cut should have value 0.7597369902683009"
# "Case57 feas problem with all lines cut should have value 6.2313456"

casesFP = [
	Dict("file" => "../data/case9.m", 
		"name" => "case9", 
		"expectedvalue" => 3.451901668361513, 
		"inactive_indices" => [1, 4, 7]
	),
	Dict("file" => "../data/case30.m", 
		"name" => "case30", 
		"expectedvalue" => 0.937, 
		"inactive_indices" => [8, 9, 10, 40]
	),
	Dict("file" => "../data/case57.m", 
		"name" => "case57", 
		"expectedvalue" => 0.7597369902683009, 
		"inactive_indices" => [41, 80]
	)
]

casesRMinmax = [
	Dict("file" => "../data/case9.m", 
		"name" => "case9", 
		"expectedvalue" => 3.4519, 
		"attack_budget" => 3,
		"inactive_indices" => [1,4,7],
		"protected_indices" => []
	),
	Dict("file" => "../data/case9.m", 
		"name" => "case9", 
		"expectedvalue" => 3.4519, 
		"attack_budget" => 3,
		"inactive_indices" => [],
		"protected_indices" => [2,3,5,6,8,9]
	),
	Dict("file" => "../data/case30.m", 
		"name" => "case30", 
		"expectedvalue" => 1.971, 
		"attack_budget" => 4,
		"inactive_indices" => [],
		"protected_indices" => []
	),
	Dict("file" => "../data/case57.m", 
		"name" => "case57", 
		"expectedvalue" => 2.768, 
		"attack_budget" => 2,
		"inactive_indices" => [],
		"protected_indices" => []
	),
	Dict("file" => "../data/case57.m", 
		"name" => "case57", 
		"expectedvalue" => 4.511, 
		"attack_budget" => 4,
		"inactive_indices" => [],
		"protected_indices" => []
	),
	Dict("file" => "../data/case118.m", 
		"name" => "case118", 
		"expectedvalue" => 2.645, 
		"attack_budget" => 2,
		"inactive_indices" => [],
		"protected_indices" => []
	),
	Dict("file" => "../data/case118.m", 
		"name" => "case118", 
		"expectedvalue" => 4.051, 
		"attack_budget" => 4,
		"inactive_indices" => [],
		"protected_indices" => []
	),
]


function getTestcasesFP()
	pm_datas = []

	for i in 1:length(casesFP)
	#Set Default Input
		case = casesFP[i]["file"]
		casename = casesFP[i]["name"]
		expect = casesFP[i]["expectedvalue"]
		inactive_branches = casesFP[i]["inactive_indices"]
		protected_branches = casesFP[i]["protected_index"]
		K = length(attack_budget)
		pm_data = PowerModels.parse_file(case)
		pm_data["attacker_budget"]=K ###Adding another key and entry
		pm_data["inactive_branches"]=inactive_branches ###Adding another key and entry
		pm_data["protected_branches"]=protected_branches ###Adding another key and entry
	        for (k,line) in pm_data["branch"]
			if line["index"] in inactive_branches
				line["br_status"]=0
			end			
	        end		
	    push!(pm_datas, 
	    	Dict("pm_data" => pm_data, 
	    		"K" => K, 
	    		"expectedvalue" => expect, 
	    		"name" => casename,
	    		"inactive_indices" => inactive_branches
	    	))
	end
	return pm_datas
end

function getTestcasesRMinmax()
	pm_datas = []

	for i in 1:length(casesRMinmax)
	#Set Default Input
		case = casesRMinmax[i]["file"]
		casename = casesRMinmax[i]["name"]
		expect = casesRMinmax[i]["expectedvalue"]
		protected_branches = casesRMinmax[i]["protected_indices"]
		inactive_branches = casesRMinmax[i]["inactive_indices"]
		K = casesRMinmax[i]["attack_budget"]-length(inactive_branches)
		pm_data = PowerModels.parse_file(case)
		pm_data["attacker_budget"]=K ###Adding another key and entry
		pm_data["protected_branches"]=protected_branches ###Adding another key and entry
		pm_data["inactive_branches"]=inactive_branches ###Adding another key and entry
	    push!(pm_datas, 
	    	Dict("pm_data" => pm_data, 
	    		"K" => K, 
	    		"expectedvalue" => expect, 
	    		"name" => casename,
	    		"protected_indices" => protected_branches,
	    		"inactive_indices" => inactive_branches
	    	))
	end
	return pm_datas
end
