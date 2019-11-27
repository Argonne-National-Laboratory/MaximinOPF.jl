# "Case9 feas problem with no lines cut should have value 0.0"
# "Case9 feas problem with lines 1,4,7 cut should have value 3.451901668361513"
# "Case9 feas problem with lines 1--9 cut should have value 4.6"
# "Case30 feas problem with no lines cut should have value 0.0"
# "Case30 feas problem with lines 8,9,10,40 cut should have value 0.937"
# "Case30 feas problem with all lines cut should have value 2.57327375"
# "Case57 feas problem with no lines cut should have value 0.0"
# "Case57 feas problem with lines 41,80 cut should have value 0.7597369902683009"
# "Case57 feas problem with all lines cut should have value 6.2313456"

cases = [
	Dict("file" => "../data/case9.m", 
		"name" => "case9", 
		"expectedvalue" => 3.451901668361513, 
		"cutindex" => [1, 4, 7]
	),
	Dict("file" => "../data/case30.m", 
		"name" => "case30", 
		"expectedvalue" => 0.937, 
		"cutindex" => [8, 9, 10, 40]
	),
	Dict("file" => "../data/case57.m", 
		"name" => "case57", 
		"expectedvalue" => 0.7597369902683009, 
		"cutindex" => [41, 80]
	)
]


function getTestcases()
	pm_datas = []

	for i in 1:length(cases)
	#Set Default Input
		case = cases[i]["file"]
		casename = cases[i]["name"]
		expect = cases[i]["expectedvalue"]
		lineindexs = cases[i]["cutindex"]
		K = length(lineindexs)
		pm_data = PowerModels.parse_file(case)
	    for (k,line) in pm_data["branch"]
			if line["index"] in lineindexs
				line["br_status"]=0
			end			
	    end		
	    push!(pm_datas, 
	    	Dict("pm_data" => pm_data, 
	    		"K" => K, 
	    		"expectedvalue" => expect, 
	    		"name" => casename,
	    		"cutindex" => lineindexs
	    	))
	end
	return pm_datas
end
