using PowerModels


testcases = [
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SOCWRConicPowerModel,
 	"name" => "case9SOCWR1",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => []
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SOCWRConicPowerModel,
 	"name" => "case9SOCWR2",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [1,2,3,4,5,6,7,8,9],
 	"protected_indices" => []
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SOCWRConicPowerModel,
 	"name" => "case9SOCWR3",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => [1,2,3,4,5,6,7,8,9]
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SOCWRConicPowerModel,
 	"name" => "case9SOCWR4",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [1,4],
 	"protected_indices" => []
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SOCWRConicPowerModel,
 	"name" => "case9SOCWR5",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => [2,3]
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SOCWRConicPowerModel,
 	"name" => "case9SOCWR6",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [1,4,7],
 	"protected_indices" => [2,3,5,6,8,9]
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SOCWRConicPowerModel,
 	"name" => "case9SOCWR7",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [1],
 	"protected_indices" => [2]
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SDPWRMPowerModel,
 	"name" => "case9SDPWR1",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => []
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SDPWRMPowerModel,
 	"name" => "case9SDPWR2",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [1,2,3,4,5,6,7,8,9],
 	"protected_indices" => []
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SDPWRMPowerModel,
 	"name" => "case9SDPWR3",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => [1,2,3,4,5,6,7,8,9]
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SDPWRMPowerModel,
 	"name" => "case9SDPWR4",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [1,4],
 	"protected_indices" => []
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SDPWRMPowerModel,
 	"name" => "case9SDPWR5",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => [2,3]
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SDPWRMPowerModel,
 	"name" => "case9SDPWR6",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [1,4,7],
 	"protected_indices" => [2,3,5,6,8,9]
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SDPWRMPowerModel,
 	"name" => "case9SDPWR7",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [1],
 	"protected_indices" => [2]
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SparseSDPWRMPowerModel,
 	"name" => "case9SparseSDPWR1",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => []
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SparseSDPWRMPowerModel,
 	"name" => "case9SparseSDPWR2",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [1,2,3,4,5,6,7,8,9],
 	"protected_indices" => []
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SparseSDPWRMPowerModel,
 	"name" => "case9SparseSDPWR3",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => [1,2,3,4,5,6,7,8,9]
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SparseSDPWRMPowerModel,
 	"name" => "case9SparseSDPWR4",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [1,4],
 	"protected_indices" => []
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SparseSDPWRMPowerModel,
 	"name" => "case9SparseSDPWR5",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => [2,3]
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SparseSDPWRMPowerModel,
 	"name" => "case9SparseSDPWR6",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [1,4,7],
 	"protected_indices" => [2,3,5,6,8,9]
	),
	Dict(
	"file" => "../data/case9.m", 
	"PMOption" => SparseSDPWRMPowerModel,
 	"name" => "case9SparseSDPWR7",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [1],
 	"protected_indices" => [2]
	),
	Dict(
	"file" => "../data/case30.m", 
	"PMOption" => SparseSDPWRMPowerModel,
 	"name" => "case30K4SparseSDPWR0",  	
 	"attack_budget" => 4,
 	"inactive_indices" => [],
 	"protected_indices" => []
	),
	Dict(
	"file" => "../data/case30.m", 
	"PMOption" => SparseSDPWRMPowerModel,
 	"name" => "case30K4SparseSDPWR1",  	
 	"attack_budget" => 0,
 	"inactive_indices" => [8,9,10,40],
 	"protected_indices" => []
	)
]
