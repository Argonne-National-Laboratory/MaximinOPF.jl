using PowerModels

testcases = [
	Dict(
	"file" => "../data/case9.m", 
 	"name" => "case9-1",
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => []
	),
	Dict(
	"file" => "../data/case9.m", 
	"name" => "case9-2",
 	"attack_budget" => 3,
 	"inactive_indices" => [1,2,3,4,5,6,7,8,9],
 	"protected_indices" => []
	),
	Dict(
	"file" => "../data/case9.m", 
	"name" => "case9-3",
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => [1,2,3,4,5,6,7,8,9]
	),
	Dict(
	"file" => "../data/case9.m", 
	"name" => "case9-4",
 	"attack_budget" => 3,
 	"inactive_indices" => [1,4],
 	"protected_indices" => []
	),
	Dict(
	"file" => "../data/case9.m", 
	"name" => "case9-5",
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => [2,3]
	),
	Dict(
	"file" => "../data/case9.m", 
	"name" => "case9-6",
 	"attack_budget" => 3,
 	"inactive_indices" => [1,4,7],
 	"protected_indices" => [2,3,5,6,8,9]
	),
	Dict(
	"file" => "../data/case9.m", 
	"name" => "case9-7",
 	"attack_budget" => 3,
 	"inactive_indices" => [1],
 	"protected_indices" => [2]
	)
]
