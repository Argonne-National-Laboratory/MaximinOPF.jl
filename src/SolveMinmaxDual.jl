include("../../MaximinOPF/src/MaximinOPF.jl")
using JuMP
using PowerModels
using Ipopt
using Mosek
using MosekTools
#using SCIP
using LinearAlgebra
using Arpack
using Printf
PowerModels.silence()

testcase = Dict(
	"file" => "data/case30.m", 
 	"name" => "case30K4",  	
 	"attack_budget" => 4,
 	"inactive_indices" => [],
 	"protected_indices" => []
	)

pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry

io = open("output.txt", "w")
for pm_form in MaximinOPF.conic_supported_pm
  println("Formulating and solving the form ",pm_form)
  model,pm=MaximinOPF.SolveMinmaxDual(pm_data,pm_form, with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0))
  println(io,"Optimal value using powerform ", pm_form, " is: ",JuMP.objective_value(model), " with status ",JuMP.termination_status(model))
end
close(io)

