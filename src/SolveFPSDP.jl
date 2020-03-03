include("../../MaximinOPF/src/MaximinOPF.jl")
using JuMP
using PowerModels
using ProxSDP
#using Ipopt
#using Mosek
#using MosekTools
#using SCIP
using LinearAlgebra
#using Arpack
#using Printf
PowerModels.silence()

testcase = Dict(
	"file" => "data/case30.m", 
	"PMOption" => SparseSDPWRMPowerModel,
 	"name" => "case30K4SDP",  	
 	"attack_budget" => 0,
 	"inactive_indices" => [8,9,10,40],
 	"protected_indices" => []
	)

pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry

io = open(string(testcase["name"],".out"), "w")

# sdp relaxations
sdp_relax=[SDPWRMPowerModel, SparseSDPWRMPowerModel]
for pm_form in sdp_relax
  println("Formulating and solving the form ",pm_form)
  pm=MaximinOPF.PF_FeasModel(pm_data, pm_form )
  set_optimizer(pm.model, with_optimizer(ProxSDP.Optimizer, log_verbose=true, tol_primal = 1e-6, tol_dual = 1e-6 ))
  optimize!(pm.model)
  println(io,"Optimal value using powerform ", pm_form, " is: ",JuMP.objective_value(pm.model), " with status ",JuMP.termination_status(pm.model))
  println(io,"Solve time is: ",JuMP.solve_time(pm.model) ) 
end
close(io)

