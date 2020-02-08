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
	"file" => "data/case57.m", 
 	"name" => "case57K2",  	
 	"attack_budget" => 2,
 	"inactive_indices" => [],
 	"protected_indices" => []
	)

function solveMaxmin(pm_data,form,optimizer)
  model = MaximinOPF.MaximinOPFModel(pm_data,form)
  if optimizer==Mosek.Optimizer
    JuMP.set_optimizer(model,with_optimizer(optimizer,MSK_IPAR_LOG=0))
  else
    JuMP.set_optimizer(model,with_optimizer(optimizer))
  end
  JuMP.optimize!(model)
  status=JuMP.termination_status(model)
  if status != OPTIMAL
    println("FLAGGING: Solve status=",status)
  end
  return model 
end #end of function

pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry

io = open("output.txt", "w")
for pm_form in MaximinOPF.conic_supported_pm
  if !(pm_form in MaximinOPF.sdp_pm)
    println("Formulating and solving the form ",pm_form)
    model=solveMaxmin(pm_data,pm_form, Mosek.Optimizer)
    println(io,"Optimal value using powerform ", pm_form, " is: ",JuMP.objective_value(model), " with status ",JuMP.termination_status(model))
  end
end
close(io)

