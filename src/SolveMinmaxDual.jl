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

function solveNodeMinmaxDual(pm_data,form,optimizer)
  pm = MaximinOPF.MinimaxOPFModel(pm_data, form)
  model = MaximinOPF.DualizeMinmaxModel(pm)  
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
  println("Formulating and solving the form ",pm_form)
  model=solveNodeMinmaxDual(pm_data,pm_form, Mosek.Optimizer)
  println(io,"Optimal value using powerform ", pm_form, " is: ",JuMP.objective_value(model), " with status ",JuMP.termination_status(model))
end
close(io)

