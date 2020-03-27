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
	"file" => "data/case300.m", 
 	"name" => "case300K4",  	
 	"attack_budget" => 4,
 	"inactive_indices" => [],
 	"protected_indices" => []
	)

function solveMaxmin(pm_data,pm_form,optimizer)
  model = MaximinOPF.MaximinOPFModel(pm_data,pm_form)
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

for pm_form in [SOCWRConicPowerModel] ### MaximinOPF.conic_supported_pm
  if !(pm_form in MaximinOPF.sdp_pm)
    io = open(string(testcase["name"],"SOC.out"), "w")
    start_time = time_ns()
    println("Formulating and solving the form ",pm_form)
    model=solveMaxmin(pm_data,pm_form, Mosek.Optimizer)
    branch_ids=sort(collect(pm_data["undecided_branches"]))
    println(io,"Optimal value using powerform ", pm_form, " is: ",JuMP.objective_value(model), " with status ",JuMP.termination_status(model))
    x_soln_str = ""
    x_soln = Dict{Int64,Float64}()
    for l in branch_ids
        x_var = variable_by_name(model,"x[$l]_1")
        x_val = JuMP.value(x_var)
        if x_val > 1.0-1.0e-8
	        x_val = 1
            x_soln_str = string(x_soln_str," $l")
        elseif x_val < 1.0e-8
	        x_val = 0
        end
        x_soln[l] = x_val
    end
    end_time = time_ns()
    runtime = (end_time-start_time)/1e9
    println(io,"Final best solution: ",x_soln_str," with value ",JuMP.objective_value(model))
    println(io,"Runtime: ",runtime)
    close(io)
  end
end


