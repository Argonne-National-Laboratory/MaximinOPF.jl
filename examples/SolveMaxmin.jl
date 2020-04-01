include("../../MaximinOPF/src/MaximinOPF.jl")
using PowerModels
using JuMP
using Mosek
using MosekTools
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
 
soc_solver=with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0)

for pm_form in [SOCWRConicPowerModel] ### MaximinOPF.conic_supported_pm
    if !(pm_form in MaximinOPF.sdp_pm)
        start_time = time_ns()
        println("Formulating and solving the form ",pm_form)
        maxmin_model=MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=true)
        JuMP.set_optimizer(maxmin_model,soc_solver)
        JuMP.optimize!(maxmin_model)
        branch_ids=sort(collect(pm_data["undecided_branches"]))
        println("Optimal value using powerform ", pm_form, " is: ",JuMP.objective_value(maxmin_model), " with status ",JuMP.termination_status(maxmin_model))
        x_soln_str = ""
        x_soln = Dict{Int64,Float64}()
        for l in branch_ids
            x_var = variable_by_name(maxmin_model,"x[$l]_1")
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
        println("Final best solution: ",x_soln_str," with value ",JuMP.objective_value(maxmin_model))
        println("Runtime: ",runtime)
    end
end


