include("../../MaximinOPF/src/MaximinOPF.jl")
using PowerModels
using JuMP
using SCS
using Ipopt
using ProxSDP
PowerModels.silence()

testcase = Dict(
	"file" => "data/case9.m", 
	"PMOption" => SparseSDPWRMPowerModel,
 	"name" => "case9K3",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => []
	)

pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry

sdp_solver=with_optimizer(ProxSDP.Optimizer, log_verbose=false, tol_primal=1e-6, tol_dual=1e-6 )
conic_solver=with_optimizer(SCS.Optimizer,verbose=0)
ip_solver=with_optimizer(Ipopt.Optimizer,print_level=0)

# nonconvex AC forms
nonconvex_ac=[ACPPowerModel, ACRPowerModel, ACTPowerModel]
for pm_form in nonconvex_ac
    println("Formulating and solving the form ",pm_form)
    pm=MaximinOPF.MinimaxOPFModel(pm_data, pm_form)
    JuMP.set_optimizer(pm.model,ip_solver)
    JuMP.optimize!(pm.model)
    println("\tOptimal value using powerform ", pm_form, " is: ",JuMP.objective_value(pm.model), " with status ",JuMP.termination_status(pm.model))
end

# linear approximations
linear_approx=[DCPPowerModel, DCMPPowerModel, NFAPowerModel]
for pm_form in linear_approx
    println("Formulating and solving the form ",pm_form)
    pm=MaximinOPF.MinimaxOPFModel(pm_data, pm_form)
    JuMP.set_optimizer(pm.model,conic_solver)
    JuMP.optimize!(pm.model)
    println("\tOptimal value using powerform ", pm_form, " is: ",JuMP.objective_value(pm.model), " with status ",JuMP.termination_status(pm.model))
end

# quadratic approximations
quadratic_approx=[DCPLLPowerModel, LPACCPowerModel]
for pm_form in quadratic_approx
    println("Formulating and solving the form ",pm_form)
    pm=MaximinOPF.MinimaxOPFModel(pm_data, pm_form)
    JuMP.set_optimizer(pm.model,ip_solver)
    JuMP.optimize!(pm.model)
    println("\tOptimal value using powerform ", pm_form, " is: ",JuMP.objective_value(pm.model), " with status ",JuMP.termination_status(pm.model))
end
# quadratic relaxations
quadratic_relax=[SOCWRPowerModel, SOCBFPowerModel, QCRMPowerModel, QCLSPowerModel]
for pm_form in quadratic_relax
    println("Formulating and solving the form ",pm_form)
    pm=MaximinOPF.MinimaxOPFModel(pm_data, pm_form)
    JuMP.set_optimizer(pm.model,ip_solver)
    JuMP.optimize!(pm.model)
    println("\tOptimal value using powerform ", pm_form, " is: ",JuMP.objective_value(pm.model), " with status ",JuMP.termination_status(pm.model))
end
quad_conic_relax=[SOCWRConicPowerModel, SOCBFConicPowerModel]
for pm_form in quad_conic_relax
    println("Formulating and solving the form ",pm_form)
    pm=MaximinOPF.MinimaxOPFModel(pm_data, pm_form)
    JuMP.set_optimizer(pm.model,conic_solver)
    JuMP.optimize!(pm.model)
    println("\tOptimal value using powerform ", pm_form, " is: ",JuMP.objective_value(pm.model), " with status ",JuMP.termination_status(pm.model))
end
# sdp relaxations
sdp_relax=[SDPWRMPowerModel, SparseSDPWRMPowerModel]
for pm_form in sdp_relax
    println("Formulating and solving the form ",pm_form)
    pm=MaximinOPF.MinimaxOPFModel(pm_data, pm_form)
    JuMP.set_optimizer(pm.model,sdp_solver)
    JuMP.optimize!(pm.model)
    println("\tOptimal value using powerform ", pm_form, " is: ",JuMP.objective_value(pm.model), " with status ",JuMP.termination_status(pm.model))
end

