include("../../MaximinOPF/src/MaximinOPF.jl")
using JuMP
using PowerModels
using SCS
using ProxSDP
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

### "sdp_pm=[ SDPWRMPowerModel, SparseSDPWRMPowerModel ]"

conic_solver=with_optimizer(SCS.Optimizer,verbose=0,max_iters=30000)
sdp_solver=with_optimizer(SCS.Optimizer,verbose=0,max_iters=30000)
sparse_sdp_solver=with_optimizer(SCS.Optimizer,verbose=0,max_iters=30000)
#sdp_solver=with_optimizer(ProxSDP.Optimizer, max_iter=100000, log_verbose=false, initial_beta=0.005, full_eig_decomp=false )
#sparse_sdp_solver=with_optimizer(ProxSDP.Optimizer, max_iter=100000,log_verbose=false,initial_beta=0.005,full_eig_decomp=true)

for pm_form in MaximinOPF.conic_supported_pm
    println("\nFormulating and solving the dual minmax problem for the form ",pm_form)
    dual_minmax_model=MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=false)
    if pm_form == SDPWRMPowerModel
        JuMP.set_optimizer(dual_minmax_model,sdp_solver)
    elseif pm_form == SparseSDPWRMPowerModel
        JuMP.set_optimizer(dual_minmax_model,sparse_sdp_solver)
    else 
        JuMP.set_optimizer(dual_minmax_model,conic_solver)
    end
    JuMP.optimize!(dual_minmax_model)
    println("\tOptimal value using powerform ", pm_form, " is: ",JuMP.objective_value(dual_minmax_model), " with status ",JuMP.termination_status(dual_minmax_model))
end

