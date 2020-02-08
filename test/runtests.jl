using MaximinOPF
using PowerModels
using Mosek
using MosekTools
using JuMP

include("testcases.jl")
PowerModels.silence()


for j in 1:length(testcases)
    pm_data = PowerModels.parse_file( testcases[j]["file"] )
    powerform = testcases[j]["PMOption"]
    pm_data["name"]=testcases[j]["name"]
    pm_data["attacker_budget"] = testcases[j]["attack_budget"] ###Adding another key and entry
    pm_data["inactive_branches"] = testcases[j]["inactive_indices"] ###Adding another key and entry
    pm_data["protected_branches"] = testcases[j]["protected_indices"] ###Adding another key and entry

    #Create JUMP Model
    maxmin_model = MaximinOPF.MaximinOPFModel(pm_data, powerform) 
    if occursin("SparseSDP", testcases[j]["name"])
      MaximinOPF.write_to_cbf_scip(maxmin_model,pm_data["name"])
    end

    println(string("Start Solving: ", testcases[j]["name"]))

    if occursin("SOC", testcases[j]["name"])
        
        set_optimizer(maxmin_model,with_optimizer(Mosek.Optimizer))  
        result = @elapsed JuMP.optimize!(maxmin_model)
        #Print Result   
        status=JuMP.termination_status(maxmin_model)
        println("Time taken to solve is: ", result, " with status ",status,".")     
    else
        println(string("Solver is not selected for ", testcases[j]["name"]))

    end
end

