using MaximinOPF
using PowerModels
using JuMP
using SCS

include("testcases.jl")
PowerModels.silence()


for j in 1:length(testcases)
    pm_data = PowerModels.parse_file( testcases[j]["file"] )
    pm_form = testcases[j]["PMOption"]
    pm_data["name"]=testcases[j]["name"]
    pm_data["attacker_budget"] = testcases[j]["attack_budget"] ###Adding another key and entry
    pm_data["inactive_branches"] = testcases[j]["inactive_indices"] ###Adding another key and entry
    pm_data["protected_branches"] = testcases[j]["protected_indices"] ###Adding another key and entry
    pm_data["pm_form"] = pm_form

    #Create JUMP Model
    maxmin_model = MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=false) # "SCS does not have mixed-integer support"

    if occursin("SOC", testcases[j]["name"])
        println(string("Start Solving: ", testcases[j]["name"]))
        set_optimizer(maxmin_model,SCS.Optimizer)  
        result = @elapsed JuMP.optimize!(maxmin_model)
        #Print Result   
        status=JuMP.termination_status(maxmin_model)
        println("Time taken to solve is: ", result, " with status ",status,".")     
    elseif occursin("SDP", testcases[j]["name"])
        # "Prepare the input for scipsdp, convert MOI to CBF"
        println(string("Preparing cbf input file for the application of scipsdp to: ", testcases[j]["name"]))
        MaximinOPF.write_to_cbf_scip(maxmin_model,testcases[j]["name"])
    end
end

