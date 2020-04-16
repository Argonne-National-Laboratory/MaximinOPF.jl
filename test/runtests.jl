using MaximinOPF
using PowerModels
using Mosek
using MosekTools
using JuMP
#using SCIP

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
    maxmin_model = MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=true, io=devnull, rm_rsoc=true, rm_therm_line_lim=false)

    println(string("Start Solving: ", testcases[j]["name"]))

    if occursin("SOC", testcases[j]["name"])
        
        set_optimizer(maxmin_model,Mosek.Optimizer)  
        result = @elapsed JuMP.optimize!(maxmin_model)
        #Print Result   
        status=JuMP.termination_status(maxmin_model)
        println("Time taken to solve is: ", result, " with status ",status,".")     
    elseif occursin("SDP", testcases[j]["name"])
        #Convert MOI to CBF
        #Solve with SCIP
        # maxmin_model = MaximinOPF.SOCtoPSD(maxmin_model)
        MaximinOPF.write_to_cbf_scip(maxmin_model,testcases[j]["name"])

        # set_optimizer(maxmin_model,SCIP.Optimizer)  
        # result = @elapsed JuMP.optimize!(maxmin_model)
        # #Print Result   
        # status=JuMP.termination_status(maxmin_model)
        # println("Time taken to solve is: ", result, " with status ",status,".")     
        println(string("Solver is not selected for ", testcases[j]["name"]))
    else
        println(string("Solver is not selected for ", testcases[j]["name"]))

    end
end

