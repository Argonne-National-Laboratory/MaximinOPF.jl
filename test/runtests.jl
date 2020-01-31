using MaximinOPF
using Ipopt
using PowerModels
using Mosek
using MosekTools
using JuMP
#using SDPA

include("testcases.jl")

PowerModels.silence()

supportedPMOptions = Dict{String,Dict{String,Any}}(
	"SOC"=> Dict{String,Any}("FORM"=>SOCWRConicPowerModel,"OPT"=>Mosek.Optimizer),
	"SDP"=> Dict{String,Any}("FORM"=>SparseSDPWRMPowerModel,"OPT"=>Mosek.Optimizer),
	"ACR"=> Dict{String,Any}("FORM"=>ACRPowerModel,"OPT"=>Ipopt.Optimizer)
)

"
#	SOCWRConicPowerModel, # MoSek
#	SDPWRMPowerModel, # MoSek
	SparseSDPWRMPowerModel, # MoSek, requires minor editing of PowerModels wrm.jl code as of 5 Jan 2020, open issue with PowerModels is in progress.
#	SOCWRPowerModel, # Not Mosek	
#	QCRMPowerModel # Not Mosek
#	SOCBFPowerModel, # Error constraint_ohms_yt_from()	
"
function evaluateModelOptVals(expected_value, ev_type, model, tol, io=Base.stdout)
	passes=true
	if !haskey(expected_value,ev_type)
		result_str=string("UNKNOWN expected value. Obj value otherwise is: ",JuMP.objective_value(model),".")
		passes=true
    	else
		discrep = abs(JuMP.objective_value(model)-expected_value)
    		if discrep < tol
			result_str=string("PASSED: Obj value is: ",JuMP.objective_value(model),", and the expected value was: ",expected_value,".")
			passes=true
    		else
			result_str=string("FAILED: Obj value is: ",JuMP.objective_value(model),", but the expected value was: ",expected_value,".")
			passes=false
    		end
    	end
	println(io,result_str)
	return passes,result_str
end

for k in keys(test_cases)
  powerform = supportedPMOptions[k[1]]["FORM"] #PowerModel Options

  pm_data = PowerModels.parse_file( string("../data/",k[2],".m") )
  pm_data["name"]=test_cases[k]["name"]
  pm_data["attacker_budget"] = test_cases[k]["attacker_budget"] ###Adding another key and entry to pm_data
  pm_data["inactive_branches"] = test_cases[k]["inactive_indices"] ###Adding another key and entry to pm_data
  pm_data["protected_branches"] = test_cases[k]["protected_indices"] ###Adding another key and entry to pm_data

  f_name = string(pm_data["name"],".txt")
  io=open(f_name, "w")

  #Create JUMP Model
    #pf_model_pm = MaximinOPF.PF_FeasModel(pm_data, powerform)
    pf_model_pm = MaximinOPF.MinimaxOPFModel(pm_data, powerform)
    pm_data["undecided_branches"]= filter(l->!(l in pm_data["protected_branches"] || l in pm_data["inactive_branches"]), ids(pf_model_pm,pf_model_pm.cnw,:branch)) 
			###Adding another key and entry
    pf_minmax_model = pf_model_pm.model
    if k[1] != "ACR"
      pf_maxmin_model = MaximinOPF.Minmax_to_Maxmin(pf_model_pm)
     #MaximinOPF.write_to_cbf(pf_maxmin_model,test_cases[k]["name"])
      MaximinOPF.write_to_cbf_scip(pf_maxmin_model,test_cases[k]["name"])
    end


  #Print Model Status		
    println(io,"inactive_branches: ",pm_data["inactive_branches"])
    println(io,"protected_branches: ",pm_data["protected_branches"])
    println(io,"undecided_branches: ",pm_data["undecided_branches"])
    println(io,"Print Model")
    println(io,"\n**********BEGIN MINMAX MODEL**********")
    println(io,pf_minmax_model)
    println(io,"\n**********END MINMAX MODEL**********")
    #println(io,"\n**********BEGIN MAXMIN MODEL**********")
    #println(io,pf_maxmin_model)
    #println(io,"\n**********END MAXMIN MODEL**********")

  #Solve Model with PowerModels Solution Builder
    opt_model = pf_minmax_model
    opt_expected = "Minmax"
   #opt_model = pf_maxmin_model
    set_optimizer(opt_model,supportedPMOptions[k[1]]["OPT"])
    result = @elapsed JuMP.optimize!(opt_model)
    println("Solving the ",opt_expected," problem ",test_cases[k]["name"],"...")
    status=JuMP.termination_status(opt_model)
    println(io,"Time taken to solve is: ", result, " with status ",status,".")
    test_cases[k]["results"]=evaluateModelOptVals(test_cases[k]["expected_values"],opt_expected, opt_model ,1e-3,io)
    close(io)
end
		
println("Test results: ")
for k in keys(test_cases)
   println("Test ",test_cases[k]["name"]," ",test_cases[k]["results"])
end
