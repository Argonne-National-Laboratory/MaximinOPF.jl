using MaximinOPF
using Ipopt
using PowerModels
using Mosek
using MosekTools
using JuMP

include("testcases.jl")

supportedPMOptions = [	
#	SOCWRConicPowerModel, # MoSek
#	SDPWRMPowerModel, # MoSek
	SparseSDPWRMPowerModel # MoSek, requires minor editing of PowerModels wrm.jl code as of 5 Jan 2020, open issue with PowerModels is in progress.
#	SOCWRPowerModel, # Not Mosek	
#	QCRMPowerModel # Not Mosek
#	SOCBFPowerModel, # Error constraint_ohms_yt_from()	
]

function evaluateModel(expectedvalue, model, tol)
    println("Protected branches: ",model.data["protected_branches"])
    println("Inactive branches: ",model.data["inactive_branches"])
    println("All other branches: ",setdiff(ids(model, :branch),model.data["protected_branches"]))
    lk = []
    println("Printing dual values x:")
    for l in setdiff(ids(model, :branch),model.data["protected_branches"])
      if abs(JuMP.dual(con(model, model.cnw, model.ccnd)[:x][l])) > tol
        push!(lk, l)	        
      end
    end
    println(lk)

    lkx = []
    psival = model.data["attacker_budget"]*JuMP.value(var(model, model.cnw, model.ccnd)[:u_K])
    psival += sum( JuMP.value(var(model, model.cnw, model.ccnd)[:u_ord_aux][l]) for l in setdiff(ids(model, :branch),model.data["protected_branches"])  )
    println("uK: ",JuMP.value(var(model, model.cnw, model.ccnd)[:u_K])," psival: ",psival)
    for l in setdiff(ids(model, :branch),model.data["protected_branches"])
      branch = ref(model, model.cnw, :branch, l)
      f_bus = branch["f_bus"]
      t_bus = branch["t_bus"]
      f_idx = (l, f_bus, t_bus)
      t_idx = (l, t_bus, f_bus)
      upf0 = JuMP.value(var(model, model.cnw, model.ccnd,:up_br)[f_idx,0])
      upt0 = JuMP.value(var(model, model.cnw, model.ccnd,:up_br)[t_idx,0])
      uqf0 = JuMP.value(var(model, model.cnw, model.ccnd,:uq_br)[f_idx,0])
      uqt0 = JuMP.value(var(model, model.cnw, model.ccnd,:uq_br)[t_idx,0])
      upf1 = JuMP.value(var(model, model.cnw, model.ccnd,:up_br)[f_idx,1])
      upt1 = JuMP.value(var(model, model.cnw, model.ccnd,:up_br)[t_idx,1])
      uqf1 = JuMP.value(var(model, model.cnw, model.ccnd,:uq_br)[f_idx,1])
      uqt1 = JuMP.value(var(model, model.cnw, model.ccnd,:uq_br)[t_idx,1])
      push!(lkx, (l,upf0+upt0+uqf0+uqt0,-upf1-upt1-uqf1-uqt1,upf0+upt0+uqf0+uqt0-upf1-upt1-uqf1-uqt1))	        
      println((l,upf0+upt0+uqf0+uqt0,-upf1-upt1-uqf1-uqt1,upf0+upt0+uqf0+uqt0-upf1-upt1-uqf1-uqt1))	        
    end
    # return 1
    
end

testresults = []
for i in 1:length(supportedPMOptions)
	powerform = supportedPMOptions[i] #PowerModel Options
	for j in 1:length(casesConic)
		pm_data = PowerModels.parse_file( casesConic[j]["file"] )
		pm_data["name"]=casesConic[j]["name"]
		pm_data["attacker_budget"] = casesConic[j]["attack_budget"] ###Adding another key and entry
		pm_data["inactive_branches"] = casesConic[j]["inactive_indices"] ###Adding another key and entry
		pm_data["protected_branches"] = casesConic[j]["protected_indices"] ###Adding another key and entry


		f_name_base="output_"
		f_name = string(f_name_base,"_",casesConic[j]["name"],"_", pm_data["attacker_budget"],"_",j,".txt")
		io=open(f_name, "w")
		
		#Create JUMP Model
		#pf_model_pm = MaximinOPF.PF_FeasModel(pm_data, powerform)
		pf_model_pm = MaximinOPF.MinimaxOPFModel(pm_data, powerform)
		pm_data["undecided_branches"]= filter(l->!(l in pm_data["protected_branches"] || l in pm_data["inactive_branches"]), ids(pf_model_pm,pf_model_pm.cnw,:branch)) 
			###Adding another key and entry
		#pf_model = pf_model_pm.model
		pf_model = MaximinOPF.DualizeModel(pf_model_pm)
		
		#Print Model Status		
		println(io,"inactive_branches: ",pm_data["inactive_branches"])
		println(io,"protected_branches: ",pm_data["protected_branches"])
		println(io,"undecided_branches: ",pm_data["undecided_branches"])
		println(io,"Print Model")
		println(io,pf_model)

		#Solve Model with PowerModels Solution Builder
#=
		println("Start Solving")
		if i > 3 
			result = @elapsed JuMP.optimize!(pf_model,with_optimizer(Ipopt.Optimizer))
		else
			result = @elapsed JuMP.optimize!(pf_model,with_optimizer(Mosek.Optimizer))
		end
		
		#Print Result
		#evaluateModel(expect, pf_model_pm, 0.001)
		status=JuMP.termination_status(pf_model)
		println(io,"Time taken to solve is: ", result, " with status ",status,".")
		println(io,"The optimal value is: ",JuMP.objective_value(pf_model)," versus the expected value of ",casesConic[j]["expected_values"]["SDP_Minmax"])
=#
		close(io)
		
	end
end

println(testresults)
println("Printing results of the tests:")
for kk=1:length(testresults)
  print("Test ",kk," ")
  if testresults[kk]["testresult"]==1
    println("passed.")
  else
    println("FAILED!")
  end
end
