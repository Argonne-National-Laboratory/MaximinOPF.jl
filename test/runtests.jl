using MaximinOPF
using Ipopt
using PowerModels
using Mosek
using MosekTools
using JuMP
#using SDPA

include("testcases.jl")

PowerModels.silence()

supportedPMOptions = [	
#	SOCWRConicPowerModel, # MoSek
#	SDPWRMPowerModel, # MoSek
	SparseSDPWRMPowerModel, # MoSek, requires minor editing of PowerModels wrm.jl code as of 5 Jan 2020, open issue with PowerModels is in progress.
#	SOCWRPowerModel, # Not Mosek	
#	QCRMPowerModel # Not Mosek
#	SOCBFPowerModel, # Error constraint_ohms_yt_from()	
]

function evaluateMinmaxModel(expected_value, pm, tol, io)
    println(io,"Protected branches: ",pm.data["protected_branches"])
    println(io,"Inactive branches: ",pm.data["inactive_branches"])
    println(io,"All other branches: ",setdiff(ids(pm, :branch),pm.data["protected_branches"]))
    println(io,"Printing dual values x:")
    for l in pm.data["undecided_branches"]
      x_val=abs(JuMP.dual(con(pm, pm.cnw, pm.ccnd)[:x][l])) 
      if x_val > tol
        println(io,"(",l,",",x_val,")")
      end
    end

    lkx = []
    psival = pm.data["attacker_budget"]*JuMP.value(var(pm, pm.cnw, pm.ccnd)[:u_K])
    psival += sum( JuMP.value(var(pm, pm.cnw, pm.ccnd)[:u_ord_aux][l]) for l in pm.data["undecided_branches"] ) 
    println("uK: ",JuMP.value(var(pm, pm.cnw, pm.ccnd)[:u_K])," psival: ",psival)
#=
    for l in ids(pm, :branch)
      branch = ref(pm, pm.cnw, :branch, l)
      f_bus = branch["f_bus"]
      t_bus = branch["t_bus"]
      f_idx = (l, f_bus, t_bus)
      t_idx = (l, t_bus, f_bus)
      upf0 = JuMP.value(var(pm, pm.cnw, pm.ccnd,:up_br)[f_idx,0])
      upt0 = JuMP.value(var(pm, pm.cnw, pm.ccnd,:up_br)[t_idx,0])
      uqf0 = JuMP.value(var(pm, pm.cnw, pm.ccnd,:uq_br)[f_idx,0])
      uqt0 = JuMP.value(var(pm, pm.cnw, pm.ccnd,:uq_br)[t_idx,0])
      upf1 = JuMP.value(var(pm, pm.cnw, pm.ccnd,:up_br)[f_idx,1])
      upt1 = JuMP.value(var(pm, pm.cnw, pm.ccnd,:up_br)[t_idx,1])
      uqf1 = JuMP.value(var(pm, pm.cnw, pm.ccnd,:uq_br)[f_idx,1])
      uqt1 = JuMP.value(var(pm, pm.cnw, pm.ccnd,:uq_br)[t_idx,1])

      upf1x = min(abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr_disc_ub][l])), abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr_disc_lb][l]))) 
      upt1x = min(abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to_disc_ub][l])), abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to_disc_lb][l]))) 
      uqf1x = min(abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr_disc_ub][l])), abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr_disc_lb][l]))) 
      uqt1x = min(abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to_disc_ub][l])), abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to_disc_lb][l]))) 

      upf0x = min(abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr_ub][l])), abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr_lb][l]))) 
      upt0x = min(abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to_ub][l])), abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to_lb][l]))) 
      uqf0x = min(abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr_ub][l])), abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr_lb][l]))) 
      uqt0x = min(abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to_ub][l])), abs(JuMP.value(con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to_lb][l]))) 

      upf0 -= upf0x
      upt0 -= upt0x
      uqf0 -= uqf0x
      uqt0 -= uqt0x
      upf1 -= upf1x
      upt1 -= upt1x
      uqf1 -= uqf1x
      uqt1 -= uqt1x

      push!(lkx, (l,upf0+upt0+uqf0+uqt0,-upf1-upt1-uqf1-uqt1,upf0+upt0+uqf0+uqt0-upf1-upt1-uqf1-uqt1))	        
      println(io,"x subgrad value: ",(l,upf0+upt0+uqf0+uqt0-upf1-upt1-uqf1-uqt1))	        
      println(io,"Comparing u: ",(l,(upf0,upt0,uqf0,uqt0),(upf1,upt1,uqf1,uqt1,)))	        
    end
=#
    
    
    discrep = abs(JuMP.objective_value(pm.model)-expected_value)
    if discrep < tol
	println(io,"The test for ",pm.data["name"]," PASSES: Obj value is: ",JuMP.objective_value(pm.model),", and the expected value was: ",expected_value,".")
    else
	println(io,"The test for ",pm.data["name"]," FAILED: Obj value is: ",JuMP.objective_value(pm.model),", but the expected value was: ",expected_value,".")
    end
    
end

testresults = []
for i in 1:length(supportedPMOptions)
	powerform = supportedPMOptions[i] #PowerModel Options
	ncases = length(casesConic)
	#for j in 1:1
	for j in 1:(ncases)
	#for j = ncases:ncases
		pm_data = PowerModels.parse_file( casesConic[j]["file"] )
		pm_data["name"]=casesConic[j]["name"]
		pm_data["attacker_budget"] = casesConic[j]["attack_budget"] ###Adding another key and entry
		pm_data["inactive_branches"] = casesConic[j]["inactive_indices"] ###Adding another key and entry
		pm_data["protected_branches"] = casesConic[j]["protected_indices"] ###Adding another key and entry


		f_name_base="output_"
		f_name = string(f_name_base,casesConic[j]["name"],"_", pm_data["attacker_budget"],"_",j,".txt")
		io=open(f_name, "w")
		
		#Create JUMP Model
		#pf_model_pm = MaximinOPF.PF_FeasModel(pm_data, powerform)
		pf_model_pm = MaximinOPF.MinimaxOPFModel(pm_data, powerform)
		pm_data["undecided_branches"]= filter(l->!(l in pm_data["protected_branches"] || l in pm_data["inactive_branches"]), ids(pf_model_pm,pf_model_pm.cnw,:branch)) 
			###Adding another key and entry
		pf_minmax_model = pf_model_pm.model
		#pf_maxmin_model = MaximinOPF.DualizeModel(pf_model_pm)
		
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
		opt_expected = "SDP_Minmax"
		println("Start Solving")
		if i > 3 
			set_optimizer(opt_model,Ipopt.Optimizer)
			#result = @elapsed JuMP.optimize!(opt_model,with_optimizer(Ipopt.Optimizer))
		else
			set_optimizer(opt_model,Mosek.Optimizer)
			#set_optimizer(opt_model,SDPA.Optimizer)
			#result = @elapsed JuMP.optimize!(opt_model,with_optimizer(Mosek.Optimizer))
			#result = @elapsed JuMP.optimize!(opt_model,with_optimizer(SDPA.Optimizer))
		end
		result = @elapsed JuMP.optimize!(opt_model)
		
		#Print Result
		#evaluateModel(expect, pf_model_pm, 0.001)
		status=JuMP.termination_status(opt_model)
		println(io,"Time taken to solve is: ", result, " with status ",status,".")
                evaluateMinmaxModel(casesConic[j]["expected_values"][opt_expected], pf_model_pm, 1e-3,io)
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
