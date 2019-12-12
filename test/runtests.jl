using MaximinOPF
using Ipopt
using PowerModels
using Mosek
using MosekTools
using JuMP

include("testcases.jl")

supportedPMOptions = [	
	SOCWRConicPowerModel, # MoSek
#	SDPWRMPowerModel, # MoSek
#	SOCWRPowerModel, # Not Mosek	
#	QCRMPowerModel # Not Mosek
	#SOCBFPowerModel, # Error constraint_ohms_yt_from()	
	#SparseSDPWRMPowerModel # Error variable_voltage()
]

function evaluateMinMax(expectedvalue, model, tol)
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

function evaluateResult(expectedvalue, solvedvalue, tol)
	if expectedvalue - solvedvalue < tol
		return 1
	else
		return -1
	end
end


testresults = []
for i in 1:length(supportedPMOptions)
	# pm_datas = getTestcasesFP()
	pm_datas = getTestcasesMinmax()
	powerfrom = supportedPMOptions[i] #PowerModel Options
	#println(length(pm_datas))
	for j in 1:length(pm_datas)
		casename = pm_datas[j]["name"]
		expect = pm_datas[j]["expectedvalue"]		
		pm_data = pm_datas[j]["pm_data"]
		nLineAttacked = pm_datas[j]["K"]
		protected_indices = pm_datas[j]["protected_indices"]
		lineindexs = pm_datas[j]["inactive_indices"]
		
		pm_data["attacker_budget"]=nLineAttacked ###Adding another key and entry
		pm_data["inactive_branches"]=lineindexs ###Adding another key and entry
		pm_data["protected_branches"]=protected_indices ###Adding another key and entry
		#Create PowerModels Model
		model = MaximinOPF.MaximinOPFModel(pm_data, powerfrom, nLineAttacked)
		#println(model)
		#println("Print PowerModels Model")
		#println(model.model)

		#Solve Model with PowerModels Solution Builder
		println("Start Solving")
		if i > 2
			result = optimize_model!(model, with_optimizer(Ipopt.Optimizer))
		else
			result = optimize_model!(model, with_optimizer(Mosek.Optimizer))
		end
		#println(result["objective"])
		
		evaluateMinMax(expect, model, 0.001)
		testresult = evaluateResult(expect, result["objective"], 0.001)
		
		push!(testresults, 
			Dict(
				"testresult" => testresult,
				"casename" => casename,
				"protected_branches" => protected_indices,
				"inactive_branches" => lineindexs,
				"expectedvalue" => expect,
				"solvedvalue" => result["objective"]
				))
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
