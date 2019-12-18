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


for i in 1:length(supportedPMOptions)
	# pm_datas = getTestcasesFP()
	pm_datas = getTestcasesMinmax()
	powerfrom = supportedPMOptions[i] #PowerModel Options
	
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
		
		#Create JUMP Model
		maximin_model = MaximinOPF.MaximinOPFModel(pm_data, powerfrom, nLineAttacked)

		#Solve Model with PowerModels Solution Builder
		println("Start Solving")
		if i > 2
			result = @elapsed JuMP.optimize!(maximin_model,with_optimizer(Ipopt.Optimizer))
		else
			result = @elapsed JuMP.optimize!(maximin_model,with_optimizer(Mosek.Optimizer))
		end
		
	end
end

