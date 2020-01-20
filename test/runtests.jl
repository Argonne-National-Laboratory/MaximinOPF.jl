using MaximinOPF
using Ipopt
using PowerModels
using Mosek
using MosekTools
using JuMP

include("testcases.jl")
PowerModels.silence()
supportedPMOptions = [	
	SOCWRConicPowerModel, # MoSek
#	SDPWRMPowerModel, # MoSek
#	SparseSDPWRMPowerModel, # MoSek, requires minor editing of PowerModels wrm.jl code as of 5 Jan 2020, open issue with PowerModels is in progress.
#	SOCWRPowerModel, # Not Mosek	
#	QCRMPowerModel # Not Mosek
#	SOCBFPowerModel, # Error constraint_ohms_yt_from()	
]

testresults = []
for i in 1:length(supportedPMOptions)
	powerform = supportedPMOptions[i] #PowerModel Options	
	for j in 1:length(testcases)
		pm_data = PowerModels.parse_file( testcases[j]["file"] )
		pm_data["name"]=testcases[j]["name"]
		pm_data["attacker_budget"] = testcases[j]["attack_budget"] ###Adding another key and entry
		pm_data["inactive_branches"] = testcases[j]["inactive_indices"] ###Adding another key and entry
		pm_data["protected_branches"] = testcases[j]["protected_indices"] ###Adding another key and entry
		
		#Create JUMP Model
		maxmin_model = MaximinOPF.MaximinOPFModel(pm_data, powerform)	

		println("Start Solving")
		if i > 3 
			result = @elapsed JuMP.optimize!(maxmin_model,with_optimizer(Ipopt.Optimizer))
		else
			result = @elapsed JuMP.optimize!(maxmin_model,with_optimizer(Mosek.Optimizer))			
		end

		#Print Result		
		status=JuMP.termination_status(maxmin_model)
		println("Time taken to solve is: ", result, " with status ",status,".")		
		


	end
end

