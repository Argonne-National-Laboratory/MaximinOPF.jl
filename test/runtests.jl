using MaximinOPF
using Ipopt
using PowerModels
using Mosek
using MosekTools

supportedCases = [
	SOCWRPowerModel, # Not Mosek	
	QCRMPowerModel, # Not Mosek
	SOCWRConicPowerModel, # MoSek
	SDPWRMPowerModel#, # MoSek
	#SOCBFPowerModel, # Error constraint_ohms_yt_from()	
	#SparseSDPWRMPowerModel # Error variable_voltage()
]
for i in 1:length(supportedCases)
	#Set Default Input
	case = "../data/case9.m"
	powerfrom = supportedCases[i] #SOCWRPowerModel
	nLineAttacked = i
	#Create PowerModels Model
	model = MaximinOPF.MaximinOPFModel(case, powerfrom, nLineAttacked)
	# println("Print PowerModels Model")
	# println(model.model)

	#Solve Model with PowerModels Solution Builder
	println("Start Solving")
	if i < 3
		result = optimize_model!(model, with_optimizer(Ipopt.Optimizer))
	else
		result = optimize_model!(model, with_optimizer(Mosek.Optimizer))
	end
	println(result["objective"])
end