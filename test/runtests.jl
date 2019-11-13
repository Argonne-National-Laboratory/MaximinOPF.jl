using MaximinOPF
using Ipopt
using PowerModels
using Mosek
using MosekTools

supportedCases = [
	SOCWRPowerModel,
	SOCWRConicPowerModel,
	QCRMPowerModel,
	SOCBFPowerModel,
	SDPWRMPowerModel,
	SparseSDPWRMPowerModel
]
for c in supportedCases
	#Set Default Input
	case = "../data/case9.m"
	powerfrom = c #SOCWRPowerModel
	nLineAttacked = 1
	#Create PowerModels Model
	model = MaximinOPF.MaximinOPFModel(case, powerfrom, nLineAttacked)
	# println("Print PowerModels Model")
	# println(model.model)

	#Solve Model with PowerModels Solution Builder
	println("Start Solving")
	result = optimize_model!(model, with_optimizer(Mosek.Optimizer))
	println(result)
end