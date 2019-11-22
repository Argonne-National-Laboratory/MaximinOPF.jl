using MaximinOPF
using Ipopt
using PowerModels
using Mosek
using MosekTools


"active power only models"
"abstract type AbstractActivePowerModel <: AbstractPowerModel end"

"variants that target conic solvers"
"abstract type AbstractConicModel <: AbstractPowerModel end"

"for branch flow models"
"abstract type AbstractBFModel <: AbstractPowerModel end"

"for variants of branch flow models that target QP or NLP solvers"
"abstract type AbstractBFQPModel <: AbstractBFModel end"

"for variants of branch flow models that target conic solvers"
"abstract type AbstractBFConicModel <: AbstractBFModel end"

"Branch flow versus bus injection"
"Conic versus general QP/NL" "NOTE: Any algorithm using Dualization.jl must use a ConicModel"

supportedCases = [
	SOCWRPowerModel, # Not Mosek
	SOCWRConicPowerModel, # MoSek
	QCRMPowerModel, # Not Mosek
	SDPWRMPowerModel, # MoSek
	SOCBFPowerModel, # Error constraint_ohms_yt_from()	
	SparseSDPWRMPowerModel # Error variable_voltage()
]
#for i in 1:4
	#Set Default Input
	case = "../data/case9.m"
	#powerfrom = supportedCases[1] #SOCWRPowerModel
	powerfrom = supportedCases[2] #SOCWRConicPowerModel
	println(powerfrom)
	nLineAttacked = 1
	#Create PowerModels Model
	model = MaximinOPF.MaximinOPFModel(case, powerfrom, nLineAttacked)
	# println("Print PowerModels Model")
	# println(model.model)

	#Solve Model with PowerModels Solution Builder
	println("Start Solving")
	#if i < 4
	#	result = optimize_model!(model, with_optimizer(Ipopt.Optimizer))
	#else
		result = optimize_model!(model, with_optimizer(Mosek.Optimizer))
	#end
	println(result)
	println("The optimal value is: ",result["objective"])
#end
