using MaximinOPF
using Ipopt
using PowerModels
using Mosek
using MosekTools

supportedCases = [	
	SOCWRConicPowerModel, # MoSek
#	SDPWRMPowerModel, # MoSek
#	SOCWRPowerModel, # Not Mosek	
#	QCRMPowerModel # Not Mosek
	#SOCBFPowerModel, # Error constraint_ohms_yt_from()	
	#SparseSDPWRMPowerModel # Error variable_voltage()
]
# "Case9 feas problem with no lines cut should have value 0.0"
# "Case9 feas problem with lines 1,4,7 cut should have value 3.451901668361513"
# "Case9 feas problem with lines 1--9 cut should have value 4.6"
# "Case30 feas problem with no lines cut should have value 0.0"
# "Case30 feas problem with lines 8,9,10,40 cut should have value 0.937"
# "Case30 feas problem with all lines cut should have value 2.57327375"
# "Case57 feas problem with no lines cut should have value 0.0"
# "Case57 feas problem with lines 41,80 cut should have value "
for i in 1:length(supportedCases)
	#Set Default Input
	#case = "../data/case9.m"
	#case = "../data/case30.m"
	case = "../data/case57.m"
	pm_data = PowerModels.parse_file(case)
	print("Cutting lines: ")
        for (k,line) in pm_data["branch"]
	  #if line["index"]==1 || line["index"]==4 || line["index"]==7
	  #if line["index"]==8 || line["index"]==9 || line["index"]==10 || line["index"]==40
	  #if line["index"]==8 || line["index"]==9 || line["index"]==10 || line["index"]==40
	  if line["index"]==41 || line["index"]==80
	    line["br_status"]=0
	    print(" ",line["index"])
	  end
        end
	print("\n")
	powerfrom = supportedCases[i] #SOCWRPowerModel
	nLineAttacked = i
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
	println(result["objective"])
end
