using MaximinOPF
using Ipopt
using PowerModels
using Mosek
using MosekTools

include("testcases.jl")

supportedPMOptions = [	
	SOCWRConicPowerModel, # MoSek
#	SDPWRMPowerModel, # MoSek
#	SOCWRPowerModel, # Not Mosek	
#	QCRMPowerModel # Not Mosek
	#SOCBFPowerModel, # Error constraint_ohms_yt_from()	
	#SparseSDPWRMPowerModel # Error variable_voltage()
]


function evaluation(expectedvalue, solvedvalue, tol)
	if expectedvalue - solvedvalue < tol
		return 1
	else
		return -1
	end
end


testresults = []
for i in 1:length(supportedPMOptions)
	#pm_datas = getTestcasesFP()
	pm_datas = getTestcasesRMinmax()
	powerfrom = supportedPMOptions[i] #PowerModel Options
	print(length(pm_datas))
	for j in 1:length(pm_datas)
		pm_data = pm_datas[j]["pm_data"]
		nLineAttacked = pm_datas[j]["K"]
		
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
		
		expect = pm_datas[j]["expectedvalue"]
		lineindexs = pm_datas[j]["cutindex"]
		casename = pm_datas[j]["name"]
		testresult = evaluation(expect, result["objective"], 0.001)
		push!(testresults, 
			Dict(
				"testresult" => testresult,
				"casename" => casename,
				"cutindex" => lineindexs,
				"expectedvalue" => expect,
				"solvedvalue" => result["objective"]
				))
	end
end

println(testresults)
