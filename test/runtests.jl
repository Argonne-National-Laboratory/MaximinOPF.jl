using MaximinOPF
using Ipopt
using PowerModels

#Set Default Input
case = "../data/case9.m"
powerfrom = ACPPowerModel #SOCWRPowerModel
Maxline = 1
#Create PowerModels Model
model = MaximinOPF.MaximinOPFModel(case, ACPPowerModel, Maxline)
println(model.model)

#Solve Model with PowerModels Solution Builder
result = optimize_model!(model, with_optimizer(Ipopt.Optimizer))
println(result)