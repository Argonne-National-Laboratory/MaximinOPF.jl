include("../../MaximinOPF/src/MaximinOPF.jl")
using JuMP
using PowerModels
using Ipopt
using Mosek
using MosekTools
#using SCIP
using LinearAlgebra
using Arpack
using Printf
PowerModels.silence()

testcase = Dict(
	"file" => "data/case9.m", 
	"PMOption" => SparseSDPWRMPowerModel,
 	"name" => "case9K3",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => []
	)

pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry

io = open("output.txt", "w")
# nonconvex AC forms
nonconvex_ac=[ACPPowerModel, ACRPowerModel, ACTPowerModel]
for pm_form in nonconvex_ac
  println("Formulating and solving the form ",pm_form)
  model,pm=MaximinOPF.SolveMinmax(pm_data,pm_form, with_optimizer(Ipopt.Optimizer))
  println(io,"Optimal value using powerform ", pm_form, " is: ",JuMP.objective_value(pm.model), " with status ",JuMP.termination_status(pm.model))
end

# linear approximations
linear_approx=[DCPPowerModel, DCMPPowerModel, NFAPowerModel]
for pm_form in linear_approx
  println("Formulating and solving the form ",pm_form)
  model,pm=MaximinOPF.SolveMinmax(pm_data,pm_form, with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0))
  println(io,"Optimal value using powerform ", pm_form, " is: ",JuMP.objective_value(pm.model), " with status ",JuMP.termination_status(pm.model))
end

# quadratic approximations
quadratic_approx=[DCPLLPowerModel, LPACCPowerModel]
for pm_form in quadratic_approx
  println("Formulating and solving the form ",pm_form)
  model,pm=MaximinOPF.SolveMinmax(pm_data,pm_form, with_optimizer(Ipopt.Optimizer))
  println(io,"Optimal value using powerform ", pm_form, " is: ",JuMP.objective_value(pm.model), " with status ",JuMP.termination_status(pm.model))
end
# quadratic relaxations
quadratic_relax=[SOCWRPowerModel, SOCBFPowerModel, QCRMPowerModel, QCLSPowerModel]
for pm_form in quadratic_relax
  println("Formulating and solving the form ",pm_form)
  model,pm=MaximinOPF.SolveMinmax(pm_data,pm_form, with_optimizer(Ipopt.Optimizer))
  println(io,"Optimal value using powerform ", pm_form, " is: ",JuMP.objective_value(pm.model), " with status ",JuMP.termination_status(pm.model))
end
quad_conic_relax=[SOCWRConicPowerModel, SOCBFConicPowerModel]
for pm_form in quad_conic_relax
  println("Formulating and solving the form ",pm_form)
  model,pm=MaximinOPF.SolveMinmax(pm_data,pm_form, with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0))
  println(io,"Optimal value using powerform ", pm_form, " is: ",JuMP.objective_value(pm.model), " with status ",JuMP.termination_status(pm.model))
end
# sdp relaxations
sdp_relax=[SDPWRMPowerModel, SparseSDPWRMPowerModel]
for pm_form in sdp_relax
  println("Formulating and solving the form ",pm_form)
  model,pm=MaximinOPF.SolveMinmax(pm_data,pm_form, with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0))
  println(io,"Optimal value using powerform ", pm_form, " is: ",JuMP.objective_value(pm.model), " with status ",JuMP.termination_status(pm.model))
end
close(io)

#=
W=Dict{Int,Float64}()
WR=Dict{Int,Float64}()
WI=Dict{Int,Float64}()
for i in ids(pm, :bus)
  W[i] = JuMP.value(var(pm, pm.cnw, :w)[i])
end
for l in ids(pm, :branch)
  branch = ref(pm, pm.cnw, :branch, l)
  f_bus = branch["f_bus"]
  t_bus = branch["t_bus"]
  WR[l] = JuMP.value(var(pm, pm.cnw, :wr)[(f_bus, t_bus)])
  WI[l] = JuMP.value(var(pm, pm.cnw, :wi)[(f_bus, t_bus)])
end
println("W: ", W)
println("WR: ",WR)
println("WI: ",WI)


fr_bus=Dict{Int,Int}()
to_bus=Dict{Int,Int}()
for l in ids(pm,:branch)
  branch = ref(pm, pm.cnw, :branch, l)
  f_bus = branch["f_bus"]
  t_bus = branch["t_bus"]
  fr_bus[l]=f_bus
  to_bus[l]=t_bus
end
println(fr_bus)
println(to_bus)

println(pm.data["bus"])
=#


#=
m_volt=JuMP.Model()
@variable(m_volt, vr[i in ids(pm, :bus)])
@variable(m_volt, vi[i in ids(pm, :bus)])
JuMP.@constraint(m_volt,[i in ids(pm,:ref_buses)], vi[i]==0)

@NLobjective(m_volt, Min,
	sum( (vr[i]*vr[i]+vi[i]*vi[i] - W[i])^2 for i in ids(pm,:bus)) 
	+ sum( (vr[to_bus[l]]*vr[fr_bus[l]]+vi[fr_bus[l]]*vi[to_bus[l]] - WR[l])^2   for l in ids(pm,:branch))
	+ sum( (vr[to_bus[l]]*vi[fr_bus[l]]-vr[fr_bus[l]]*vi[to_bus[l]] - WI[l])^2   for l in ids(pm,:branch))
)
JuMP.set_optimizer(m_volt,with_optimizer(Ipopt.Optimizer))
JuMP.optimize!(m_volt)
status=JuMP.termination_status(m_volt)
println("Optimal value with status ",status," for the bus voltage finding subproblem is: ",JuMP.objective_value(m_volt))
=#

#=

m_sdp=JuMP.Model()
nbuses=length(ids(pm,:bus))
@variable(m_sdp,wmatr[i=1:nbuses,j=1:nbuses])
@variable(m_sdp,wmati[i=1:nbuses,j=1:nbuses])
@constraint(m_sdp, CDiagR[i in ids(pm,:bus)], wmatr[i,i] == W[i])
@constraint(m_sdp, CDiagI[i in ids(pm,:bus)], wmati[i,i] == 0)
@constraint(m_sdp, CR[l in ids(pm,:branch)], wmatr[fr_bus[l],to_bus[l]] == WR[l])
@constraint(m_sdp, CI1[l in ids(pm,:branch)], wmati[fr_bus[l],to_bus[l]] == WI[l])
@constraint(m_sdp, CI2[l in ids(pm,:branch)], wmati[to_bus[l],fr_bus[l]] == -WI[l])
@constraint(m_sdp, [wmatr wmati; -wmati wmatr] in JuMP.PSDCone())
@objective(m_sdp, Min, sum(wmatr[i,i]^2 for i in ids(pm,:bus)) )
JuMP.set_optimizer(m_sdp,with_optimizer(Mosek.Optimizer))
JuMP.optimize!(m_sdp)
status=JuMP.termination_status(m_sdp)
#println("Optimal solution with status ",status," for the bus voltage finding subproblem is: ")
WVal=zeros(2*nbuses,2*nbuses)
for i in ids(pm,:bus)
  for j in ids(pm,:bus)
    WVal[i,j] = JuMP.value(wmatr[i,j])
    WVal[j,i] = WVal[i,j]
    WVal[nbuses+i,nbuses+j] = WVal[i,j]
    WVal[nbuses+j,nbuses+i] = WVal[i,j]
    #print(@sprintf(" %.2f",JuMP.value(wmatr[i,j])))
  end
  #print("\n")
end
for i in ids(pm,:bus)
  for j in ids(pm,:bus)
    WVal[nbuses+i,j] = JuMP.value(wmati[i,j])
    WVal[nbuses+j,i] = -WVal[nbuses+i,j]
    WVal[i,nbuses+j] = WVal[nbuses+j,i]
    WVal[j,nbuses+i] = WVal[nbuses+i,j]
    #print(@sprintf(" %.2f",JuMP.value(wmati[i,j])))
  end
  #print("\n")
end
#println(round.(WVal; digits=3))

neigs=10
E = eigs(WVal, which=:LM, nev=neigs)
#println("Eigenvalues: ",E[1])
#println("Eigenvectors: ",E[2])
println("A candidate bus voltage solution is: ")
for i in ids(pm,:bus)
  vR[i]= sum(sqrt(E[1][k])*E[2][i,k] for k=1:1)
  vI[i]= sum(sqrt(E[1][k])*E[2][nbuses+i,k] for k=1:1)
end
println(vR)
println(vI)
println("The voltage magnitudes are: ")
for i in ids(pm,:bus)
  print(" ",sqrt(vR[i]^2+vI[i]^2))
end
println()


bus_idx=ids(pm,:bus)
vR=Dict{Int,Float64}()
vI=Dict{Int,Float64}()
testcase["PMOption"]=QCRMPowerModel
for kk=1:neigs
  if E[1][kk] > 1e-4
    for i in bus_idx
     #vR[i]=2*rand()-1.0
     #vI[i]=2*rand()-1.0
      vR[i]=sqrt(E[1][kk])*E[2][i,kk]
      vI[i]=sqrt(E[1][kk])*E[2][nbuses+i,kk]
    end
    pm=solveNodeMinmax(pm_data,testcase["PMOption"], Ipopt.Optimizer) # ;vRInit=vR,vIInit=vI)
    ACRVal=JuMP.objective_value(pm.model)
    if ACRVal-SDPVal < 1e-3
      println("At iteration ",kk,", gap with the SDPVal has been closed: ACR global opt val is: ",ACRVal)
      break
    end
  end
end
=#
