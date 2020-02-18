include("../../MaximinOPF/src/MaximinOPF.jl")
using JuMP, MathOptInterface
using PowerModels
using Ipopt
using Mosek
using MosekTools
#using SCIP
using LinearAlgebra
using Arpack
using Printf
PowerModels.silence()


function solveMaxminECP(pm_data,pm_form,pm_optimizer,io=Base.stdout)
  #println(io,"Formulating and solving the form ",pm_form, " for problem ",pm_data["name"], " with attack budget K=",pm_data["attacker_budget"],".")
  MAX_N_COLS=100
  time_Start = time_ns()
  #feas_pm = MaximinOPF.PF_FeasModel(pm_data, pm_form)
  #println(io,feas_pm.model)

  minmax_pm=MaximinOPF.MinimaxOPFModel(pm_data,pm_form)
  minmax=Dict{String,Any}()
  minmax["branch_ids"] = ids(minmax_pm, minmax_pm.cnw, :branch)
  minmax["model"],minmax["con_dict"]=MaximinOPF.ConvertModelToDualizableForm(minmax_pm.model)
  minmax["psd_con_expr"] = Dict{Int64,Any}()
  minmax["psd_dict"] = minmax["con_dict"]["PSD"]
  for kk=1:length(minmax["psd_dict"])
    minmax_psd_con_kk = constraint_by_name(minmax["model"],minmax["psd_dict"][kk])
    minmax["psd_con_expr"][kk]=constraint_object( minmax_psd_con_kk).func
  end


  maxmin=Dict{String,Any}()
  maxmin["branch_ids"] = ids(minmax_pm, minmax_pm.cnw, :branch)
  maxmin["x_soln"]=Dict{Int64,Float64}()
  for l in maxmin["branch_ids"]
      maxmin["x_soln"][l] = 0
  end
  best_x_soln=copy(maxmin["x_soln"])
  bestLB=-1e20


  feas=Dict{String,Any}()
  feas["branch_ids"] = ids(minmax_pm, minmax_pm.cnw, :branch)
  feas["model"], maprefs = copy_model(minmax["model"])
  feas["arcs"]= ref(minmax_pm,minmax_pm.cnw,:arcs)
  feas["psd_dict"]=minmax["psd_dict"] #shallow copy OK
  for l in feas["branch_ids"]
    delete(feas["model"],constraint_by_name(feas["model"],"x[$l]"))
  end
  JuMP.set_optimizer(feas["model"],pm_optimizer)
  updateFeasObj(feas,maxmin["x_soln"])
  JuMP.optimize!(feas["model"])
  FP_Val=JuMP.objective_value(feas["model"])
  println("FP Objective value: ",FP_Val," with status: ",JuMP.termination_status(feas["model"]))
  if FP_Val > bestLB
      bestLB=FP_Val
      best_x_soln=copy(maxmin["x_soln"])
      println(io,"Updating incumbent solution: ")
      printX(best_x_soln,io)
      println(io,"The incumbent value is now: ",bestLB,".")
  end
  WVal=Dict{Int64,Any}()
  WVal[0]=Dict{Int64,Any}()
  feas["psd_con_expr"] = Dict{Int64,Any}()
  for kk=1:length(feas["psd_dict"])
    feas["psd_con_expr"][kk]=constraint_object(  constraint_by_name(feas["model"],feas["psd_dict"][kk])  ).func
    WVal[0][kk]=trunc.( (JuMP.value.(feas["psd_con_expr"][kk]))[:];digits=10)
    jj=1
    for nn=1:length(WVal[0][kk])
	if nn==(Int)((jj+1)*jj/2)
	  WVal[0][kk][nn]=1
	  jj += 1
	else
	  WVal[0][kk][nn]=0
	end
    end
  end
  #println("WVal[0]", WVal[0])
  



  maxmin["model"] = MaximinOPF.DualizeMinmaxModel(minmax["model"])
  println(io,maxmin["model"])
    ##Convert to add integrality constraints, remove explicit statement of PSD constraints to be replaced with cutting planes
  ###Collect PSD matrix expressions
  psd_con=all_constraints(maxmin["model"], Array{VariableRef,1}, MathOptInterface.PositiveSemidefiniteConeTriangle)
  n_psd_con=length(psd_con)
  maxmin["psd_mat_expr"] = Dict{Int64,Dict{Int64,String}}()
  for kk=1:n_psd_con
    maxmin["psd_mat_expr"][kk] = Dict{Int64,String}()
    psd_con_expr=constraint_object(psd_con[kk]).func
    n_el=length(psd_con_expr)
    for nn=1:n_el
      maxmin["psd_mat_expr"][kk][nn] = name(psd_con_expr[nn])
    end
  end
  
  computeNewConstraintWVal(maxmin, WVal[0];io=io)



  ##### Now remove the explicit PSD constraints
  for kk=1:length(psd_con)
    delete(maxmin["model"],psd_con[kk])
  end
  ##### Fix the x's to be binary
  for l in maxmin["branch_ids"]
      JuMP.fix(variable_by_name(maxmin["model"],"x[$l]_1"),0;force=true)
      JuMP.unfix(variable_by_name(maxmin["model"],"x[$l]_1"))
      JuMP.set_binary(variable_by_name(maxmin["model"],"x[$l]_1"))
  end

  JuMP.set_optimizer(maxmin["model"],pm_optimizer)
  resolveMP(maxmin;io=io)

  updateFeasObj(feas,maxmin["x_soln"])
  JuMP.optimize!(feas["model"])
  FP_Val=JuMP.objective_value(feas["model"])
  println("FP Objective value: ",FP_Val," with status: ",JuMP.termination_status(feas["model"]))

  println("Iteration 0: UB: ",maxmin["UB_Val"]," with status: ",JuMP.termination_status(maxmin["model"]))

  for ii=1:20
    computeNewConstraintEig(maxmin;io=Base.stdout)

    resolveMP(maxmin;io=io)
    println("Iteration $ii: UB: ",maxmin["UB_Val"]," with status: ",JuMP.termination_status(maxmin["model"]))
  
    updateFeasObj(feas,maxmin["x_soln"])
    JuMP.optimize!(feas["model"])
    FP_Val=JuMP.objective_value(feas["model"])
    if FP_Val > bestLB
      bestLB=FP_Val
      best_x_soln=copy(maxmin["x_soln"])
      println(io,"Updating incumbent solution: ")
      printX(best_x_soln,io)
      println(io,"The incumbent value is now: ",bestLB,".")
    end
    println("FP Objective value: ",FP_Val," with status: ",JuMP.termination_status(feas["model"]))

    if maxmin["UB_Val"] - bestLB < 1e-3
	println(io,"Terminating due to closure of gap. Terminal iteration: ",ii)
	break
    end
  end

#=
  for kk in 1:14

    FP_Val=JuMP.objective_value(feas["model"])
    println(io,"Optimal value for FP is: ",FP_Val, " with status ",JuMP.termination_status(feas["model"]))

    computeNewConstraintWVal(maxmin,feas;io=io)
    resolveMP(maxmin,io)
    #computeNewConstraintEig(maxmin;io=io)
    println(io,"Optimal value for MAXMIN is: ",maxmin["UB_Val"], " with status ",JuMP.termination_status(maxmin["model"]))
    print(io,"With x solution: ")
    printX(maxmin["x_soln"],io)

    println(io,"Iteration $kk: UB: ",maxmin["UB_Val"], " and best LB: ",bestLB,".")


  end
=#

end

#=
function tri_sqr_map(nn::Int)
  nn_copy = nn
  ii,jj=1,1
  for jj=1:nn
    if nn_copy-jj > 0
	nn_copy -= jj
	jj += 1
    else
	ii = nn_copy
    end
  end
  return ii,jj
end
function sqr_tri_map(ii::Int,jj::Int)
  nn=(Int)((jj)*(jj-1)/2)
  nn += ii
  return nn
end
=#

function updateFeasObj(feas,x_soln)
  @objective(feas["model"],Min, 
    sum( (1-x_soln[a[1]])*(
		variable_by_name(feas["model"],string("0_up_br1[",a,",",0,"]")) 
		+ variable_by_name(feas["model"],string("0_uq_br1[",a,",",0,"]")) 
		+ variable_by_name(feas["model"],string("0_up_br1[",a,",",1,"]")) 
		+ variable_by_name(feas["model"],string("0_uq_br1[",a,",",1,"]")) 
	 ) 
            + x_soln[a[1]]*(
		variable_by_name(feas["model"],string("0_up_br0[",a,",",0,"]")) 
		+ variable_by_name(feas["model"],string("0_uq_br0[",a,",",0,"]")) 
		+ variable_by_name(feas["model"],string("0_up_br0[",a,",",1,"]")) 
		+ variable_by_name(feas["model"],string("0_uq_br0[",a,",",1,"]")) 
	      ) 
	for a in feas["arcs"])
  )

end

function computeNewConstraintWVal(maxmin, WVal;io=Base.stdout)
  for kk in keys(maxmin["psd_mat_expr"])
      n_el=length(maxmin["psd_mat_expr"][kk])
      PSD_Vec=[variable_by_name(maxmin["model"],maxmin["psd_mat_expr"][kk][nn]) for nn in 1:n_el]
	
      SG = Array{Float64,1}(undef,n_el)
      jj=1
      for nn = 1:n_el
	 if nn==(Int)((jj+1)*(jj)/2)
	     SG[nn] = WVal[kk][nn]
	     jj += 1
	 else
	     SG[nn] = 2*WVal[kk][nn] 
	 end
      end
      @constraint(maxmin["model"], sum(SG[nn]*PSD_Vec[nn] for nn=1:n_el)  >= 0)
  end

end

function computeNewConstraintEig(maxmin;io=Base.stdout)
  for kk in keys(maxmin["psd_mat_expr"])
      n_el=length(maxmin["psd_mat_expr"][kk])
      ncols=(Int)((-1+sqrt(1+8*n_el))/2)
      PSD_Vec=[variable_by_name(maxmin["model"],maxmin["psd_mat_expr"][kk][nn]) for nn in 1:n_el]
      PSD0 = JuMP.value.(PSD_Vec) 
      psd_mat = zeros(ncols,ncols)

      nn=1
      for jj=1:ncols
        for ii=1:jj
	 psd_mat[ii,jj] = PSD0[nn] 
	 psd_mat[jj,ii] = psd_mat[ii,jj]
	 nn += 1
	end
      end
#println(io,"psd_mat_$kk: ",psd_mat)
      E = eigs(psd_mat, which=:SR, nev=2)
      E[2][:,1:2] = round.(E[2][:,1:2],digits=6)
      #psd_mat = round.(psd_mat,digits=6)
#println(io,"eigenval_$kk: ", E[1])
      VVstar=[E[2][i,1]*E[2][j,1] for i=1:ncols,j=1:ncols]
	
#println(io,"VVstar_$kk: ",VVstar)
#println(io,"Frob IP $kk: ", sum(VVstar[i,j]*psd_mat[i,j] for i=1:ncols,j=1:ncols))
      SG = Array{Float64,1}(undef,n_el)
      nn=1
      for j=1:ncols
          for i=1:j
	      if i==j
	        SG[nn] = VVstar[i,j] 
	      else
	        SG[nn] = VVstar[i,j]+VVstar[j,i] 
	      end
	      nn += 1
          end
      end
#println(io,"psd_mat: ",psd_mat)
#println(io,"PSD0: ",PSD0)
#println(io,"PSDVal: ",JuMP.value.(PSD_Vec))

      #@constraint(model, sum(SG[nn]*(PSD_Vec[nn]-PSD0[nn]) for nn=1:n_el) + E[1][1]  >= 0)
#println(io,"Frob IPx $kk: ", JuMP.value( sum(SG[nn]*PSD_Vec[nn] for nn=1:n_el)) )
      cref=@constraint(maxmin["model"], sum(SG[nn]*PSD_Vec[nn] for nn=1:n_el)  >= 0)
#println("cref: ",cref)
  end
end


function resolveMP(maxmin;io=Base.stdout)
  JuMP.optimize!(maxmin["model"])
  maxmin["UB_Val"]=JuMP.objective_value(maxmin["model"])
  for l in keys(maxmin["x_soln"])
    maxmin["x_soln"][l]=min(1,max(JuMP.value(variable_by_name(maxmin["model"],"x[$l]_1")),0))
  end
#println(io,"Checking x_soln: ",maxmin["x_soln"])
end

function printX(x_soln::Dict{Int64,Float64},io=Base.stdout)
  for l in keys(x_soln)
      if x_soln[l] > 0.99
	 print(io," ",l)
      end
  end
  print(io,"\n")
end




testcase = Dict(
	"file" => "data/case30.m", 
 	"name" => "case30K4",  	
 	"attack_budget" => 4,
 	"inactive_indices" => [],
 	"protected_indices" => []
	)

pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry

pm_optimizer=with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=4)

io = open("output.txt", "w")
#io = Base.stdout
pm_form=SparseSDPWRMPowerModel
solveMaxminECP(pm_data,pm_form,pm_optimizer,io)
close(io)
