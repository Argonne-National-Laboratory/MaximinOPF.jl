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
  println(io,"Minmax: ",minmax["model"])


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
    #println("WVal[0][kk]", WVal[0][kk])
  end
  



  maxmin["model"] = MaximinOPF.DualizeMinmaxModel(minmax["model"])
  println(io,maxmin["model"])
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
  
#=
  int_relaxed_model = copy_model(maxmin["model"])
  JuMP.set_optimizer(int_relaxed_model,pm_optimizer)
  JuMP.optimize!(int_relaxed_model)
  int_relaxed_UB=JuMP.objective_value(int_relaxed_model)

  JuMP.set_optimizer(maxmin["model"],pm_optimizer)
  resolveMP(maxmin;io=io)
  int_relaxed_UB = maxmin["UB_Val"]
  println("Integer-relaxed UB: ",int_relaxed_UB," with status: ",JuMP.termination_status(maxmin["model"]))

  @variable(maxmin["model"], obj_var <= int_relaxed_UB)
  obj_func = JuMP.objective_function(maxmin["model"])
  @constraint(maxmin["model"], obj_con, obj_func - obj_var >= 0)
  @objective(maxmin["model"], Max, obj_var)
=#

  
  ##Convert to add integrality constraints, remove explicit statement of PSD constraints to be replaced with cutting planes
  ##### Now remove the explicit PSD constraints
  n_psd_con = length(psd_con)
  for kk=1:length(psd_con)
    delete(maxmin["model"],psd_con[kk])
  end
  computeNewConstraintWVal(maxmin, WVal[0];io=io)

  ##### Fix the x's to be binary
  for l in maxmin["branch_ids"]
      #JuMP.fix(variable_by_name(maxmin["model"],"x[$l]_1"),0;force=true)
      #JuMP.unfix(variable_by_name(maxmin["model"],"x[$l]_1"))
      if has_lower_bound(variable_by_name(maxmin["model"],"x[$l]_1"))
        JuMP.delete_lower_bound(variable_by_name(maxmin["model"],"x[$l]_1"))
      end
      if has_upper_bound(variable_by_name(maxmin["model"],"x[$l]_1"))
        JuMP.delete_upper_bound(variable_by_name(maxmin["model"],"x[$l]_1"))
      end
      JuMP.set_binary(variable_by_name(maxmin["model"],"x[$l]_1"))
  end

  JuMP.set_optimizer(maxmin["model"],pm_optimizer)
  resolveMP(maxmin; io=io)
  #computeNewConstraintEig(maxmin;io=io)


  print("Iteration 0: UB: ",maxmin["UB_Val"]," with status: ",JuMP.termination_status(maxmin["model"])," with x solution: ")
    printX(maxmin["x_soln"])

  for ii=1:40
    println(io,"Iteration $ii")
    updateFeasObj(feas,maxmin["x_soln"])
    JuMP.optimize!(feas["model"])
    FP_Val=JuMP.objective_value(feas["model"])
    if FP_Val > bestLB
      bestLB=FP_Val
      best_x_soln=copy(maxmin["x_soln"])
      println(io,"Updating incumbent solution: ")
      printX(best_x_soln,io)
      println(io,"The incumbent value is now: ",bestLB,".")
      if io!=Base.stdout
        println("Updating incumbent solution: ")
        printX(best_x_soln)
        println("The incumbent value is now: ",bestLB,".")
      end
    end
    WVal[ii]=Dict{Int64,Any}()
    for kk=1:length(feas["psd_dict"])
      WVal[ii][kk]=trunc.( (JuMP.value.(feas["psd_con_expr"][kk]))[:];digits=10)
    end

    computeNewConstraintWVal(maxmin, WVal[ii];io=io)
    resolveMP(maxmin;io=io)
    #computeNewConstraintEig(maxmin;io=io)
    println(io,"FP Objective value: ",FP_Val," with status: ",JuMP.termination_status(feas["model"])," and bestLB is: ",bestLB)
    if io!=Base.stdout
      println("FP Objective value: ",FP_Val," with status: ",JuMP.termination_status(feas["model"])," and bestLB is: ",bestLB)
    end
    print(io,"Iteration $ii: UB: ",maxmin["UB_Val"]," status: ",JuMP.termination_status(maxmin["model"])," with x solution: ")
    if io!=Base.stdout
      print("Iteration $ii: UB: ",maxmin["UB_Val"]," status: ",JuMP.termination_status(maxmin["model"])," with x solution: ")
    end
    printX(maxmin["x_soln"],io)
    if io!=Base.stdout
      printX(maxmin["x_soln"])
    end

    if maxmin["UB_Val"] - bestLB < 1e-3
	println(io,"Terminating due to closure of gap. Terminal iteration: ",ii)
        if io!=Base.stdout
          printX(maxmin["x_soln"])
        end
	break
    end
  end
  time_End = (time_ns()-time_Start)/1e9
  println("Finishing after ", time_End," seconds.")
  print("Best value: ", bestLB, " with solution: ")
  printX(best_x_soln)
end

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
  #println("New feas objective: ", JuMP.objective_function(feas["model"]) )
end

function computeNewConstraintWVal(maxmin, WVal;io=Base.stdout)
  n_cuts = 0
  for kk in keys(maxmin["psd_mat_expr"])
    n_psd_vars=length(maxmin["psd_mat_expr"][kk])
    PSD_Vec=[variable_by_name(maxmin["model"],maxmin["psd_mat_expr"][kk][nn]) for nn in 1:n_psd_vars]
    jj=1
    SG = Array{Float64,1}(undef,n_psd_vars)
    for nn=1:n_psd_vars
      if nn == (Int)(((jj+1)*jj)/2)
        #JuMP.set_lower_bound(variable_by_name(maxmin["model"],maxmin["psd_mat_expr"][kk][nn]),0)
	SG[nn] = WVal[kk][nn]
	jj += 1
      else
	SG[nn] = 2*WVal[kk][nn]
      end
    end
    cref=@constraint(maxmin["model"], sum(trunc(SG[nn];digits=6)*PSD_Vec[nn] for nn=1:n_psd_vars)  >= -1e-6)
    n_cuts += 1
    if JuMP.termination_status(maxmin["model"]) != OPTIMIZE_NOT_CALLED
      PSD0 = JuMP.value.(PSD_Vec)
      IP_SG_PSD0 = triFrobIP(SG,PSD0)
      println(io,"<SG,PSD0>: ",IP_SG_PSD0)
    end
  end
  println(io,"WVal: Adding $n_cuts cuts.")
end

###Inputs are symmetric matrices in triangular form
function triFrobIP(tri_mat1::Array{Float64,1},tri_mat2::Array{Float64,1})
   n_el=length(tri_mat1)
   @assert(n_el==length(tri_mat2))
   result = 0.0
   jj=1
   
   for nn=1:n_el
      if nn == (Int)(((jj+1)*jj)/2)
	result += tri_mat1[nn]*tri_mat2[nn]
	jj += 1
      else
	result += 2*tri_mat1[nn]*tri_mat2[nn]
      end
    end
    return result
end
###Input is symmetric matrix in triangular form
function computeEigenDecomp(tri_mat::Array{Float64,1})
   n_el=length(tri_mat)
   ncols=(Int)((-1+sqrt(1+8*n_el))/2)
   psd_mat = zeros(ncols,ncols)
   nn=1
   sq_tri_map=Dict{Tuple{Int,Int},Int}()
   tri_sq_map=Dict{Int,Tuple{Int,Int}}()
   for jj=1:ncols
      for ii=1:jj
	 tri_sq_map[nn]=(ii,jj)
	 sq_tri_map[(ii,jj)]=nn
	 psd_mat[ii,jj] = tri_mat[nn] 
	 psd_mat[jj,ii] = psd_mat[ii,jj]
	 nn += 1
      end
   end
   E = eigs(psd_mat, which=:SR, nev=(ncols-1))
   return E, sq_tri_map, tri_sq_map
end

function computeNewConstraintEig(maxmin;io=Base.stdout)
  n_cuts = 0
  for kk in keys(maxmin["psd_mat_expr"])
      println(io,"$kk th PSD constraint")
      n_el = length(maxmin["psd_mat_expr"][kk])
      PSD_Vec=[variable_by_name(maxmin["model"],maxmin["psd_mat_expr"][kk][nn]) for nn in 1:n_el]
      PSD0 = JuMP.value.(PSD_Vec) 
#println(io,"PSD0[$kk}: ",PSD0)
      E, sq_tri_map, tri_sq_map = computeEigenDecomp(PSD0)
      n_eigs = length(E[1])
println(io,"eigenval_$kk: ", E[1])
	
     for mm in 1:n_eigs # filter(mmm->(E[1][mmm]<0), 1:n_eigs)
        SG = Array{Float64,1}(undef,n_el)
        for nn=1:n_el
	  (ii,jj)=tri_sq_map[nn]
	  if ii==jj
	     SG[nn] = E[2][ii,mm]*E[2][jj,mm]  
	  else
	     SG[nn] = E[2][ii,mm]*E[2][jj,mm] + E[2][jj,mm]*E[2][ii,mm]   
	  end
        end
        cref=@constraint(maxmin["model"], sum(trunc(SG[nn];digits=6)*PSD_Vec[nn] for nn=1:n_el)  >= -1e-2)
        n_cuts += 1
     end
  end
  println("Eigs: Adding $n_cuts cuts.")
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
