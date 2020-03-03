include("../src/MaximinOPF.jl")
using JuMP, MathOptInterface
using PowerModels
using Ipopt
using Mosek
using MosekTools
using CPLEX
#using SCIP
using LinearAlgebra
using SparseArrays
using Arpack
using Printf
PowerModels.silence()

function convertSOCtoPSD(model::JuMP.Model)
    #BRIDGE SOC CONSTRAINTS
    model_rsoc_moi = MOI.Utilities.Model{Float64}()
    rsoc_bridged_model = MOI.Bridges.Constraint.SOCtoPSD{Float64}(model_rsoc_moi)
    MOI.copy_to(rsoc_bridged_model,backend(model))

    model_psd_moi = MOI.Utilities.Model{Float64}()
    psd_bridged_model = MOI.Bridges.Constraint.RSOCtoPSD{Float64}(model_psd_moi)
    MOI.copy_to(psd_bridged_model,model_rsoc_moi)
    model_psd = JuMP.Model()
    MOI.copy_to(backend(model_psd),model_psd_moi)
    return model_psd
end

function solveMaxminECP(pm_data,pm_form,pm_optimizer,io=Base.stdout)
    #println(io,"Formulating and solving the form ",pm_form, " for problem ",pm_data["name"], " with attack budget K=",pm_data["attacker_budget"],".")
    MAX_N_COLS=100
    MAX_N_ITER=40
    using_eigs = true
    time_Start = time_ns()
    feas_pm = MaximinOPF.PF_FeasModel(pm_data, pm_form)
    #println(io,feas_pm.model)

    maxmin=Dict{String,Any}()
    maxmin["branch_ids"] = pm_data["undecided_branches"]
    println(maxmin["branch_ids"])
    maxmin["x_soln"]=Dict{Int64,Float64}()
    maxmin["x_soln_01"]=spzeros(Int64,maximum(maxmin["branch_ids"]))
    for l in maxmin["branch_ids"]
        maxmin["x_soln"][l] = 0
    end
    best_x_soln=copy(maxmin["x_soln"])
    bestLB=-1e20

    base_maxmin = MaximinOPF.MaximinOPFModel(pm_data, pm_form; rm_rsoc=true, rm_therm_line_lim=true)
    maxmin["model"] = convertSOCtoPSD(base_maxmin)
    JuMP.set_optimizer(maxmin["model"],with_optimizer(CPLEX.Optimizer))
    JuMP.set_parameter(maxmin["model"],"CPXPARAM_ScreenOutput",0)
    #JuMP.@constraint(maxmin["model"],objective_function(maxmin["model"], AffExpr) <= 1e4)
    gatherPSDConInfo(maxmin)


    feas=Dict{String,Any}()
    feas["branch_ids"] = pm_data["undecided_branches"]
    feas["x_soln"]=maxmin["x_soln"] #Shallow copy should be OK here
    feas["opt_val"]=0

    feas["model"],feas["ref_map"]=copy_model(maxmin["model"])
    JuMP.set_optimizer(feas["model"],pm_optimizer)
    gatherPSDConInfo(feas)
    relax_integrality(feas)

    removePSD_Constraints(maxmin) 
    lb_ids,ub_ids=add_artificial_var_bds(maxmin["model"]; bd_mag=1e1, io=io)
    println(io,maxmin["model"])
    #relax_integrality(maxmin)


    if using_eigs
        maxmin["cuts"]=Dict{SparseVector{Int64,Int64},Dict{Int64,Any}}()
	    add_initial_cuts(maxmin;io=io)
    else
        cut_pairs = Dict{Int64,Tuple{Any,Any}}()
	    add_initial_cuts(maxmin,feas,cut_pairs;io=io)
    end
    

    if !using_eigs
        cut_via_dual_feas(feas,maxmin,cut_pairs; add_cut=true, io=io)
    end
    resolveMP(maxmin; io=io)
    println(io,"Iteration 0: ")
    println("Iteration 0: ")
    print("\tUB: ",maxmin["UB_Val"]," with status: ",JuMP.termination_status(maxmin["model"])," with x solution: ")
    printX(maxmin["x_soln"])
    if !using_eigs
        println("\tLB: ",feas["opt_val"], " with status: ",JuMP.termination_status(feas["model"]),", best value is: ",bestLB)
    end

    new_cut = true
    if using_eigs
        new_cut = computeNewConstraintEig(maxmin;io=io)
    end


    for ii=1:MAX_N_ITER
        if !new_cut
	        println("Terminating at iteration $ii since the relaxed solution is feasible.")
            break
        end
        if !using_eigs
            cut_via_dual_feas(feas,maxmin,cut_pairs; add_cut=true, io=io)
            if feas["opt_val"] > bestLB
                bestLB=feas["opt_val"]
                best_x_soln=copy(maxmin["x_soln"])
                println(io,"Updating incumbent solution: ")
                printX(best_x_soln,io)
                println(io,"The incumbent value is now: ",bestLB,".")
                println("Updating incumbent solution: ")
                printX(best_x_soln)
                println("The incumbent value is now: ",bestLB,".")
            end
        end
        #artificial_lb_maxmin_vars = JuMP.all_variables(maxmin["model"])[lb_ids]
        #artificial_ub_maxmin_vars = JuMP.all_variables(maxmin["model"])[ub_ids]
        oldUB = maxmin["UB_Val"]
        resolveMP(maxmin; io=io)
        println(io,"Iteration $ii: ")
        println("Iteration $ii: ")
        print("\tUB: ",maxmin["UB_Val"]," with status: ",JuMP.termination_status(maxmin["model"])," with x solution: ")
        printX(maxmin["x_soln"])
        if !using_eigs
            println("\tLB: ",feas["opt_val"], " with status: ",JuMP.termination_status(feas["model"]),", best value is: ",bestLB)
        end
        if using_eigs
            new_cut = computeNewConstraintEig(maxmin;io=io)
        end
        if maxmin["UB_Val"] - bestLB < 1e-3
	        println("Terminating at iteration $ii since the relaxed solution is feasible.")
	        print("Final attack solution: ")
	        printX(best_x_soln)
	        println("With value: ", bestLB)
	        break
        end
        if oldUB + 1e-3 < maxmin["UB_Val"]
	        println("FLAGGING: oldUB ",oldUB, " < ", maxmin["UB_Val"],", somethings might be wrong.")
	        #break
        end
    end
    if !using_eigs
        cut_pair_keys = sort(collect(keys(cut_pairs)))
        for cc in cut_pair_keys
            println(io,"Cut $cc: ",cut_pairs[cc][1])
        end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Finishing after ", time_End," seconds.")
    for kk in keys(maxmin["cuts"])
        println(io,"For Bin ",kk,":")
        for ii in keys(maxmin["cuts"][kk])
            println(io,"\tSub-Bin ",ii, " has ",length(maxmin["cuts"][kk][ii])," cuts.")
        end
    end
end

function add_artificial_var_bds(model::JuMP.Model; bd_mag=1e6, io=Base.stdout)
    all_vars = JuMP.all_variables(model)
    n_vars = length(all_vars)
    var_ids = collect(1:n_vars)
    artificial_lb_var_ids=[]
    artificial_ub_var_ids=[]
    for vv in var_ids
        if !has_lower_bound(all_vars[vv])
            set_lower_bound(all_vars[vv], -bd_mag)
            #@constraint(model, all_vars[vv] >= -bd_mag)
            push!(artificial_lb_var_ids,vv)
            println(io,"Artificial LB of ",-bd_mag," for variable ",all_vars[vv])
        end
        if !has_upper_bound(all_vars[vv])
            set_upper_bound(all_vars[vv], bd_mag)
            #@constraint(model, all_vars[vv] <= bd_mag)
            push!(artificial_ub_var_ids,vv)
            println(io,"Artificial UB of ", bd_mag," for variable ",all_vars[vv])
        end
    end
    return artificial_lb_var_ids, artificial_ub_var_ids
end

function gatherPSDConInfo(model_dict)
    con_types=list_of_constraint_types(model_dict["model"])
    n_con_types=length(con_types)
    model_dict["psd_con"] = Dict{Int64,Dict{String,Any}}()
    nn = 0
    for cc=1:n_con_types
        if con_types[cc][2]==MathOptInterface.PositiveSemidefiniteConeTriangle
            psd_con = all_constraints(model_dict["model"], con_types[cc][1], con_types[cc][2]) 
            n_psd_con = length(psd_con)
            for kk=1:n_psd_con
                nn += 1
                model_dict["psd_con"][nn] = Dict{String,Any}()
                model_dict["psd_con"][nn]["ref"] = psd_con[kk]
                model_dict["psd_con"][nn]["expr"] = @expression(model_dict["model"], constraint_object(psd_con[kk]).func)
                model_dict["psd_con"][nn]["set"] = constraint_object(psd_con[kk]).set
                n_els = length(model_dict["psd_con"][nn]["expr"])
                model_dict["psd_con"][nn]["expr_val"] = zeros(n_els)
                model_dict["psd_con"][nn]["dual_val"] = zeros(n_els)
                model_dict["psd_con"][nn]["sg"] = Array{Float64,1}(undef, n_els)
            end
	    end
    end
end

function add_initial_cuts(maxmin,feas,cut_pairs; io=Base.stdout)
    n_expr = length(maxmin["psd_con"])
    @assert n_expr == length(feas["psd_con"])
    for nn=1:n_expr
        psd_expr=(maxmin["psd_con"][nn]["expr"],feas["psd_con"][nn]["expr"])
        n_els = length(psd_expr[1])
        @assert n_els == length(psd_expr[2])
	    jj,ii=1,1
        SG = Array{Float64,1}(undef,n_els)
        for mm in 1:n_els
            if typeof(psd_expr[1][mm])==VariableRef
                if JuMP.has_upper_bound(psd_expr[1][mm])
                    JuMP.delete_upper_bound( psd_expr[1][mm] ) ### aggregated UB constraint added later in function
                end
            end
            if ii==jj
                SG[mm] = 1
                if typeof(psd_expr[1][mm])==VariableRef
                    JuMP.set_lower_bound(psd_expr[1][mm], 0)
                else
                    JuMP.@constraint(maxmin["model"],psd_expr[1][mm] >= 0)
                end
                JuMP.@constraint(feas["model"],psd_expr[2][mm] >= 0)
                jj += 1
                ii = 1
            else
                SG[mm] = 2
                if typeof(psd_expr[1][mm])==VariableRef
                    JuMP.set_lower_bound( psd_expr[1][mm], -1e3 )
                else
                    JuMP.@constraint(maxmin["model"],psd_expr[1][mm] >= -1e3)
                end
                ii += 1
            end
        end
        JuMP.@constraint(maxmin["model"], sum( SG[mm]*psd_expr[1][mm] for mm in 1:n_els) <= 1e3)
        JuMP.@constraint(feas["model"], sum( SG[mm]*psd_expr[2][mm] for mm in 1:n_els) <= 1e3)
    end
    #delete_inactive_cuts(maxmin,feas,cut_pairs; io=io)
end

function add_initial_cuts(model_dict; io=Base.stdout)
    n_expr = length(model_dict["psd_con"])
    for nn=1:n_expr
        psd_expr=model_dict["psd_con"][nn]["expr"]
        n_els = length(psd_expr)
	    jj,ii=1,1
        SG = Array{Float64,1}(undef,n_els)
        for mm in 1:n_els
            if typeof(psd_expr[mm])==VariableRef
                if JuMP.has_upper_bound(psd_expr[mm])
                    JuMP.delete_upper_bound( psd_expr[mm] ) ### aggregated UB constraint added later in function
                end
            end
            if ii==jj
                SG[mm] = 1
                if typeof(psd_expr[mm])==VariableRef
                    JuMP.set_lower_bound(psd_expr[mm], 0)
                else
                    JuMP.@constraint(model_dict["model"],psd_expr[mm] >= 0)
                end
                jj += 1
                ii = 1
            else
                SG[mm] = 2
                if typeof(psd_expr[mm])==VariableRef
                    JuMP.set_lower_bound( psd_expr[mm], -1e3 )
                else
                    JuMP.@constraint(model_dict["model"],psd_expr[mm] >= -1e3)
                end
                ii += 1
            end
        end
        JuMP.@constraint(model_dict["model"], sum( SG[mm]*psd_expr[mm] for mm in 1:n_els) <= 1e3)
    end
end

function removePSD_Constraints(model_dict) 
    n_expr = length(model_dict["psd_con"])
    for nn=1:n_expr
        delete(model_dict["model"],model_dict["psd_con"][nn]["ref"])
        delete!(model_dict["psd_con"][nn],"ref")
    end
end

function relax_integrality(model_dict)
    for l in model_dict["branch_ids"]
	    if JuMP.is_integer(variable_by_name(model_dict["model"],"x[$l]_1"))
            JuMP.unset_integer(variable_by_name(model_dict["model"],"x[$l]_1"))
            JuMP.set_lower_bound(variable_by_name(model_dict["model"],"x[$l]_1"),0)
            JuMP.set_upper_bound(variable_by_name(model_dict["model"],"x[$l]_1"),1)
	    else
            JuMP.set_lower_bound(variable_by_name(model_dict["model"],"x[$l]_1"),0)
            JuMP.set_upper_bound(variable_by_name(model_dict["model"],"x[$l]_1"),1)
	    end
    end
end

function enforce_integrality(model_dict)
    for l in model_dict["branch_ids"]
        if has_lower_bound(variable_by_name(model_dict["model"],"x[$l]_1"))
            delete_lower_bound(variable_by_name(model_dict["model"],"x[$l]_1"))
        end
        if has_upper_bound(variable_by_name(model_dict["model"],"x[$l]_1"))
            delete_upper_bound(variable_by_name(model_dict["model"],"x[$l]_1"))
        end
        JuMP.set_integer(variable_by_name(model_dict["model"],"x[$l]_1"))
        set_lower_bound(variable_by_name(model_dict["model"],"x[$l]_1"),0)
        set_upper_bound(variable_by_name(model_dict["model"],"x[$l]_1"),1)
    end
end

function convertVecToSG(Vec::Array{Float64,1},SG::Array{Float64,1})
	n_el=length(Vec)
	@assert n_el == length(SG)
	ii,jj=1,1
    for nn=1:n_el
        if ii==jj
	        SG[nn] = Vec[nn]
	        jj += 1
            ii = 1
        else
	        SG[nn] = 2*Vec[nn]
            ii += 1
	    end
    end
end
function convertSGToVec(SG::Array{Float64,1},Vec::Array{Float64,1})
	n_el=length(Vec)
	@assert n_el == length(SG)
	jj=1
        for nn=1:n_el
          if nn == (Int)(((jj+1)*jj)/2)
	    Vec[nn] = SG[nn]
	    jj += 1
          else
	    Vec[nn] = 0.5*SG[nn]
	  end
	end
end

function cut_via_dual_feas(feas,maxmin,cut_pairs; add_cut=true, fix_x=true, io=Base.stdout)
    n_con=length(feas["psd_con"])
    if fix_x
      for l in feas["branch_ids"]
        JuMP.fix(variable_by_name(feas["model"],"x[$l]_1"),maxmin["x_soln"][l];force=true)
      end
    end
    JuMP.optimize!(feas["model"])
    for nn=1:n_con
	    feas["psd_con"][nn]["dual_val"][:] = JuMP.dual.(feas["psd_con"][nn]["ref"])[:]
    end
    
    feas["opt_val"] = JuMP.objective_value(feas["model"])
    println(io,"Optimal value of the dual feas problem is: ",feas["opt_val"], " with status: ",JuMP.termination_status(feas["model"]))

    if add_cut
      for nn=1:n_con
	    n_el=length(feas["psd_con"][nn]["dual_val"])
        convertVecToSG(feas["psd_con"][nn]["dual_val"],maxmin["psd_con"][nn]["sg"])
	    cut_pairs[ length(cut_pairs) ] = (
            @constraint(maxmin["model"], sum(maxmin["psd_con"][nn]["sg"][ii]*maxmin["psd_con"][nn]["expr"][ii] for ii=1:n_el)  >= 0),
            @constraint(feas["model"],   sum(maxmin["psd_con"][nn]["sg"][ii]*feas["psd_con"][nn]["expr"][ii] for ii=1:n_el)  >= 0)
        )
      end
    end

    if fix_x
      for l in feas["branch_ids"]
        JuMP.unfix(variable_by_name(feas["model"],"x[$l]_1"))
      end
    end
end

#=
function eval_cut_con(model_dict; io=Base.stdout)
  println(io,"Evaluating cut values: ")
  for cc in keys(model_dict["cuts"])
      println(io,"$cc => ",trunc(JuMP.value(model_dict["cuts"][cc]);digits=5))
  end
end
function eval_cut_dual(model_dict; io=Base.stdout)
  println(io,"Evaluating cut dual values: ")
  inactive_cut_ids=[]
  for cc in keys(model_dict["cuts"])
    if abs(JuMP.dual(model_dict["cuts"][cc])) < 1e-6
       push!(inactive_cut_ids,cc) 
    end
    println(io,"$cc => ",trunc(JuMP.dual(model_dict["cuts"][cc]);digits=5))
  end
  println(io,"Inactive cut indices: ",inactive_cut_ids)
  println(io,length(inactive_cut_ids)," sectiones delendae sunt ex ",length(model_dict["cuts"]), " sectionibus.")
  return inactive_cut_ids
end

function delete_inactive_cuts(maxmin,feas,cut_pairs; io=Base.stdout)
    relax_integrality(maxmin)
    JuMP.optimize!(maxmin["model"])
    println(io,"Solve status: ",JuMP.termination_status(maxmin["model"]))
    println(io,"Int-relaxed maxmin obj val:: ",JuMP.objective_value(maxmin["model"]))
    del_ids=[]
    for cc in keys(cut_pairs)
        if abs(JuMP.dual(cut_pairs[cc][1])) < 1e-6
            push!(del_ids,cc) 
        end
    end
    for dd in del_ids
        JuMP.delete(maxmin["model"],cut_pairs[dd][1])
        JuMP.delete(feas["model"],cut_pairs[dd][2])
    end
    old_list = copy(cut_pairs)
    empty!(cut_pairs)
    for kk in keys(old_list)
        if JuMP.is_valid(maxmin["model"],old_list[kk][1]) 
            @assert JuMP.is_valid(feas["model"],old_list[kk][2])
            cut_pairs[ length(cut_pairs) ] = old_list[kk]
        end
    end
    enforce_integrality(maxmin)
    return del_ids
end
=#

#=
function delete_inactive_cuts(model_dict,del_ids=[]; io=Base.stdout)
    if length(del_ids)==0
        del_ids = eval_cut_dual(model_dict; io=io)
    end
    for dd in del_ids
        JuMP.delete(model_dict["model"],model_dict["cuts"][dd])
    end
    nn=0
    old_list = copy(model_dict["cuts"])
    empty!(model_dict["cuts"])
    for kk in keys(old_list)
        if JuMP.is_valid(model_dict["model"],old_list[kk])
            model_dict["cuts"][nn]=old_list[kk]
            nn += 1
        end
    end
    return del_ids
end
=#

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
    @assert 2*n_el == (ncols+1)*ncols
    psd_mat = zeros(ncols,ncols)
    ii,jj=1,1
    for nn=1:n_el
        if ii==jj
	        psd_mat[ii,jj] = tri_mat[nn] 
            jj += 1; ii = 1
        else
	        psd_mat[ii,jj] = tri_mat[nn] 
	        psd_mat[jj,ii] = tri_mat[nn] 
            ii += 1
        end
    end
    E = eigs(psd_mat, which=:SR, nev=(ncols-1))
    return E
end

function computeNewConstraintEig(maxmin;io=Base.stdout)
    new_cut = false
    
    if !haskey(maxmin["cuts"],maxmin["x_soln_01"])
        maxmin["cuts"][ maxmin["x_soln_01"] ]=Dict{Int64,Any}()
    end
    cut_bin = maxmin["cuts"][ maxmin["x_soln_01"] ]
    for kk in keys(maxmin["psd_con"])
        if !haskey(cut_bin,kk)
            cut_bin[kk]=[]
        end
        println(io,"maxmin status $kk: ",JuMP.termination_status(maxmin["model"]))
        PSD_Vec = maxmin["psd_con"][kk]["expr"]
        PSD0 = maxmin["psd_con"][kk]["expr_val"]
        #println(io,"PSD0[$kk}: ",PSD0)
        E = computeEigenDecomp(PSD0)
        n_eigs = length(E[1])
        println(io,"eigenval_$kk: ", E[1])
        n_el = length(maxmin["psd_con"][kk]["expr"])
        SG = zeros(n_el)
        if E[1][1] < -1e-6
            new_cut = true
            #neg_eig_sum = sum(E[1][mm] for mm in filter(mmm->(E[1][mmm]) < 0, 1:n_eigs))
            for mm in 1:1 ###filter(mmm->(E[1][mmm] < 0), 1:n_eigs)
                #SG[1:n_el] .= 0
                ii,jj=1,1
                for nn=1:n_el
	                if ii==jj
	                    SG[nn] = E[2][ii,mm]*E[2][jj,mm]  
                        jj += 1
                        ii = 1
	                else
	                    SG[nn] = (E[2][ii,mm]*E[2][jj,mm] + E[2][jj,mm]*E[2][ii,mm]) 
                        ii += 1
	                end
                end
                push!( cut_bin[kk], @constraint(maxmin["model"], sum(SG[nn]*PSD_Vec[nn] for nn=1:n_el) >= 0))
            end
            #maxmin["cuts"][ length(maxmin["cuts"]) ] = @constraint(maxmin["model"], sum(SG[nn]*(PSD_Vec[nn]-PSD0[nn]) for nn=1:n_el) + neg_eig_sum   >= 0)
            #maxmin["cuts"][ length(maxmin["cuts"]) ] = @constraint(maxmin["model"], sum(SG[nn]*PSD_Vec[nn] for nn=1:n_el) >= 0)
        end
    end
    println(io,"End of computeNewConstraintEig() call.")
    return new_cut
end


function resolveMP(maxmin;io=Base.stdout)
  JuMP.optimize!(maxmin["model"])
  for kk in keys(maxmin["psd_con"])
    maxmin["psd_con"][kk]["expr_val"][:] = JuMP.value.(maxmin["psd_con"][kk]["expr"])[:]
  end
  maxmin["UB_Val"]=JuMP.objective_value(maxmin["model"])
  maxmin["x_soln_01"]=spzeros(Int64,maximum(maxmin["branch_ids"]))
  maxmin["x_soln"]=Dict{Int64,Float64}()
  for l in maxmin["branch_ids"]
    x_var = variable_by_name(maxmin["model"],"x[$l]_1")
    x_val = JuMP.value(x_var)
    if x_val > 1.0-1.0e-8
	    x_val = 1
        maxmin["x_soln_01"][l] = 1
    elseif x_val < 1.0e-8
	    x_val = 0
    end
    maxmin["x_soln"][l] = x_val
  end
end

function printX(x_soln::Dict{Int64,Float64},io=Base.stdout)
  for l in keys(x_soln)
      if x_soln[l] > 0.01
	if x_soln[l] > 0.99
	 print(io," ",l)
	else
	 print(io," ($l,",trunc(x_soln[l];digits=3),")")
	end
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

pm_optimizer=with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=8)

io = open(string(testcase["name"],".out"), "w")
#io = Base.stdout
#sdp_relax=[SDPWRMPowerModel, SparseSDPWRMPowerModel]
pm_form=SparseSDPWRMPowerModel
#pm_form=SDPWRMPowerModel
solveMaxminECP(pm_data,pm_form,pm_optimizer,io)
close(io)
