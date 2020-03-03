include("../../MaximinOPF/src/MaximinOPF.jl")
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

### Assume model already has a subproblem solver attached to it
### Return an aggregated cut for the PSD constraint, and a Dictionary of new constraints
function prepare_to_solve_PSD_via_ProxPt(model::JuMP.Model; io=Base.stdout)
    model_info=Dict{String,Any}()
    model_info["linobj_expr"] = objective_function(model, AffExpr)
    model_info["obj_sense"] = objective_sense(model)
    model_info["opt_val"]=0
    model_info["model"]=model

    gatherPSDConInfo(model_info)
    removePSD_Constraints(model_info)
    lb_ids,ub_ids=add_artificial_var_bds(model_info["model"]; bd_mag=1e1, io=io)  ### Necesse est ut problema relaxatum esset delimitum


    #model_info["cuts"]=Dict{SparseVector{Int64,Int64},Dict{Int64,Any}}()
    model_info["new_cuts"]=Dict{Int64,Vector{Any}}()
    for kk in keys(model_info["psd_con"])
        model_info["new_cuts"][kk] = []
    end
    @objective(model, model_info["obj_sense"], 
        model_info["linobj_expr"] + sum( (model_info["psd_con"][kk][nn] - 0)^2 for (kk,nn) in (keys(model_info["psd_con"]),keys(model_info["psd_con"][kk]))) )

	add_psd_initial_cuts(model_info;io=io)
    
    JuMP.optimize!(model_info["model"])
    for kk in keys(model_info["psd_con"])
        model_info["psd_con"][kk]["expr_val"][:] = JuMP.value.(model_info["psd_con"][kk]["expr"])[:]
    end
    model_info["opt_val"]=JuMP.objective_value(model_info["model"])

    println(io,"Iteratione 0: ")
    println("Iteratione 0: ")
    print("\tUB est: ", model_info["opt_val"]," in statu: ",JuMP.termination_status(model_info["model"]))
    new_cut = computeNewConstraintEig(model_info;io=io)

    for ii=1:MAX_N_ITER
        if !new_cut
	        println("Iteratione $ii terminante, quia solutio relaxata est factibilis.")
            break
        end

        JuMP.optimize!(model_info["model"])
        for kk in keys(model_info["psd_con"])
            model_info["psd_con"][kk]["expr_val"][:] = JuMP.value.(model_info["psd_con"][kk]["expr"])[:]
        end
        oldUB = model_info["opt_val"]
        model_info["opt_val"]=JuMP.objective_value(model_info["model"])

        println(io,"Iteration $ii: ")
        println("Iteration $ii: ")
        print("\tUB: ",model_info["UB_Val"]," in statu: ",JuMP.termination_status(model_info["model"])," cum solutione x: ")
        new_cut = computeNewConstraintEig(maxmin;io=io)
        if oldUB + 1e-3 < maxmin["UB_Val"]
	        println("FLAGGING: oldUB ",oldUB, " < ", maxmin["UB_Val"],", aliquid esset erratum.")
	        #break
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
    return model_info
end

function continue_to_solvePSD_via_ProxPt(model_info::Dict{String,Any}; io=Base.stdout)
    empty!(model_info["new_cuts"])

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
                model_dict["psd_con"][nn]["sg"] = Array{Float64,1}(undef, n_els)
            end
	    end
    end
end

function add_psd_initial_cuts(model_dict; io=Base.stdout)
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
