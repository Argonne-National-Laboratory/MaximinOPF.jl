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

### "Assum the model has linear and PSD constraints"
### "Assume model already has a subproblem solver attached to it, "
###    "and that problem modifications will not otherwise cause reference to expressions and constraints to break"
function prepare_to_solve_PSD_via_ProxPt(model::JuMP.Model; io=Base.stdout)
    time_Start = time_ns()

    model_info=Dict{String,Any}()
    model_info["linobj_expr"] = objective_function(model, AffExpr)
    model_info["obj_sense"] = objective_sense(model)
    model_info["opt_val"]=1e20
    model_info["lin_objval"]=1e20
    model_info["lin_objval_ctr"]=1e20
    model_info["model"]=model

    gatherPSDConInfo(model_info)
    removePSD_Constraints(model_info)
    #lb_ids,ub_ids=add_artificial_var_bds(model_info["model"]; bd_mag=1e1, io=io)  ### Necesse est ut problema relaxatum esset delimitum
    JuMP.@constraint(model_info["model"], obj_ub, model_info["linobj_expr"] <= 1e3)
    JuMP.@constraint(model_info["model"], obj_lb, model_info["linobj_expr"] >= -1e3)

	add_psd_initial_cuts(model_info;io=io)
    
    time_End = (time_ns()-time_Start)/1e9
    println("Initializing finished after ", time_End," seconds.")
    return model_info
end

 # "This function should be implemented to be recallable," 
 ### "so that a user can call this function multiple times to get ever more accurate solution information."
 ### "Return aggregated cuts for the PSD constraints, and a Dictionary of new constraints"
function solve_PSD_via_ProxPt(model_info::Dict{String,Any}; io=Base.stdout)
    MAX_N_ITER = 10 
    PSD=model_info["psd_con"]
    for kk in keys(model_info["psd_con"])
        empty!( PSD[kk]["new_cuts"])
    end
    n_del=0

    @objective( model_info["model"], model_info["obj_sense"], model_info["linobj_expr"] ) 
    JuMP.optimize!( model_info["model"] ) 
    for kk in keys(model_info["psd_con"])
        PSD[kk]["expr_val"][:] = JuMP.value.(PSD[kk]["expr"])[:]
        for nn in keys(PSD[kk]["new_cuts"])
            PSD[kk]["new_cuts"][nn]["val"] = JuMP.value( PSD[kk]["new_cuts"][nn]["ref"] )
            PSD[kk]["new_cuts"][nn]["dual_val"] = JuMP.dual( PSD[kk]["new_cuts"][nn]["ref"] ) ###"Assumptions under which dual values exist"
        end
    end
    model_info["lin_objval"]=JuMP.objective_value(model_info["model"])
    model_info["lin_objval_ctr"] = model_info["lin_objval"]

    exists_new_cut = computeNewConstraintEig(model_info;io=io)
    for kk in keys(PSD)
        PSD[kk]["expr_val_ctr"][:] = PSD[kk]["expr_val"][:] 
        PSD[kk]["eig_val_ctr"] = PSD[kk]["eig_val"]
    end

    prox_sign = 1
    if model_info["obj_sense"]==MOI.MAX_SENSE
        prox_sign = -1
    end
    prox_t = 10
    @objective( model_info["model"], model_info["obj_sense"], model_info["linobj_expr"] 
        + 0.5* prox_sign * prox_t * sum( sum( ( PSD[kk]["expr"][nn] - PSD[kk]["expr_val_ctr"][nn] )^2 for nn in keys(PSD[kk]["expr"])) for kk in keys(PSD))
    )

    for ii=0:MAX_N_ITER
        if !exists_new_cut
	        println("Sub-Iteratione $ii terminante, quia solutio relaxata est factibilis.")
            break
        end

        JuMP.optimize!( model_info["model"] ) ### "Replace later with more simple solving technique, guaranteed to converge, like steepest descent, Newton-like?"
        for kk in keys(model_info["psd_con"])
            PSD[kk]["expr_val"][:] = JuMP.value.(PSD[kk]["expr"])[:]
            for nn in keys(PSD[kk]["new_cuts"])
                PSD[kk]["new_cuts"][nn]["val"] = JuMP.value( PSD[kk]["new_cuts"][nn]["ref"] )
                PSD[kk]["new_cuts"][nn]["dual_val"] = JuMP.dual( PSD[kk]["new_cuts"][nn]["ref"] ) ###"Assumptions under which dual values exist"
            end
        end
        oldUB = model_info["opt_val"]
        model_info["opt_val"]=JuMP.objective_value(model_info["model"])
        model_info["lin_objval"]=JuMP.value(model_info["linobj_expr"])
        model_info["solve_status"]=JuMP.termination_status(model_info["model"])

        println("\tProx obj value: ",model_info["lin_objval"]," in statu: ", model_info["solve_status"])

        ### "TODO: incorporate ssc update of proximal center"
#### "How does this relate to ALM, what happens if all constraints are relaxed in an ALM framework; ADMM??"
        exists_new_cut = computeNewConstraintEig(model_info;io=io)
        mp_obj_val = model_info["lin_objval"] 
        obj_val_ctr = model_info["lin_objval_ctr"] + 100*sum( PSD[kk]["eig_val_ctr"] for kk in keys(PSD) ) 
        obj_val =     model_info["lin_objval"]     + 100*sum( PSD[kk]["eig_val"]     for kk in keys(PSD) ) 
        if (obj_val - obj_val_ctr) >= 0.1*(mp_obj_val - obj_val_ctr)
            model_info["lin_objval_ctr"] = model_info["lin_objval"]
            for kk in keys(PSD)
                PSD[kk]["expr_val_ctr"][:] = PSD[kk]["expr_val"][:] 
                PSD[kk]["eig_val_ctr"] = PSD[kk]["eig_val"]
            end
            @objective( model_info["model"], model_info["obj_sense"], model_info["linobj_expr"] 
                + 0.5* prox_sign * prox_t * sum( sum( ( PSD[kk]["expr"][nn] - PSD[kk]["expr_val_ctr"][nn] )^2 for nn in keys(PSD[kk]["expr"])) for kk in keys(PSD))
            )
        end
println("\t\tMin eig value is: ", minimum( PSD[kk]["eig_val"] for kk in keys(PSD)), ",\tmin center eig value is: ", minimum( PSD[kk]["eig_val_ctr"] for kk in keys(PSD)))
        ### "the above update of the objective function will be conditional once SSC is implemented"
    end
    ###TODO: Aggregate the set of new cuts into an aggregate
#=
    for kk in keys(PSD)
        #println("$kk cut vals: ")
        for nn in keys(PSD[kk]["new_cuts"])
        #print(" ", has_duals( model_info["model"] ))
            if abs( PSD[kk]["new_cuts"][nn]["val"] ) >= 1e-4
                JuMP.delete(model_info["model"], PSD[kk]["new_cuts"][nn]["ref"] )
                n_del += 1
            end
        end
        #print("\n")
    end
=#
    println("\t\tNumber of cuts deleted: ",n_del)
end

function add_artificial_var_bds(model::JuMP.Model; bd_mag=1e3, io=Base.stdout)
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
                model_dict["psd_con"][nn]["all_cuts"] = Dict{Int64,Dict{String,Any}}() ### "Will be productive later"
                model_dict["psd_con"][nn]["new_cuts"] = Dict{Int64,Dict{String,Any}}()
                model_dict["psd_con"][nn]["ref"] = psd_con[kk]
                model_dict["psd_con"][nn]["expr"] = @expression(model_dict["model"], constraint_object(psd_con[kk]).func)
                model_dict["psd_con"][nn]["set"] = constraint_object(psd_con[kk]).set
                n_els = length(model_dict["psd_con"][nn]["expr"])
                model_dict["psd_con"][nn]["expr_val"] = zeros(n_els)
                model_dict["psd_con"][nn]["expr_val_ctr"] = zeros(n_els)
                model_dict["psd_con"][nn]["proj_expr_val"] = zeros(n_els)
                model_dict["psd_con"][nn]["sg"] = Array{Float64,1}(undef, n_els)
                model_dict["psd_con"][nn]["eig_val"] = 0.0
                model_dict["psd_con"][nn]["ssc_val"] = 0.0
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

function computeNewConstraintEig(model_info;io=Base.stdout)
    new_cut = false
    
    for kk in keys(model_info["psd_con"])
        PSD = model_info["psd_con"][kk]
        PSD_Vec = PSD["expr"]
        PSD0 = PSD["expr_val"]
        # PSD["proj_expr_val"][:] = PSD["expr_val"][:]

        E = computeEigenDecomp(PSD0)
        n_eigs = length(E[1])
        println(io,"eigenval_$kk: ", E[1])
        n_el = length(PSD["expr"])
        SG = zeros(n_el)
        PSD["eig_val"] = E[1][1]
        if E[1][1] < -1e-6
            new_cut = true
            for mm in 1:1 # filter(mmm->(E[1][mmm] < 0), 1:n_eigs)
                ii,jj=1,1
                for nn=1:n_el
                    # PSD["proj_expr_val"][nn] -= E[1][mm]*E[2][ii,mm]*E[2][jj,mm]
	                if ii==jj
	                    SG[nn] = E[2][ii,mm]*E[2][jj,mm]  
                        jj += 1
                        ii = 1
	                else
	                    SG[nn] = (E[2][ii,mm]*E[2][jj,mm] + E[2][jj,mm]*E[2][ii,mm]) 
                        ii += 1
	                end
                    PSD["sg"][nn] = SG[nn]
                end
                ncuts=length(PSD["new_cuts"])
                PSD["new_cuts"][ncuts] = Dict{String,Any}()
                PSD["new_cuts"][ncuts]["ref"] = @constraint(model_info["model"], sum( SG[nn]*PSD_Vec[nn] for nn=1:n_el ) >= 0) 
            end
        end
    end
    return new_cut
end

function ProjectToPSD(model_info;io=Base.stdout)
    is_nontrivial_projection = false 
    for kk in keys(model_info["psd_con"])
        PSD_Vec = model_info["psd_con"][kk]["expr"]
        PSD0 = model_info["psd_con"][kk]["expr_val"]
        model_info["psd_con"][kk]["proj_expr_val"][:] = model_info["psd_con"][kk]["expr_val"][:]

        E = computeEigenDecomp(PSD0)
        n_eigs = length(E[1])
        println(io,"eigenval_$kk: ", E[1])
        n_el = length(model_info["psd_con"][kk]["expr"])
        SG = zeros(n_el)
        if E[1][1] < -1e-6
            is_nontrivial_projection = true
            new_cut = true
            for mm in filter(mmm->(E[1][mmm] < 0), 1:n_eigs)
                ii,jj=1,1
                for nn=1:n_el
                    model_info["psd_con"][kk]["proj_expr_val"][nn] -= E[1][mm]*E[2][ii,mm]*E[2][jj,mm]
                end
            end
        end
    end
    return is_nontrivial_projection
end
