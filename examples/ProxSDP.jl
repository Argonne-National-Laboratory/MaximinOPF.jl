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

### "Synopsis:"
#### "How does this relate to ALM, what happens if all constraints are relaxed in an ALM framework; ADMM??"
####    "Because we do not have a closed-form evaluation of the minimum eigenvalue function, we cannot apply ALM"
####    "Rather, the ALM approach is replaced with a proximal bundle method (PBM) approach; there are different PBM approaches to evaluate"
####    "Also we can project a Matrix to a PSD matrix w.r.t. the Frobenius inner product using eigendecomposition, so ADMM is applicable and is implemented also."

### "Assum the model has linear and PSD constraints"
### "Assume model already has a subproblem solver attached to it, "
###    "and that problem modifications will not otherwise cause reference to expressions and constraints to break"
function prepare_to_solve_PSD_via_ProxPt(model::JuMP.Model, model_info::Dict{String,Any}=Dict{String,Any}(); io=devnull)
    time_Start = time_ns()

    @expression(model, linobj_expr, objective_function(model, AffExpr))
    model_info["model"] = model
    model_info["opt_val"]=1e20
    model_info["lin_objval"]=1e20
    model_info["lin_objval_ctr"]=1e20
    model_info["cp_model"]=model
    model_info["prox_t_min"] = 1e-3
    model_info["prox_t_max"] = 1e3

    model_info["obj_sense"] = objective_sense(model)
    model_info["prox_sign"] = 1
    if model_info["obj_sense"]==MOI.MAX_SENSE
        model_info["prox_sign"] = -1
    end

    gatherPSDConInfo(model_info) ### "Sets the 'psd_info' entry of model_info"
	add_psd_initial_cuts(model_info;io=io)
    removePSD_Constraints(model_info["psd_info"])
    
    time_End = (time_ns()-time_Start)/1e9
    println("Initializing finished after ", time_End," seconds.")
    return model_info
end

#=
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
function add_artificial_var_bds(model_info::Dict{String,Any}; bd_mag=1e0, io=Base.stdout)
    model = model_info["model"]
    all_vars = JuMP.all_variables(model)
    n_vars = length(all_vars)
    var_ids = collect(1:n_vars)
    artificial_lb_var_ids=[]
    artificial_ub_var_ids=[]
    for vv in var_ids
        if has_lower_bound(all_vars[vv])
            if !has_upper_bound(all_vars[vv])
                set_upper_bound(all_vars[vv], bd_mag)
                push!(artificial_ub_var_ids,vv)
            end
        elseif has_upper_bound(all_vars[vv])
            if !has_lower_bound(all_vars[vv])
                set_lower_bound(all_vars[vv], -bd_mag)
                push!(artificial_lb_var_ids,vv)
            end
        end
    end
    model_info["artificial_lb_var_ids"] = artificial_lb_var_ids
    model_info["artificial_ub_var_ids"] = artificial_ub_var_ids
end
=#
function add_artificial_var_bds(model_info; bd_mag=1, io=devnull)
    linobj_expr = objective_function(model_info["model"], AffExpr)
    for tt in keys(linobj_expr.terms)
        if objective_sense(model_info["model"])==MOI.MAX_SENSE
            if linobj_expr.terms[tt] > 0
                if !has_upper_bound(tt)
                    set_upper_bound(tt,1)
                end
            elseif linobj_expr.terms[tt] < 0
                if !has_lower_bound(tt)
                    set_lower_bound(tt,-1)
                end
            end
        elseif objective_sense(model_info["model"])==MOI.MIN_SENSE
            if linobj_expr.terms[tt] > 0
                if !has_lower_bound(tt)
                    set_lower_bound(tt,-1)
                end
            elseif linobj_expr.terms[tt] < 0
                if !has_upper_bound(tt)
                    set_upper_bound(tt,1)
                end
            end
        else
            println("Minimization sense is abnormal: ", objective_sense(model_info["model"]) )
        end
    end
end

function gatherPSDConInfo(model_info::Dict{String,Any})
    model_info["psd_info"] = gatherPSDConInfo(model_info["model"])
end

function gatherPSDConInfo(model::JuMP.Model)
    con_types=list_of_constraint_types(model)
    n_con_types=length(con_types)
    PSD = Dict{Tuple{Int64,Int64},Dict{String,Any}}()
    psd_con_type_ids = filter( cc->(con_types[cc][2]==MathOptInterface.PositiveSemidefiniteConeTriangle), 1:n_con_types)
    for cc in psd_con_type_ids
        psd_con = all_constraints(model, con_types[cc][1], con_types[cc][2]) 
        n_con = length(psd_con)
        for nn in 1:n_con
            PSD[cc,nn] = Dict{String,Any}()
            PSD[cc,nn]["cref"] = psd_con[nn]
            PSD[cc,nn]["expr"] = constraint_object(psd_con[nn]).func
            PSD[cc,nn]["vec_len"] = length(PSD[cc,nn]["expr"])
        end
    end
    @expression(model, psd_expr[kk in keys(PSD), mm=1:PSD[kk]["vec_len"]], PSD[kk]["expr"][mm] )
    for kk in keys(PSD)
            PSD[kk]["model"] = model
            PSD[kk]["all_cuts"] = Dict{Int64,Dict{String,Any}}() ### "Will be productive later"
            PSD[kk]["new_cuts"] = Dict{Int64,Dict{String,Any}}()
            vec_len = PSD[kk]["vec_len"] 
            PSD[kk]["expr_val"] = zeros(vec_len)
            PSD[kk]["dual_val"] = zeros(vec_len)
            PSD[kk]["expr_val_ctr"] = zeros(vec_len)
            PSD[kk]["old_proj_expr_val"] = zeros(vec_len)
            PSD[kk]["proj_expr_val"] = zeros(vec_len)
            PSD[kk]["C"] = zeros(vec_len) ### Dual solution
            PSD[kk]["prim_res"] = zeros(vec_len)
            PSD[kk]["prim_res_norm"] = 0.0
            PSD[kk]["dual_res"] = zeros(vec_len)
            PSD[kk]["dual_res_norm"] = 0.0
            PSD[kk]["prox_t"] = 1
            PSD[kk]["sg"] = zeros(vec_len)
            PSD[kk]["ij_pairs"] = Array{Tuple{Int64,Int64},1}(undef, vec_len)
            PSD[kk]["ii"] = zeros(Int64,vec_len)
            PSD[kk]["jj"] = zeros(Int64,vec_len)
            PSD[kk]["diag_ids"] = Dict{Int64,Int64}()
            PSD[kk]["min_eigval"] = 0.0
            PSD[kk]["ssc_val"] = 0.0
            PSD[kk]["ip"] = ones(vec_len) ### "coefficients to aid in computing Frobenius inner product"
            PSD[kk]["is_off_diag"] = zeros(vec_len) ### "coefficients to aid in computing Frobenius inner product"
	        jj,ii=1,1
            for mm in 1:vec_len
                PSD[kk]["ij_pairs"][mm]=(ii,jj)
                PSD[kk]["ii"][mm]=ii
                PSD[kk]["jj"][mm]=jj
                if ii==jj
                    PSD[kk]["diag_ids"][jj] = mm
                    PSD[kk]["ip"][mm] = 1
                    PSD[kk]["is_off_diag"][mm] = 0
                    jj += 1
                    ii = 1
                else
                    PSD[kk]["ip"][mm] = 2
                    PSD[kk]["is_off_diag"][mm] = 1
                    ii += 1
                end
            end
            PSD[kk]["ncols"] = PSD[kk]["jj"][ PSD[kk]["vec_len"] ]
            PSD[kk]["psd_mat"] = zeros( PSD[kk]["ncols"], PSD[kk]["ncols"] )
            PSD[kk]["diag_els"] = filter(mmm->(PSD[kk]["is_off_diag"][mmm] == 0),1:PSD[kk]["vec_len"])
            PSD[kk]["off_diag_els"] = filter(mmm->(PSD[kk]["is_off_diag"][mmm] == 1),1:PSD[kk]["vec_len"])
            PSD[kk]["quad_terms"] = AffExpr(0.0)
            PSD[kk]["quad_term_vals"] = 0.0
            PSD[kk]["Lagr_terms"] = AffExpr(0.0)
            PSD[kk]["Lagr_term_vals"] = 0.0
    end
    return PSD
end

function add_psd_initial_cuts(model_info; bdmag=1e2, io=Base.stdout)
    model=model_info["model"]
    psd_expr = model[:psd_expr]
    PSD = model_info["psd_info"]
    JuMP.@constraint( model, psd_lbs[kk in keys(PSD), mm in 1:PSD[kk]["vec_len"]], psd_expr[kk,mm] >= -bdmag*PSD[kk]["is_off_diag"][mm] )
    JuMP.@constraint( model, psd_ubs[kk in keys(PSD), mm in 1:PSD[kk]["vec_len"]], sum( psd_expr[kk,mm] for mm in 1:PSD[kk]["vec_len"]) <= bdmag )

#=
    JuMP.@constraint( psd_info["cp_model"], psd_diag_lbs[kk in keys(PSD), mm in PSD[kk]["diag_els"]], psd_expr[kk,mm] >= 0.0 )
    JuMP.@constraint( psd_info["cp_model"], psd_diag_ubs[kk in keys(PSD), mm in PSD[kk]["diag_els"]], psd_expr[kk,mm] <= bdmag )

    JuMP.@constraint( psd_info["cp_model"], psd_off_diag_con1[kk in keys(PSD), mm in PSD[kk]["off_diag_els"]], 
              psd_expr[kk, PSD[kk]["diag_ids"][ PSD[kk]["ij_pairs"][mm][1] ]  ] + 2*psd_expr[kk,mm] 
            + psd_expr[kk, PSD[kk]["diag_ids"][ PSD[kk]["ij_pairs"][mm][2] ]  ] >= 0.0 )
    JuMP.@constraint( psd_info["cp_model"], psd_off_diag_con2[kk in keys(PSD), mm in PSD[kk]["off_diag_els"]], 
              psd_expr[kk, PSD[kk]["diag_ids"][ PSD[kk]["ij_pairs"][mm][1] ]  ] - 2*psd_expr[kk,mm] 
            + psd_expr[kk, PSD[kk]["diag_ids"][ PSD[kk]["ij_pairs"][mm][2] ]  ] >= 0.0 )
=#
#=
    for kk in psd_info["psd_con_keys"]
        n_els = psd_info["psd_con"][kk]["vec_len"] 
        IP = psd_info["psd_con"][kk]["ip"]
        for mm in 1:n_els
            JuMP.@constraint(psd_info["cp_model"], psd_expr_lbs[kk psd_expr[kk,mm] >= 0)
            if IP[mm]==1
                JuMP.@constraint(psd_info["cp_model"],psd_expr[kk,mm] >= 0)
            else
                JuMP.@constraint(psd_info["cp_model"],psd_expr[kk,mm] >= -1e3)
            end
        end
    end
    JuMP.@constraint(psd_info["cp_model"], sum( sum( PSD[kk]["ip"][mm]*psd_expr[kk,mm] for mm in 1:PSD[kk]["vec_len"]) for kk in keys(PSD) ) >= -1e3)
=#
end

function removePSD_Constraints(PSD) 
    for kk in keys(PSD)
        delete(PSD[kk]["model"],PSD[kk]["cref"])
    end
end
 # "This function should be implemented to be recallable," 
 ### "so that a user can call this function multiple times to get ever more accurate solution information."
 ### "Return aggregated cuts for the PSD constraints, and a Dictionary of new constraints"
function solve_PSD_via_ADMM(model_info::Dict{String,Any}; max_n_iter=100, prox_t=1, io=Base.stdout)
    model_info["model"] = model
    psd_expr = model[:psd_expr]
    prox_sign = model_info["prox_sign"] 
    PSD=model_info["psd_info"]

    for kk in keys(PSD)
        PSD[kk]["prox_t"] = prox_t
        PSD[kk]["proj_expr_val"][:] .= 0
        PSD[kk]["C"][:] .= 0
    end
    for ii=0:max_n_iter
        for kk in keys(PSD)
            PSD[kk]["quad_terms"] = sum( PSD[kk]["ip"][nn]*( psd_expr[kk,nn] - PSD[kk]["proj_expr_val"][nn] )^2 for nn in 1:PSD[kk]["vec_len"]) 
            PSD[kk]["Lagr_terms"] = sum( PSD[kk]["ip"][nn]*PSD[kk]["C"][nn]*( psd_expr[kk,nn] - PSD[kk]["proj_expr_val"][nn] )  for nn in 1:PSD[kk]["vec_len"]) 
        end
        @objective( model, model_info["obj_sense"], model[:linobj_expr]  
            - prox_sign * sum( PSD[kk]["Lagr_terms"] for kk in keys(PSD))
            + 0.5* prox_sign * sum( PSD[kk]["prox_t"] * PSD[kk]["quad_terms"] for kk in keys(PSD) )
        )
        try
            JuMP.optimize!( model) 
        catch exc
            println(exc)
	        println("Catching solve error, breaking from subiteration loop.")
            break
        end
        pd_res_val = 0.0
        for kk in keys(PSD)
            for nn=1:PSD[kk]["vec_len"]
                PSD[kk]["expr_val"][nn] = JuMP.value(psd_expr[kk,nn])
            end
            PSD[kk]["quad_term_vals"] = JuMP.value(PSD[kk]["quad_terms"])
            pd_res_val += PSD[kk]["quad_term_vals"] 
        end
        model_info["prox_val"]=JuMP.objective_value(model)
        model_info["solve_status"]=JuMP.termination_status(model)

        is_nontrivial_projection = PSDProjections(model_info;io=io)
        prim_res = trunc(sqrt(sum( PSD[kk]["prim_res_norm"]^2 for kk in keys(PSD)));digits=4)
        dual_res = trunc(sqrt(sum( PSD[kk]["dual_res_norm"]^2 for kk in keys(PSD)));digits=4)
        println("\tsub-Iter $ii prox obj value: ",model_info["prox_val"]," in statu: ", model_info["solve_status"]," p_res=",prim_res," d_res=",dual_res )
        if sqrt(pd_res_val) < 1e-4
	        println("Sub-Iteratione $ii terminante, quia solutio relaxata est factibilis.")
            break
        else
#=
            if mod(ii,100)==40
                for kk in keys(model_info["psd_con"])
                    if PSD[kk]["prim_res_norm"] > 10*PSD[kk]["dual_res_norm"]
                        PSD[kk]["prox_t"] *= 2
                    elseif PSD[kk]["dual_res_norm"] > 10*PSD[kk]["prim_res_norm"]
                        PSD[kk]["prox_t"] /= 2
                    end
                    PSD[kk]["prox_t"] = min( max(model_info["prox_t_min"],PSD[kk]["prox_t"]), model_info["prox_t_max"])
                end
            end
=#
        end
    end
    
end

function solveProxPtSP(model_info)
    PSD=model_info["psd_con"]
    psd_expr = model_info["cp_model"][:psd_expr]
    for kk in keys(model_info["psd_con"])
        PSD[kk]["quad_terms"] = sum( PSD[kk]["ip"][nn]*( psd_expr[kk,nn] - PSD[kk]["expr_val_ctr"][nn] )^2 for nn in 1:PSD[kk]["vec_len"])
    end
    @objective( model_info["cp_model"], model_info["obj_sense"], model_info["cp_model"][:linobj_expr]  
        + 0.5* model_info["prox_sign"] * sum( PSD[kk]["prox_t"] * PSD[kk]["quad_terms"] for kk in keys(PSD) )
    )
    JuMP.optimize!( model_info["cp_model"] ) ### "Replace later with more simple solving technique, guaranteed to converge, like steepest descent, Newton-like?"

    for kk in keys(model_info["psd_con"])
        for nn=1:PSD[kk]["vec_len"]
            PSD[kk]["expr_val"][nn] = JuMP.value.(psd_expr[kk,nn])
        end
        for nn in keys(PSD[kk]["new_cuts"])
            PSD[kk]["new_cuts"][nn]["val"] = JuMP.value( PSD[kk]["new_cuts"][nn]["ref"] )
            PSD[kk]["new_cuts"][nn]["dual_val"] = JuMP.dual( PSD[kk]["new_cuts"][nn]["ref"] ) ###"Assumptions under which dual values exist"
        end
    end
    oldUB = model_info["opt_val"]
    model_info["opt_val"]=JuMP.objective_value(model_info["cp_model"])
    model_info["lin_objval"]=JuMP.value(model_info["cp_model"][:linobj_expr])
    model_info["solve_status"]=JuMP.termination_status(model_info["cp_model"])
end

function conditionallyUpdateProxCenter(model_info)
    PSD=model_info["psd_info"]
    mp_obj_val = model_info["lin_objval"] 
    obj_val_ctr = model_info["lin_objval_ctr"] + 100*sum( PSD[kk]["eigval_ctr"] for kk in keys(PSD) ) 
    obj_val =     model_info["lin_objval"]     + 100*sum( PSD[kk]["min_eigval"]     for kk in keys(PSD) ) 
    if (obj_val - obj_val_ctr) >= 0.1*(mp_obj_val - obj_val_ctr) 
        model_info["lin_objval_ctr"] = model_info["lin_objval"]
        for kk in keys(PSD)
            PSD[kk]["expr_val_ctr"][:] = PSD[kk]["expr_val"][:] 
            PSD[kk]["eigval_ctr"] = PSD[kk]["min_eigval"]
        end
    end
end

 # "This function should be implemented to be recallable," 
 ### "so that a user can call this function multiple times to get ever more accurate solution information."
 ### "Return aggregated cuts for the PSD constraints, and a Dictionary of new constraints"
function solve_PSD_via_ProxPt(model_info::Dict{String,Any}; max_n_iter=100, prox_t=1, io=Base.stdout)
    psd_expr = model_info["model"][:psd_expr]
    prox_sign = model_info["prox_sign"] 
    PSD=model_info["psd_con"]

    for kk in keys(PSD)
        PSD[kk]["prox_t"] = prox_t
        empty!( PSD[kk]["new_cuts"])
    end

    for kk in keys(PSD)
        PSD[kk]["expr_val_ctr"][:] .= 0 
        PSD[kk]["eigval_ctr"] = 0
        PSD[kk]["C"][:] .= 0
    end
    model_info["lin_objval_ctr"] = model_info["prox_sign"]*1e20

    for ii=0:max_n_iter
        try
            solveProxPtSP(model_info)
        catch exc
            println(exc)
            println("Catching solve error, breaking from subiteration loop.")
            break
        end

        println("\tsub-Iter $ii Prox obj value: ",model_info["lin_objval"]," in statu: ", model_info["solve_status"])

        ### "TODO: incorporate ssc update of proximal center"
        total_negeigs_sum = PSDSubgradient(PSD;io=io)
        for kk in keys(PSD)
            ncuts=length(PSD[kk]["new_cuts"])
            PSD[kk]["new_cuts"][ncuts] = Dict{String,Any}()
            PSD[kk]["new_cuts"][ncuts]["ref"] = @constraint(PSD[kk]["model"], sum( PSD[kk]["ip"][nn]*PSD[kk]["sg"][nn]*psd_expr[kk,nn] for nn=1:PSD[kk]["vec_len"] ) >= 0) 
        end
        if prox_t > 0
            conditionallyUpdateProxCenter(model_info)
        end
        println("\t\tMin eig value is: ", minimum( PSD[kk]["min_eigval"] for kk in keys(PSD)), 
            ",\tmin center eig value is: ", minimum( PSD[kk]["eigval_ctr"] for kk in keys(PSD)))
        aggregateNewCuts(model_info)
        if total_negeigs_sum < -1e-4
	        println("Sub-Iteratione $ii terminante, quia solutio relaxata est factibilis.")
            break
        end
    end
end

function aggregateNewCuts(model_info)
    ###TODO: Aggregate the set of new cuts into an aggregate
    n_del = 0
    PSD=model_info["psd_info"]["psd_con"]
    for kk in keys(PSD)
#=
        for nn in keys(PSD[kk]["new_cuts"])
        #print(" ", has_duals( model_info["model"] ))
            if abs( PSD[kk]["new_cuts"][nn]["val"] ) >= 1e-4
                JuMP.delete(model_info["model"], PSD[kk]["new_cuts"][nn]["ref"] )
                n_del += 1
            end
        end
        #print("\n")
=#
    end
    println("\t\tNumber of cuts deleted: ",n_del)

end



###Input is symmetric matrix in triangular form
function computePSDMat(PSD, tri_mat)
    for nn=1:PSD["vec_len"]
        ii,jj=PSD["ii"][nn],PSD["jj"][nn]
	    PSD["psd_mat"][ii,jj] = tri_mat[nn]
        if ii != jj
	        PSD["psd_mat"][jj,ii] = PSD["psd_mat"][ii,jj]
        end
    end
end

function PSDSubgradient(model_info; io=Base.stdout)
    model_info["total_sum_neg_eigs"] = 0
    model = model_info["model"]
    psd_expr = model[:psd_expr]
    PSD = model_info["psd_info"]
    sg_info = Dict{Tuple{Int64,Int64},Array{Float64,1}}()
    for kk in keys(PSD)
        sg_info[kk] = zeros(PSD[kk]["vec_len"])
        computePSDMat(PSD[kk], PSD[kk]["expr_val"])
        E = eigen(PSD[kk]["psd_mat"])
        eig_vals,eig_vecs = E
#=
        E = eigs(PSD[kk]["psd_mat"], which=:SR, nev=1)
        eig_vals,eig_vecs = E[1],E[2]
=#
        PSD[kk]["neg_eigs_sum"] = 0
        neg_eigs = filter(mmm->(eig_vals[mmm] < -1e-8), 1:length(eig_vals) ) 
        if length(neg_eigs) > 0
            PSD[kk]["neg_eigs_sum"] = sum( eig_vals[mm] for mm in neg_eigs )
            model_info["total_sum_neg_eigs"] += PSD[kk]["neg_eigs_sum"]
        end
        PSD[kk]["min_eigval"],min_idx = findmin(eig_vals)
        for nn=1:PSD[kk]["vec_len"]
            ii,jj=PSD[kk]["ii"][nn],PSD[kk]["jj"][nn]
	        PSD[kk]["sg"][nn] = eig_vecs[ii,min_idx]*eig_vecs[jj,min_idx]  
            sg_info[kk][nn] = PSD[kk]["sg"][nn]
        end
    end
    return sg_info
end

function PSDProjections(model_info; io=Base.stdout)
    is_nontrivial_projection = false 
    psd_expr = model_info["model"][:psd_expr]
    prox_sign = model_info["prox_sign"]
    PSD = model_info["psd_info"]
    for kk in keys(PSD)
        PSD[kk]["old_proj_expr_val"][:] = PSD[kk]["proj_expr_val"][:]
        PSD[kk]["proj_expr_val"][:] = PSD[kk]["expr_val"][:] - (1.0/PSD[kk]["prox_t"])*PSD[kk]["C"][:]
        computePSDMat(PSD[kk], PSD[kk]["expr_val"])
        E = eigen(PSD[kk]["psd_mat"])
        eig_vals,eig_vecs = E

        n_eigs = length(eig_vals)
        neg_eigs = filter(mmm->(eig_vals[mmm] < -1e-8), 1:n_eigs)
        PSD[kk]["neg_eigs_sum"] = sum( eig_vals[mm] for mm in neg_eigs )
        println(io,"eigenval_$kk: ", eig_vals)
        PSD[kk]["min_eigval"],min_eigval_idx = findmin(eig_vals)
        if PSD[kk]["min_eigval"] < 0
            is_nontrivial_projection = true
        end
        for nn=1:PSD[kk]["vec_len"]
            ii,jj=PSD[kk]["ii"][nn],PSD[kk]["jj"][nn]
            PSD[kk]["proj_expr_val"][nn] -= sum( eig_vals[mm]*eig_vecs[ii,mm]*eig_vecs[jj,mm] for mm in neg_eigs)
            PSD[kk]["C"][nn] -= PSD[kk]["prox_t"]*(PSD[kk]["expr_val"][nn]-PSD[kk]["proj_expr_val"][nn])
            PSD[kk]["prim_res"][nn] = PSD[kk]["expr_val"][nn] - PSD[kk]["proj_expr_val"][nn]
            PSD[kk]["dual_res"][nn] = PSD[kk]["prox_t"]*(PSD[kk]["proj_expr_val"][nn] - PSD[kk]["old_proj_expr_val"][nn])
        end
        PSD[kk]["prim_res_norm"] = sqrt(sum( PSD[kk]["prim_res"][nn]^2 for nn=1:PSD[kk]["vec_len"]))
        PSD[kk]["dual_res_norm"] = sqrt(sum( PSD[kk]["dual_res"][nn]^2 for nn=1:PSD[kk]["vec_len"]))
    end
    return is_nontrivial_projection
end
