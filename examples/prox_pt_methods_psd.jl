include("../src/MaximinOPF.jl")
include("psd_utils.jl")
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


 # "This function should be implemented to be recallable," 
 ### "so that a user can call this function multiple times to get ever more accurate solution information."
 ### "Return aggregated cuts for the PSD constraints, and a Dictionary of new constraints"
function solve_PSD_via_ADMM(model_info::Dict{String,Any}; 
    max_n_iter=1000, prox_t=1, prim_tol=1e-3, dual_tol=1e-3, cs_tol=1e-3, use_cuts=false, rescale=false, display_freq::Int64=10,io=Base.stdout)
    model = model_info["model"] 
    model_info["prox_t"] = prox_t
    model_info["prox_t_min"] = 1e-3
    model_info["prox_t_max"] = 1e1
    model_info["proj_time"]=zeros(max_n_iter)
    model_info["qp_solve_time"]=zeros(max_n_iter)
    model_info["iter_time"]=zeros(max_n_iter)
    psd_expr = model[:psd_expr]
    prox_sign = model_info["prox_sign"] 
    PSD=model_info["psd_info"]

    max_violation = 1e20

    for kk in keys(PSD)
        PSD[kk]["prox_t"] = prox_t
        PSD[kk]["expr_val_ctr"][:] .= 0
        PSD[kk]["proj_expr_val"][:] .= 0
        PSD[kk]["C"][:] .= 0
    end
    for ii=1:max_n_iter
        iter_start_time = time_ns()
        for kk in keys(PSD)
            PSD[kk]["quad_terms"] = sum( 
                PSD[kk]["ip"][nn]*( ( psd_expr[kk,nn] - PSD[kk]["expr_val_ctr"][nn] + PSD[kk]["C"][nn]  )^2 - PSD[kk]["C"][nn]^2 ) 
                for nn in 1:PSD[kk]["vec_len"]) 
        end
        @objective( model, model_info["obj_sense"], model[:linobj_expr]  
            + 0.5* prox_sign * prox_t*sum(
                PSD[kk]["scale_factor"]*( PSD[kk]["quad_terms"] ) for kk in keys(PSD) )
        )

        qp_start_time = time_ns()
        try
            JuMP.optimize!( model) 
        catch exc
            println(exc)
	        println("Catching solve error, breaking from subiteration loop.")
            break
        end
        qp_end_time = time_ns()
        model_info["qp_solve_time"][ii] = (qp_end_time-qp_start_time)/1e9

        for kk in keys(PSD)
            for nn=1:PSD[kk]["vec_len"]
                PSD[kk]["expr_val"][nn] = JuMP.value(psd_expr[kk,nn])
            end
            for cc in keys(PSD[kk]["cuts"])
                PSD[kk]["cuts"][cc]["dual_val"] = abs(JuMP.dual(cc))
            end
        end
        model_info["prox_val"]=JuMP.objective_value(model)
        obj_val=JuMP.value(model[:linobj_expr]) 
        model_info["solve_status"]=JuMP.termination_status(model)

        proj_start_time = time_ns()
        ADMMProjections(model_info;io=io, validate=true)
        proj_end_time = time_ns()
        model_info["proj_time"][ii] = (proj_end_time-proj_start_time)/1e9

        prim_res = trunc(sqrt(sum( PSD[kk]["prim_res_norm"]^2 for kk in keys(PSD)));digits=4)
        max_prim = trunc(maximum( PSD[kk]["prim_res_norm"] for kk in keys(PSD));digits=4)

        dual_res = trunc(prox_t*sqrt(sum( PSD[kk]["dual_res_norm"]^2 for kk in keys(PSD)));digits=4)
        max_dual = trunc(prox_t*maximum( PSD[kk]["dual_res_norm"] for kk in keys(PSD));digits=4)

        #dual_res = trunc(prox_t*model_info["<C,X>"];digits=4)
        if mod(ii,display_freq)==0 || ii==1
            println("\tIter $ii obj value: ",
                round(obj_val;digits=4), 
                " (ALObj: ",round(model_info["prox_val"];digits=4),
                "), in statu: ", model_info["solve_status"],
                " p_res=",prim_res," d_res=",dual_res, 
                ", p<C,X>=",round(prox_t*model_info["<C,X>"];digits=4),
                " prox_t=",prox_t )
        elseif max(prim_res,dual_res) < max_violation
            max_violation = max(prim_res,dual_res)
            println("****Iter $ii prox obj value: ",
                round(obj_val;digits=4), 
                " (ALObj: ",round(model_info["prox_val"];digits=4),
                "), in statu: ", model_info["solve_status"],
                " p_res=",prim_res," d_res=",dual_res, 
                ", p<C,X>=",round(prox_t*model_info["<C,X>"];digits=4),
                " prox_t=",prox_t )
        end
        if prim_res < prim_tol
            if dual_res < dual_tol 

                println("********Iter $ii prox obj value: ",
                    round(obj_val;digits=4), 
                    " (ALObj: ",round(model_info["prox_val"];digits=4),
                    " in statu: ", model_info["solve_status"],
                    " p_res=",prim_res," d_res=",dual_res, 
                    ", p<C,X>=",round(prox_t*model_info["<C,X>"];digits=4),
                    " prox_t=",prox_t )
                if prox_t*model_info["<C,X>"] < cs_tol || !use_cuts
	                println("Sub-Iteratione $ii terminante, propter solutionem relaxatam est factibilem.")
                    iter_end_time = time_ns()
                    model_info["iter_time"][ii] = (iter_end_time-iter_start_time)/1e9
                    break
                else
                    add_C_cuts(model_info; delete_inactive=false)
                    println("\t\tAdding cuts at iteration $ii to improve sastifaction of optimality criteria.")
                    for kk in keys(PSD)
                        PSD[kk]["C"][:] .= 0
                    end
                end
            else  ### "dual_res >= dual_tol"
                if rescale && model_info["prox_t"] > model_info["prox_t_min"] 
                    model_info["prox_t"] *= 0.5
                    prox_t=round(model_info["prox_t"];digits=5)
                end
                #add_or_modify_C_cuts(model_info)
                if use_cuts
                    add_C_cuts(model_info; delete_inactive=false)
                    println("\t\tAdding cuts at iteration $ii to improve sastifaction of optimality criteria.")
                    for kk in keys(PSD)
                        PSD[kk]["C"][:] .= 0
                    end
                end
            end
        elseif dual_res < dual_tol  
            if rescale && model_info["prox_t"] < model_info["prox_t_max"] && prim_res >= prim_tol && mod(ii,100)==0
                model_info["prox_t"] *= 2
                prox_t=round(model_info["prox_t"];digits=5)
                println("\t\tIncreasing prox_t at iteration $ii to ",prox_t)
                for kk in keys(PSD)
                    PSD[kk]["C"][:] /= 2
                end
            end
        end
        iter_end_time = time_ns()
        model_info["iter_time"][ii] = (iter_end_time-iter_start_time)/1e9
        if mod(ii,100)==0 
            println("\t\tIter $ii: Last 100 iterations total time, qp solve time and proj time in secs: ",
                "\t",round(sum(model_info["iter_time"][jj] for jj in (ii-99):ii);digits=4),
                "\t",round(sum(model_info["qp_solve_time"][jj] for jj in (ii-99):ii);digits=4),
                "\t",round(sum(model_info["proj_time"][jj] for jj in (ii-99):ii);digits=4))
        end
        if mod(ii,1000)==0 
            println("\t\tIter $ii: Last 1000 iterations total time, qp solve time and proj time in secs: ",
                "\t",round(sum(model_info["iter_time"][jj] for jj in (ii-999):ii);digits=4),
                "\t",round(sum(model_info["qp_solve_time"][jj] for jj in (ii-999):ii);digits=4),
                "\t",round(sum(model_info["proj_time"][jj] for jj in (ii-999):ii);digits=4))
        end
        if mod(ii,10000)==0 
            println("\t\tIter $ii: Last 10000 iterations total time, qp solve time and proj time in secs: ",
                "\t",round(sum(model_info["iter_time"][jj] for jj in (ii-9999):ii);digits=4),
                "\t",round(sum(model_info["qp_solve_time"][jj] for jj in (ii-9999):ii);digits=4),
                "\t",round(sum(model_info["proj_time"][jj] for jj in (ii-9999):ii);digits=4))
        end
    end
end

function solveProxPtSP(model_info)
    PSD=model_info["psd_info"]
    model=model_info["model"]
    psd_expr = model[:psd_expr]
    for kk in keys(PSD)
        PSD[kk]["quad_terms"] = sum( PSD[kk]["ip"][nn]*( psd_expr[kk,nn] - PSD[kk]["expr_val_ctr"][nn] )^2 for nn in 1:PSD[kk]["vec_len"])
    end
    @objective( model, model_info["obj_sense"], model[:linobj_expr]  
        + 0.5* model_info["prox_sign"] * sum( PSD[kk]["prox_t"] * PSD[kk]["quad_terms"] for kk in keys(PSD) )
    )
    JuMP.optimize!( model) ### "Replace later with more simple solving technique, guaranteed to converge, like steepest descent, Newton-like?"

    for kk in keys(PSD)
        for nn=1:PSD[kk]["vec_len"]
            PSD[kk]["expr_val"][nn] = JuMP.value.(psd_expr[kk,nn])
        end
        for nn in keys(PSD[kk]["new_cuts"])
            PSD[kk]["new_cuts"][nn]["val"] = JuMP.value( PSD[kk]["new_cuts"][nn]["ref"] )
            PSD[kk]["new_cuts"][nn]["dual_val"] = JuMP.dual( PSD[kk]["new_cuts"][nn]["ref"] ) ###"Assumptions under which dual values exist"
        end
    end
    oldUB = model_info["opt_val"]
    model_info["opt_val"]=JuMP.objective_value(model)
    model_info["lin_objval"]=JuMP.value(model[:linobj_expr])
    model_info["solve_status"]=JuMP.termination_status(model)
end


 # "This function should be implemented to be recallable," 
 ### "so that a user can call this function multiple times to get ever more accurate solution information."
 ### "Return aggregated cuts for the PSD constraints, and a Dictionary of new constraints"
function solve_PSD_via_ProxPt(model_info::Dict{String,Any}; max_n_iter=100, prox_t=1, io=Base.stdout)
    psd_expr = model_info["model"][:psd_expr]
    prox_sign = model_info["prox_sign"] 
    PSD=model_info["psd_info"]

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
        PSDSubgradient(model_info;io=io)
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
        if model_info["total_sum_neg_eigs"] > -1e-4
	        println("Sub-Iteratione $ii terminante, quia solutio relaxata est factibilis.")
            break
        end
    end
end
