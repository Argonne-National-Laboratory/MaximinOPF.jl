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
    @expression(model, linobj_expr, objective_function(model, AffExpr))
    model_info["obj_sense"] = objective_sense(model)
    model_info["opt_val"]=1e20
    model_info["lin_objval"]=1e20
    model_info["lin_objval_ctr"]=1e20
    model_info["model"]=model
    model_info["prox_t"] = 1
    model_info["prox_sign"] = 1
    if model_info["obj_sense"]==MOI.MAX_SENSE
        model_info["prox_sign"] = -1
    end

    gatherPSDConInfo(model_info)
    removePSD_Constraints(model_info)
    #lb_ids,ub_ids=add_artificial_var_bds(model_info["model"]; bd_mag=1e1, io=io)  ### Necesse est ut problema relaxatum esset delimitum
    JuMP.@constraint(model_info["model"], obj_ub, linobj_expr <= 1e3)
    JuMP.@constraint(model_info["model"], obj_lb, linobj_expr >= -1e3)

	add_psd_initial_cuts(model_info;io=io)
    
    time_End = (time_ns()-time_Start)/1e9
    println("Initializing finished after ", time_End," seconds.")
    return model_info
end
 # "This function should be implemented to be recallable," 
 ### "so that a user can call this function multiple times to get ever more accurate solution information."
 ### "Return aggregated cuts for the PSD constraints, and a Dictionary of new constraints"
function solve_PSD_via_ADMM(model_info::Dict{String,Any}; max_n_iter=100, prox_t=1, io=Base.stdout)
    psd_expr = model_info["model"][:psd_expr]
    model_info["prox_t"] = prox_t
    prox_sign = model_info["prox_sign"] 
    PSD=model_info["psd_con"]
    quad_terms = Dict{Tuple{Int64,Int64},Any}()
    Lagr_terms = Dict{Tuple{Int64,Int64},Any}()

    for kk in keys(model_info["psd_con"])
        PSD[kk]["proj_expr_val"][:] .= 0
        PSD[kk]["C"][:] .= 0
    end

    

    for ii=0:max_n_iter
        for kk in keys(model_info["psd_con"])
            quad_terms[kk] = sum( PSD[kk]["ip"][nn]*( psd_expr[kk,nn] - PSD[kk]["proj_expr_val"][nn] )^2 for nn in 1:PSD[kk]["vec_len"]) 
            Lagr_terms[kk] = sum( PSD[kk]["ip"][nn]*PSD[kk]["C"][nn]*( psd_expr[kk,nn] - PSD[kk]["proj_expr_val"][nn] )  for nn in 1:PSD[kk]["vec_len"]) 
        end
        @objective( model_info["model"], model_info["obj_sense"], model_info["model"][:linobj_expr]  
            - prox_sign * sum( Lagr_terms[kk] for kk in keys(PSD))
            + 0.5* prox_sign * prox_t * sum( quad_terms[kk] for kk in keys(PSD) )
        )
        try
            JuMP.optimize!( model_info["model"] ) 
        catch exc
            println(exc)
	        println("Catching solve error, breaking from subiteration loop.")
            break
        end
        for kk in keys(model_info["psd_con"])
            n_el = PSD[kk]["vec_len"]
            for nn=1:n_el
                PSD[kk]["expr_val"][nn] = JuMP.value(psd_expr[kk,nn])
            end
            #println("Quad $kk: ",JuMP.value(quad_terms[kk]))
            #println("Lagr $kk: ",JuMP.value(Lagr_terms[kk]))
        end
        model_info["prox_val"]=JuMP.objective_value(model_info["model"])
        model_info["solve_status"]=JuMP.termination_status(model_info["model"])

        is_nontrivial_projection = ADMMProjections(model_info;io=io)
        if mod(ii,20)==0
            println("\tProx obj value: ",model_info["prox_val"]," in statu: ", model_info["solve_status"])
        end
        if !is_nontrivial_projection
	        println("Sub-Iteratione $ii terminante, quia solutio relaxata est factibilis.")
            break
        end
    end
    
end

 # "This function should be implemented to be recallable," 
 ### "so that a user can call this function multiple times to get ever more accurate solution information."
 ### "Return aggregated cuts for the PSD constraints, and a Dictionary of new constraints"
function solve_PSD_via_ProxPt(model_info::Dict{String,Any}; io=Base.stdout)
    MAX_N_ITER = 6 
    psd_expr = model_info["model"][:psd_expr]

    PSD=model_info["psd_con"]
    for kk in keys(model_info["psd_con"])
        empty!( PSD[kk]["new_cuts"])
    end
    n_del = 0

    @objective( model_info["model"], model_info["obj_sense"], model_info["model"][:linobj_expr] ) 
    JuMP.optimize!( model_info["model"] ) 
    for kk in keys(model_info["psd_con"])
        n_el = model_info["psd_con"][kk]["vec_len"]
        for nn=1:n_el
            PSD[kk]["expr_val"][nn] = JuMP.value.(psd_expr[kk,nn])
        end
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
    prox_t = 1
    model_info["prox_t"] = prox_t
    
    @objective( model_info["model"], model_info["obj_sense"], model_info["model"][:linobj_expr]  
        + 0.5* prox_sign * prox_t * sum( sum( ( psd_expr[kk,nn] - PSD[kk]["expr_val_ctr"][nn] )^2 for nn in 1:length(PSD[kk]["expr_val"])) for kk in keys(PSD))
    )

    for ii=0:MAX_N_ITER
        if !exists_new_cut
	        println("Sub-Iteratione $ii terminante, quia solutio relaxata est factibilis.")
            break
        end

        try
            JuMP.optimize!( model_info["model"] ) ### "Replace later with more simple solving technique, guaranteed to converge, like steepest descent, Newton-like?"
        catch exc
            println(exc)
	        println("Catching solve error, breaking from subiteration loop.")
            break
        end

        for kk in keys(model_info["psd_con"])
            n_el = model_info["psd_con"][kk]["vec_len"]
            for nn=1:n_el
                PSD[kk]["expr_val"][nn] = JuMP.value.(psd_expr[kk,nn])
            end
            for nn in keys(PSD[kk]["new_cuts"])
                PSD[kk]["new_cuts"][nn]["val"] = JuMP.value( PSD[kk]["new_cuts"][nn]["ref"] )
                PSD[kk]["new_cuts"][nn]["dual_val"] = JuMP.dual( PSD[kk]["new_cuts"][nn]["ref"] ) ###"Assumptions under which dual values exist"
            end
        end
        oldUB = model_info["opt_val"]
        model_info["opt_val"]=JuMP.objective_value(model_info["model"])
        model_info["lin_objval"]=JuMP.value(model_info["model"][:linobj_expr])
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
            @objective( model_info["model"], model_info["obj_sense"], model_info["model"][:linobj_expr] 
                + 0.5* prox_sign * prox_t * sum( sum( ( psd_expr[kk,nn] - PSD[kk]["expr_val_ctr"][nn] )^2 for nn in 1:length(PSD[kk]["expr_val"])) for kk in keys(PSD))
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

function gatherPSDConInfo(model_info)
    con_types=list_of_constraint_types(model_info["model"])
    n_con_types=length(con_types)
    model_info["psd_con"] = Dict{Tuple{Int64,Int64},Dict{String,Any}}()
    model_info["psd_con_type_ids"] = filter( cc->(con_types[cc][2]==MathOptInterface.PositiveSemidefiniteConeTriangle), 1:n_con_types)
    model_info["psd_con_keys"]=[]
    psd_con = Dict{Int64,Any}()
    psd_con_expr = Dict{Tuple{Int64,Int64},Any}()
    for cc in model_info["psd_con_type_ids"]
        psd_con[cc] = all_constraints(model_info["model"], con_types[cc][1], con_types[cc][2]) 
        n_con = length(psd_con[cc])
        for nn in 1:n_con
            push!(model_info["psd_con_keys"],(cc,nn))    
            psd_con_expr[cc,nn] = constraint_object(psd_con[cc][nn]).func
        end
    end
    @expression(model_info["model"], psd_expr[kk in model_info["psd_con_keys"], nn=1:length(psd_con_expr[kk])], psd_con_expr[kk][nn] )
    for kk in model_info["psd_con_keys"]
            model_info["psd_con"][kk] = Dict{String,Any}()
            model_info["psd_con"][kk]["all_cuts"] = Dict{Int64,Dict{String,Any}}() ### "Will be productive later"
            model_info["psd_con"][kk]["new_cuts"] = Dict{Int64,Dict{String,Any}}()
            model_info["psd_con"][kk]["vec_len"] = length(psd_con_expr[kk])
            vec_len = model_info["psd_con"][kk]["vec_len"] 
            model_info["psd_con"][kk]["expr_val"] = zeros(vec_len)
            model_info["psd_con"][kk]["expr_val_ctr"] = zeros(vec_len)
            model_info["psd_con"][kk]["proj_expr_val"] = zeros(vec_len)
            model_info["psd_con"][kk]["C"] = zeros(vec_len) ### Dual solution
            model_info["psd_con"][kk]["sg"] = Array{Float64,1}(undef, vec_len)
            model_info["psd_con"][kk]["eig_val"] = 0.0
            model_info["psd_con"][kk]["ssc_val"] = 0.0
            model_info["psd_con"][kk]["ip"] = ones(vec_len) ### "coefficients to aid in computing Frobenius inner product"
	        jj,ii=1,1
            for mm in 1:vec_len
            if ii==jj
                model_info["psd_con"][kk]["ip"][mm] = 1
                jj += 1
                ii = 1
            else
                model_info["psd_con"][kk]["ip"][mm] = 2
                ii += 1
            end
        end
            
    end
end

function add_psd_initial_cuts(model_info; io=Base.stdout)
    psd_expr = model_info["model"][:psd_expr]
    for kk in model_info["psd_con_keys"]
        n_els = model_info["psd_con"][kk]["vec_len"] 
        IP = model_info["psd_con"][kk]["ip"]
        for mm in 1:n_els
            if IP[mm]==1
                JuMP.@constraint(model_info["model"],psd_expr[kk,mm] >= 0)
            else
                JuMP.@constraint(model_info["model"],psd_expr[kk,mm] >= -1e3)
            end
        end
        JuMP.@constraint(model_info["model"], sum( IP[mm]*psd_expr[kk,mm] for mm in 1:n_els) <= 1e3)
    end
end

function removePSD_Constraints(model_info) 
    con_types=list_of_constraint_types(model_info["model"])
    for cc in model_info["psd_con_type_ids"]
        psd_con = all_constraints(model_info["model"], con_types[cc][1], con_types[cc][2]) 
        n_con = length(psd_con)
        for nn in 1:n_con
            delete(model_info["model"],psd_con[nn])
        end
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
    #E = eigs(psd_mat, which=:SR, nev=ncols)
    E = eigen(psd_mat)
    return E
end

function computeNewConstraintEig(model_info;io=Base.stdout)
    new_cut = false
    psd_expr = model_info["model"][:psd_expr]
    
    for kk in keys(model_info["psd_con"])
        PSD = model_info["psd_con"][kk]
        PSD0 = PSD["expr_val"]
        # PSD["proj_expr_val"][:] = PSD["expr_val"][:]

        E = computeEigenDecomp(PSD0)
        eig_vals,eig_vecs = E
        n_eigs = length(eig_vals)
        println(io,"eigenval_$kk: ", eig_vals)
        n_el = length(PSD["expr_val"])
        SG = zeros(n_el)
        min_eigval,min_eigval_idx = findmin(eig_vals)
        PSD["eig_val"] = min_eigval
        if min_eigval < -1e-6
            new_cut = true
            mm = min_eigval_idx
            #for mm in 1:1 # filter(mmm->(E[1][mmm] < 0), 1:n_eigs)
                ii,jj=1,1
                for nn=1:n_el
	                if ii==jj
	                    SG[nn] = eig_vecs[ii,mm]*eig_vecs[2][jj,mm]  
                        jj += 1
                        ii = 1
	                else
	                    SG[nn] = (eig_vecs[ii,mm]*eig_vecs[2][jj,mm] + eig_vecs[2][jj,mm]*eig_vecs[2][ii,mm]) 
                        ii += 1
	                end
                    PSD["sg"][nn] = SG[nn]
                end
                ncuts=length(PSD["new_cuts"])
                PSD["new_cuts"][ncuts] = Dict{String,Any}()
                PSD["new_cuts"][ncuts]["ref"] = @constraint(model_info["model"], sum( SG[nn]*psd_expr[kk,nn] for nn=1:n_el ) >= 0) 
            #end
        end
    end
    return new_cut
end

function PSDProjection(model_info;io=Base.stdout)
    is_nontrivial_projection = false 
    psd_expr = model_info["model"][:psd_expr]
    prox_sign = model_info["prox_sign"] 
    prox_t = model_info["prox_t"]
    for kk in keys(model_info["psd_con"])
        PSD = model_info["psd_con"][kk]
        PSD0 = copy(PSD["expr_val"])
        n_el = PSD["vec_len"]

        E = computeEigenDecomp(PSD0)
        eig_vals,eig_vecs = E
        n_eigs = length(eig_vals)
        neg_eigs = filter(mmm->(eig_vals[mmm] < 0), 1:n_eigs)
        println(io,"eigenval_$kk: ", eig_vals)
        min_eigval,min_eigval_idx = findmin(eig_vals)
        if min_eigval < -1e-6
            is_nontrivial_projection = true
        end
        PSD["proj_expr_val"][:] = PSD0[:]
        for mm in neg_eigs
            ii,jj=1,1
            for nn=1:n_el
                PSD["proj_expr_val"][nn] -= eig_vals[mm]*eig_vecs[ii,mm]*eig_vecs[jj,mm]
	            if ii==jj
                    jj += 1
                    ii = 1
	            else
                    ii += 1
	            end
            end
        end
    end
    return is_nontrivial_projection
end

function ADMMProjections(model_info;io=Base.stdout)
    is_nontrivial_projection = false 
    psd_expr = model_info["model"][:psd_expr]
    prox_sign = model_info["prox_sign"] 
    prox_t = model_info["prox_t"]
    for kk in keys(model_info["psd_con"])
        PSD = model_info["psd_con"][kk]
        PSD0 = copy(PSD["expr_val"])
        n_el = PSD["vec_len"]
        for nn=1:n_el
            PSD0[nn] -= (1.0/prox_t)*PSD["C"][nn]
        end

        E = computeEigenDecomp(PSD0)
        eig_vals,eig_vecs = E
        n_eigs = length(eig_vals)
        neg_eigs = filter(mmm->(eig_vals[mmm] < 0), 1:n_eigs)
        println(io,"eigenval_$kk: ", eig_vals)
        min_eigval,min_eigval_idx = findmin(eig_vals)
        if min_eigval < -1e-6
            is_nontrivial_projection = true
        end
        PSD["proj_expr_val"][:] = PSD0[:]
        for mm in neg_eigs
            ii,jj=1,1
            for nn=1:n_el
                PSD["proj_expr_val"][nn] -= eig_vals[mm]*eig_vecs[ii,mm]*eig_vecs[jj,mm]
	            if ii==jj
                    jj += 1
                    ii = 1
	            else
                    ii += 1
	            end
            end
        end
        for nn=1:n_el
            PSD["C"][nn] -= prox_t*(PSD["expr_val"][nn]-PSD["proj_expr_val"][nn])
        end
    end
    return is_nontrivial_projection
end
