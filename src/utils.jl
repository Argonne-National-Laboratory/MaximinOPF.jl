using JuMP, MathOptInterface
using PowerModels

function add_artificial_var_bds(model::JuMP.Model; bd_mag=1e3, io=Base.stdout)
    all_vars = JuMP.all_variables(model)
    n_vars = length(all_vars)
    var_ids = collect(1:n_vars)
    artificial_lb_var_ids = []
    artificial_ub_var_ids = []
    for vv in var_ids
        if !is_integer(all_vars[vv]) && !is_binary(all_vars[vv])
            if !has_lower_bound(all_vars[vv]) 
                set_lower_bound(all_vars[vv], -bd_mag)
                push!(artificial_lb_var_ids,vv)
                #println(io,"Artificial LB of ",-bd_mag," for variable ",all_vars[vv])
            end
            if !has_upper_bound(all_vars[vv])
                set_upper_bound(all_vars[vv], bd_mag)
                push!(artificial_ub_var_ids,vv)
                #println(io,"Artificial UB of ", bd_mag," for variable ",all_vars[vv])
            end
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
    PSD=model_info["psd_con"]
    for kk in model_info["psd_con_keys"]
            PSD[kk] = Dict{String,Any}()
            PSD[kk]["all_cuts"] = Dict{Int64,Dict{String,Any}}() ### "Will be productive later"
            PSD[kk]["new_cuts"] = Dict{Int64,Dict{String,Any}}()
            PSD[kk]["vec_len"] = length(psd_con_expr[kk])
            vec_len = PSD[kk]["vec_len"] 
            PSD[kk]["expr_val"] = zeros(vec_len)
            PSD[kk]["expr_val_ctr"] = zeros(vec_len)
            PSD[kk]["old_proj_expr_val"] = zeros(vec_len)
            PSD[kk]["proj_expr_val"] = zeros(vec_len)
            PSD[kk]["C"] = zeros(vec_len) ### Dual solution
            PSD[kk]["prim_res"] = zeros(vec_len)
            PSD[kk]["dual_res"] = zeros(vec_len)
            PSD[kk]["sg"] = Array{Float64,1}(undef, vec_len)
            PSD[kk]["eig_val"] = 0.0
            PSD[kk]["ssc_val"] = 0.0
            PSD[kk]["ip"] = ones(vec_len) ### "coefficients to aid in computing Frobenius inner product"
            PSD[kk]["is_off_diag"] = zeros(vec_len) ### "coefficients to aid in computing Frobenius inner product"
            PSD[kk]["ij_coor"] = Vector{Tuple{Int64,Int64}}(undef,vec_len)
	        jj,ii=1,1
            for mm in 1:vec_len
                PSD[kk]["ij_coor"][mm]=(ii,jj)
                if ii==jj
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
            PSD[kk]["n_matcols"]=jj-1
            PSD[kk]["quad_terms"] = AffExpr(0.0)
            PSD[kk]["quad_term_vals"] = 0.0
            PSD[kk]["Lagr_terms"] = AffExpr(0.0)
            PSD[kk]["Lagr_term_vals"] = 0.0
    end
end

function add_psd_initial_cuts(model::JuMP.Model; io=Base.stdout)
    model_info = Dict{String,Any}()
    model_info["model"] = model
    gatherPSDConInfo(model_info)
    psd_expr = model_info["model"][:psd_expr]
    PSD=model_info["psd_con"]
    #JuMP.@constraint( model_info["model"], psd_expr_lbs[kk in keys(PSD), mm in 1:PSD[kk]["vec_len"]], psd_expr[kk,mm] >= -PSD[kk]["is_off_diag"][mm]*1e2 )
    JuMP.@constraint( model_info["model"], psd_expr_lbs[kk in keys(PSD), mm in 1:PSD[kk]["vec_len"]], psd_expr[kk,mm] >= -1e2 )
    JuMP.@constraint( model_info["model"], psd_expr_ubs[kk in keys(PSD), mm in 1:PSD[kk]["vec_len"]], psd_expr[kk,mm] <= 1e2 )
    #JuMP.@constraint( model_info["model"], psd_expr_lbs[kk in keys(PSD), mm in filter(mmm->(!PSD[kk]["is_off_diag"][mmm]),1:PSD[kk]["vec_len"])], psd_expr[kk,mm] >= 0 )

    #JuMP.@constraint(model_info["model"], sum( sum( PSD[kk]["ip"][mm]*psd_expr[kk,mm] for mm in 1:PSD[kk]["vec_len"]) for kk in keys(PSD) ) <= 1e2)
end

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

function unfix_vars(model,branch_ids)
    for l in branch_ids
        if is_fixed(variable_by_name(model,"x[$l]_1"))
            unfix(variable_by_name(model,"x[$l]_1"))
        end
        set_lower_bound(variable_by_name(model,"x[$l]_1"),0)
        set_upper_bound(variable_by_name(model,"x[$l]_1"),1)
    end
end
