module MaximinOPF
using PowerModels
using Dualization
using MathOptInterface

include("Variables.jl")
include("Objectives.jl")
include("Constraints.jl")
include("utils.jl")
greet() = print("Hello World!")

supported_pm=[  ACPPowerModel, ACRPowerModel, ACTPowerModel, 
		SOCWRPowerModel, SOCWRConicPowerModel, 
		SOCBFPowerModel, SOCBFConicPowerModel, 
		QCRMPowerModel, QCLSPowerModel, 
		SDPWRMPowerModel, SparseSDPWRMPowerModel,
                DCPPowerModel, DCMPPowerModel, NFAPowerModel,
		DCPLLPowerModel, LPACCPowerModel ]
conic_supported_pm=[ SOCWRConicPowerModel, SOCBFConicPowerModel, 
		SDPWRMPowerModel, SparseSDPWRMPowerModel,
                DCPPowerModel, DCMPPowerModel, NFAPowerModel]
sdp_pm=[ SDPWRMPowerModel, SparseSDPWRMPowerModel ]
nonconvex_pm=[ACPPowerModel, ACRPowerModel, ACTPowerModel]
not_supported_pm=[IVRPowerModel]


function MaximinOPFModel(pm_data, pm_form; rm_rsoc=true, rm_therm_line_lim=false)
    println("Calling MaximinOPFModel() with powerform: ",pm_form)
    if pm_form in conic_supported_pm 
        minmax_pm = MinimaxOPFModel(pm_data, pm_form)
	if rm_rsoc
	  removeRSOC(minmax_pm)
	end
	if pm_form in conic_supported_pm && rm_therm_line_lim
	  removeThermalLineLimits(minmax_pm)
	end

        #Test Out
        io = open(string(pm_data["name"],".out"), "w")
        println(io,"name: ", pm_data["name"])
        println(io,"attacker_budget: ", pm_data["attacker_budget"])
        println(io,"inactive_branches: ", pm_data["inactive_branches"])
        println(io,"protected_branches: ", pm_data["protected_branches"])
        println(io, "Model:")
        println(io, minmax_pm.model)
        close(io)

        maxmin_model = DualizeMinmaxModel(minmax_pm)
        for l in minmax_pm.data["undecided_branches"]
            if has_lower_bound(variable_by_name(maxmin_model,"x[$l]_1"))
                delete_lower_bound(variable_by_name(maxmin_model,"x[$l]_1"))
            end
            if has_upper_bound(variable_by_name(maxmin_model,"x[$l]_1"))
                delete_upper_bound(variable_by_name(maxmin_model,"x[$l]_1"))
            end
            JuMP.set_integer(variable_by_name(maxmin_model,"x[$l]_1"))
            set_lower_bound(variable_by_name(maxmin_model,"x[$l]_1"),0)
            set_upper_bound(variable_by_name(maxmin_model,"x[$l]_1"),1)
        end
	return maxmin_model
    else
        println("Model type: ",pm_form," is either nonconic and/or nonconvex and thus is not currently supported by function MaximinOPFModel().")
	println("WARNING: function MaximinOPFModel() is returning 'nothing'")
	return nothing
    end
end

function removeRSOC(pm)
    con_types=list_of_constraint_types(pm.model)
    n_con_types=length(con_types)
    for cc=1:n_con_types
        if con_types[cc][2]==MathOptInterface.RotatedSecondOrderCone
	  ### These constraints are presumably associated with defining the quadratic cost function for the OPF and are usually not needed here.
	    rsoc_con = all_constraints(pm.model, con_types[cc][1], con_types[cc][2]) 
	    n_rsoc_con = length(rsoc_con)
	    for nn=1:n_rsoc_con
		JuMP.delete(pm.model,rsoc_con[nn])
	    end
	end
    end
end

function removeThermalLineLimits(pm)
    con_types=list_of_constraint_types(pm.model)
    n_con_types=length(con_types)
    for cc=1:n_con_types
        if con_types[cc][2]==MathOptInterface.SecondOrderCone
	    #println("SOC:")
	    soc_con = all_constraints(pm.model, con_types[cc][1], con_types[cc][2]) 
	    n_soc_con = length(soc_con)
	    for nn=1:n_soc_con
		soc_expr = constraint_object(soc_con[nn]).func 
		n_vars=length( soc_expr ) 
		if n_vars == 3  ### This condition identifies the thermal line limits in SOC form
		  JuMP.delete(pm.model,soc_con[nn])
		end
	    end
	end
    end
end

"As precondition, assume that minmax_model is in dualizable form and is of a supported form"
function MaximinOPFModel(minmax_model::JuMP.Model, line_idx=[])
        maxmin_model = DualizeMinmaxModel(minmax_model)
        for l in line_idx
            if has_lower_bound(variable_by_name(maxmin_model,"x[$l]_1"))
                delete_lower_bound(variable_by_name(maxmin_model,"x[$l]_1"))
            end
            if has_upper_bound(variable_by_name(maxmin_model,"x[$l]_1"))
                delete_upper_bound(variable_by_name(maxmin_model,"x[$l]_1"))
            end
            JuMP.set_integer(variable_by_name(maxmin_model,"x[$l]_1"))
            set_lower_bound(variable_by_name(maxmin_model,"x[$l]_1"),0)
            set_upper_bound(variable_by_name(maxmin_model,"x[$l]_1"),1)
        end
	return maxmin_model
end

function MinimaxOPFModel(pm_data, powerform)
    if powerform in supported_pm
      pm = instantiate_model(pm_data, powerform, Post_PF)
      variable_ordering_auxiliary(pm)
      con(pm, pm.cnw)[:x] = Dict{Int,JuMP.ConstraintRef}()
      for l in pm.data["undecided_branches"]
        con(pm, pm.cnw)[:x][l] = constraint_abs_branch_flow_ordering(pm, l)
        JuMP.set_name(con(pm, pm.cnw)[:x][l],"x[$l]")  
      end
      objective_minmax_problem(pm)
      return pm
    else
        println("Not Supported Power Model Option. WARNING: returning 'nothing'")
	return nothing
    end
end

function SolveMinmax(pm_data,pm_form,optimizer)
  pm = MaximinOPF.MinimaxOPFModel(pm_data, pm_form)
  pm.data["undecided_branches"]= filter(l->!(l in pm_data["protected_branches"] || l in pm_data["inactive_branches"]), ids(pm,pm.cnw,:branch)) 
  model=pm.model

  SolveMinmax(model,optimizer)

  pm.data["x_vals"]=Dict{Int64,Float64}()
  for l in pm.data["undecided_branches"]
        pm.data["x_vals"][l]=JuMP.dual(con(pm, pm.cnw)[:x][l])
        pm.data["x_vals"][l] = min(1,pm.data["x_vals"][l])
        pm.data["x_vals"][l] = max(0,pm.data["x_vals"][l])
  end
  for l in pm.data["protected_branches"]
        pm.data["x_vals"][l] = 0
  end
  for l in pm.data["inactive_branches"]
        pm.data["x_vals"][l] = 1
  end
  return model, pm
end #end of function

function SolveMinmax(model::JuMP.Model,optimizer)
  JuMP.set_optimizer(model,optimizer)
  JuMP.optimize!(model)
  status=JuMP.termination_status(model)
  if status != OPTIMAL
    println("FLAGGING: Solve status=",status)
  end
end

function SolveFP(pm_data,powerform,optimizer, x_vals=Dict{Int64,Float64}() )
  pm = MaximinOPF.PF_FeasModel(pm_data, powerform, x_vals)
  JuMP.set_optimizer(pm.model,optimizer)
  JuMP.optimize!(pm.model)
  status=JuMP.termination_status(pm.model)
  if status != OPTIMAL
    println("FLAGGING: Solve status=",status)
  end
  return pm
end 

function PF_FeasModel(pm_data, powerform, x_vals=Dict{Int64,Float64}() )
    if powerform in supported_pm
        pm = instantiate_model(pm_data, powerform, Post_PF)
        for l in ids(pm,pm.cnw,:branch)
          if !haskey(x_vals,l)
	    if l in pm_data["inactive_branches"]
		x_vals[l]=1
	    else
		x_vals[l]=0
	    end
          end
        end
        objective_feasibility_problem(pm,x_vals)
        return pm
    else
	println("Do nothing")
	return nothing
    end
end


function SolveMinmaxDual(pm_data,pm_form,optimizer)
  pm = MaximinOPF.MinimaxOPFModel(pm_data, pm_form)
  pm.data["undecided_branches"]= filter(l->!(l in pm_data["protected_branches"] || l in pm_data["inactive_branches"]), ids(pm,pm.cnw,:branch)) 
  model = MaximinOPF.DualizeMinmaxModel(pm)  

  SolveMinmaxDual(model,optimizer)

  pm.data["x_vals"]=Dict{Int64,Float64}()
  for l in pm.data["undecided_branches"]
        x=variable_by_name(model,"x[$l]_1")
        pm.data["x_vals"][l]=JuMP.value(x)
        pm.data["x_vals"][l] = min(1,pm.data["x_vals"][l])
        pm.data["x_vals"][l] = max(0,pm.data["x_vals"][l])
  end
  for l in pm.data["protected_branches"]
        pm.data["x_vals"][l] = 0
  end
  for l in pm.data["inactive_branches"]
        pm.data["x_vals"][l] = 1
  end
  return model, pm, con_dict
end #end of function

function SolveMinmaxDual(model::JuMP.Model,optimizer)
  JuMP.set_optimizer(model,optimizer)
  JuMP.optimize!(model)
  status=JuMP.termination_status(model)
  if status != OPTIMAL
    println("FLAGGING: Solve status=",status)
  end
end

function DualizeMinmaxModel(minmax_model_pm::AbstractPowerModel)
    if typeof(minmax_model_pm) in conic_supported_pm 
	    dualizable_minmax_model = ConvertModelToDualizableForm(minmax_model_pm.model)
        AssignModelDefaultConstraintNames(dualizable_minmax_model)
        dualized_minmax_model = dualize(dualizable_minmax_model)
        #add_artificial_var_bds(dualized_minmax_model)  #Some solvers may require the feasible set to be bounded.
        return dualized_minmax_model
    else 
        println("DualizeMinmaxModel not supported for convex conic PowerModels, returning nothing.")
	    return nothing
    end
end


function DualizeMinmaxModel(dualizable_minmax_model::JuMP.Model)
    AssignModelDefaultConstraintNames(dualizable_minmax_model)
    dualized_minmax_model = dualize(dualizable_minmax_model)
end

function ConvertModelToDualizableForm(model::JuMP.Model)
        soc_model = MOI.Utilities.Model{Float64}()
        soc_bridged_model = MOI.Bridges.Constraint.QuadtoSOC{Float64}(soc_model)
        MOI.copy_to(soc_bridged_model,backend(model))    

        dualizable_model = JuMP.Model()
        bridged_model = MOI.Bridges.Constraint.Square{Float64}(backend(dualizable_model))
        MOI.copy_to(bridged_model,soc_model)    
	return dualizable_model
end

function AssignModelDefaultConstraintNames(model::JuMP.Model)
    type_str_map = Dict{DataType,String}(
	VariableRef=>"Var", 
	GenericAffExpr{Float64,VariableRef}=>"Expr", 
        Array{GenericAffExpr{Float64,VariableRef},1}=>"VecExpr", 
	MathOptInterface.EqualTo{Float64}=>"EQ",
	MathOptInterface.GreaterThan{Float64}=>"GE",
	MathOptInterface.LessThan{Float64}=>"LE",
	MathOptInterface.SecondOrderCone=>"SOC",
	MathOptInterface.RotatedSecondOrderCone=>"RSOC",
	MathOptInterface.PositiveSemidefiniteConeTriangle=>"PSD",
	MathOptInterface.PositiveSemidefiniteConeSquare=>"PSDSQ",
    )
        con_types=list_of_constraint_types(model)
        n_con_types=length(con_types)
	con_dict = Dict{Tuple{DataType,DataType},Dict{Int64,String}}() 

	for kk=1:n_con_types
	  con = all_constraints(model, con_types[kk][1], con_types[kk][2]) 
	  con_dict[ ( con_types[kk][1],  con_types[kk][2] )  ] = Dict{Int64,String}()
	    n_con=length(con)
	    for jj=1:n_con
		if length(JuMP.name(con[jj])) == 0
	  	  con_dict[ ( con_types[kk][1],  con_types[kk][2] )  ][jj] = string(type_str_map[con_types[kk][1]],"_",type_str_map[con_types[kk][2]],"[",jj,"]")
  		  JuMP.set_name(con[jj],con_dict[ ( con_types[kk][1] ,  con_types[kk][2] ) ][jj] )
		else
	  	  con_dict[(  con_types[kk][1],  con_types[kk][2] ) ][jj]  = JuMP.name(con[jj])
		end
	    end
	end
	return con_dict
end

function write_to_cbf(model,fn_base::String)
    JuMP.write_to_file( model, string(fn_base,".cbf"), format = MOI.FileFormats.FORMAT_CBF)
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
function write_to_cbf_scip(model,fn_base::String)
    model_psd = convertSOCtoPSD(model)
    add_artificial_var_bds(model_psd)
    add_psd_initial_cuts(model_psd)
    fname=string(fn_base,"_scip",".cbf")
    JuMP_write_to_file( model_psd, fname, format = MOI.FileFormats.FORMAT_CBF)
    ### CHANGING VER 3 to VER 2 
    fix_version=`sed -i -z 's/VER\n3/VER\n2/g' $fname`
    run(fix_version)
end

###Copied from current master version of JuMP, as of 7 Feb 2020. At some point, we can use the JuMP stable version.
function JuMP_write_to_file(
    model::Model,
    filename::String;
    format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_AUTOMATIC
)
    dest = MOI.FileFormats.Model(format = format, filename = filename)
    # We add a `full_bridge_optimizer` here because MOI.FileFormats models may not
    # support all constraint types in a JuMP model.
    bridged_dest = MOI.Bridges.full_bridge_optimizer(dest, Float64)
    MOI.copy_to(bridged_dest, backend(model))
    # `dest` will contain the underlying model, with constraints bridged if
    # necessary.
    MOI.write_to_file(dest, filename)
    return
end

function Post_PF(pm::AbstractPowerModel)
    pm.setting["output"]=Dict{String,Any}("duals"=>true)
    if !haskey(pm.data,"inactive_branches")
        pm.data["inactive_branches"] = []
    end
    if !haskey(pm.data,"protected_branches")
        pm.data["protected_branches"] = []
    end
    pm.data["all_branches"] = ids(pm,pm.cnw,:branch)
    pm.data["undecided_branches"] = filter(l->!(l in pm.data["protected_branches"] || l in pm.data["inactive_branches"]), ids(pm,pm.cnw,:branch))
    if typeof(pm) in [SOCBFPowerModel, SOCBFConicPowerModel] 
      build_opf_bf(pm)
    else
      build_opf(pm)
    end

#=
    #println("Variables: ",var(pm,pm.cnw,:w))
    w_var = var(pm,pm.cnw,:w)
    w_ids_aux = filter( iii->( has_lower_bound(w_var[iii]) || has_upper_bound(w_var[iii]) ), keys(w_var))
    #w_ids_aux = filter(iii->(has_lower_bound(w_var[iii])), keys(w_var))
    @variable(pm.model, u_w[ii in w_ids_aux] >= 0)
    for ii in w_ids_aux
        if has_lower_bound(w_var[ii])
            lb = lower_bound(w_var[ii])
            delete_lower_bound(w_var[ii])
            @constraint(pm.model, w_var[ii] - lb +  u_w[ii] >= 0)
        end
        if has_upper_bound(w_var[ii])
            ub = upper_bound(w_var[ii])
            delete_upper_bound(w_var[ii])
            @constraint(pm.model, w_var[ii] - ub -  u_w[ii] <= 0)
        end
    end

    #println("Variables: ",var(pm,pm.cnw,:wr))
    wr_var = var(pm,pm.cnw,:wr)
    wr_ids_aux = filter(iii->( has_lower_bound(wr_var[iii]) || has_upper_bound(wr_var[iii]) ), keys(wr_var))
    #wr_ids_aux = filter(iii->( has_lower_bound(wr_var[iii]) ), keys(wr_var))
    @variable(pm.model, u_wr[ii in wr_ids_aux] >= 0)
    for ii in wr_ids_aux
        if has_lower_bound(wr_var[ii])
            lb = lower_bound(wr_var[ii])
            delete_lower_bound(wr_var[ii])
            @constraint(pm.model, wr_var[ii] - lb +  u_wr[ii] >= 0)
        end
        if has_upper_bound(wr_var[ii])
            ub = upper_bound(wr_var[ii])
            delete_upper_bound(wr_var[ii])
            @constraint(pm.model, wr_var[ii] - ub -  u_wr[ii] <= 0)
        end
    end

    #println("Variables: ",var(pm,pm.cnw,:wi))
    wi_var = var(pm,pm.cnw,:wi)
    wi_ids_aux = filter( iii->( has_lower_bound(wi_var[iii]) || has_upper_bound(wi_var[iii]) ), keys(wi_var) )
    #wi_ids_aux = filter( iii->( has_lower_bound(wi_var[iii]) ), keys(wi_var) )
    @variable(pm.model, u_wi[ii in wi_ids_aux] >= 0)
    for ii in wi_ids_aux
        if has_lower_bound(wi_var[ii])
            lb = lower_bound(wi_var[ii])
            delete_lower_bound(wi_var[ii])
            @constraint(pm.model, wi_var[ii] - lb +  u_wi[ii] >= 0)
        end
        if has_upper_bound(wi_var[ii])
            ub = upper_bound(wi_var[ii])
            delete_upper_bound(wi_var[ii])
            @constraint(pm.model, wi_var[ii] - ub -  u_wi[ii] <= 0 )
        end
    end
=#

    # Add new variables, modify some existing constraints and add new ones also.
    variable_branch_flow_slacks(pm)

    con(pm, pm.cnw)[:abs_pflow_fr] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw)[:abs_pflow_to] = Dict{Int,JuMP.ConstraintRef}()
    if haskey( var(pm,pm.cnw),:q)
      con(pm, pm.cnw)[:abs_qflow_fr] = Dict{Int,JuMP.ConstraintRef}()
      con(pm, pm.cnw)[:abs_qflow_to] = Dict{Int,JuMP.ConstraintRef}()
    end

    for i in ids(pm,:bus)
	if haskey(sol(pm, pm.cnw, :bus, i),:lam_kcl_r) && !(typeof(pm) in [NFAPowerModel])
	  pbusref=sol(pm, pm.cnw, :bus, i)[:lam_kcl_r]
          JuMP.set_name(pbusref,string("p_bal[",i,"]"))
	  for a in ref(pm,pm.cnw, :bus_arcs, i)
	    up1m = var(pm,pm.cnw, :up_br1)[a,0]
	    up1p = var(pm,pm.cnw, :up_br1)[a,1]
	    JuMP.set_normalized_coefficient(pbusref,up1m,-1)
	    JuMP.set_normalized_coefficient(pbusref,up1p,1)
	  end
	end
	if haskey(sol(pm, pm.cnw, :bus, i),:lam_kcl_i) && !(typeof(pm) in [DCPPowerModel,DCMPPowerModel,DCPLLPowerModel,NFAPowerModel])
	  qbusref=sol(pm, pm.cnw, :bus, i)[:lam_kcl_i]
          JuMP.set_name(qbusref,string("q_bal[",i,"]"))
	  for a in ref(pm,pm.cnw, :bus_arcs, i)
	    uq1m = var(pm,pm.cnw, :uq_br1)[a,0]
	    uq1p = var(pm,pm.cnw, :uq_br1)[a,1]
	    JuMP.set_normalized_coefficient(qbusref,uq1m,-1)
	    JuMP.set_normalized_coefficient(qbusref,uq1p,1)
	  end
	end
	
    end
    for l in setdiff(ids(pm, :branch),pm.data["protected_branches"])
        ref_p1,ref_p3,ref_q1,ref_q3 = constraint_def_abs_flow_values(pm, l)
	if ref_p1 != nothing
          con(pm, pm.cnw)[:abs_pflow_fr][l] = ref_p1
          JuMP.set_name(ref_p1,"pi_f[$l]")  
	end
	if ref_p3 != nothing
          con(pm, pm.cnw)[:abs_pflow_to][l] = ref_p3
          JuMP.set_name(ref_p3,"pi_t[$l]")  
	end
	if ref_q1 != nothing
          con(pm, pm.cnw)[:abs_qflow_fr][l] = ref_q1
          JuMP.set_name(ref_q1,"phi_f[$l]")  
	end
	if ref_q3 != nothing
          con(pm, pm.cnw)[:abs_qflow_to][l] = ref_q3
          JuMP.set_name(ref_q3,"phi_t[$l]")  
	end
    end


end

end # module
