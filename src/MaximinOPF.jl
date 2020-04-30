module MaximinOPF
using PowerModels
using Dualization
using MathOptInterface

include("Variables.jl")
include("Objectives.jl")
include("Constraints.jl")
include("utils.jl")
greet() = print("Hello World!")

supported_pm=[  
    ACPPowerModel, ACRPowerModel, ACTPowerModel, 
    SOCWRPowerModel, SOCWRConicPowerModel, 
    SOCBFPowerModel, SOCBFConicPowerModel, 
    QCRMPowerModel, QCLSPowerModel, 
    SDPWRMPowerModel, SparseSDPWRMPowerModel,
    DCPPowerModel, DCMPPowerModel, NFAPowerModel,
    DCPLLPowerModel, LPACCPowerModel 
]
conic_supported_pm=[ 
    SOCWRConicPowerModel, SOCBFConicPowerModel, 
	SDPWRMPowerModel, SparseSDPWRMPowerModel,
    DCPPowerModel, DCMPPowerModel, NFAPowerModel
]
sdp_pm=[ SDPWRMPowerModel, SparseSDPWRMPowerModel ]
nonconvex_pm=[ACPPowerModel, ACRPowerModel, ACTPowerModel]
not_supported_pm=[IVRPowerModel]


function MaximinOPFModel(pm_data, pm_form; enforce_int=true, io=devnull, rm_therm_line_lim=false)
    println(io,"Calling MaximinOPFModel() with pm_form: ",pm_form)
    minmax_pm = MinimaxOPFModel(pm_data, pm_form; rm_therm_line_lim=rm_therm_line_lim)

    #Test Out
    println(io,"name: ", pm_data["name"])
    println(io,"attacker_budget: ", pm_data["attacker_budget"])
    println(io,"inactive_branches: ", pm_data["inactive_branches"])
    println(io,"protected_branches: ", pm_data["protected_branches"])
    println(io, "Model:")
    println(io, minmax_pm.model)

    if pm_form in conic_supported_pm 
        maxmin_model = DualizeMinmaxModel(minmax_pm)
        if enforce_int
            for l in minmax_pm.data["undecided_branches"]
                if has_lower_bound(variable_by_name(maxmin_model,"x[$l]_1"))
                    delete_lower_bound(variable_by_name(maxmin_model,"x[$l]_1"))
                end
                if has_upper_bound(variable_by_name(maxmin_model,"x[$l]_1"))
                    delete_upper_bound(variable_by_name(maxmin_model,"x[$l]_1"))
                end
                JuMP.set_integer(variable_by_name(maxmin_model,"x[$l]_1"))
                JuMP.@constraint(maxmin_model, variable_by_name(maxmin_model,"x[$l]_1") >= 0)
                JuMP.@constraint(maxmin_model, variable_by_name(maxmin_model,"x[$l]_1") <= 1)
            end
        end
	    return maxmin_model
    else
        println(Base.stderr,"WARNING: Model type ",pm_form," is either nonconic and/or nonconvex and thus is not currently supported by function MaximinOPFModel().")
	    println(Base.stderr,"\tWARNING: function MaximinOPFModel() is returning the MINMAX model and NOT the MAXMIN model.")
	    return minmax_pm.model
    end
end

function MinimaxOPFModel(pm_data, pm_form; rm_therm_line_lim=false)
    if pm_form in supported_pm
        pm = instantiate_model(pm_data, pm_form, Post_PF)
        variable_ordering_auxiliary(pm)
        con(pm, pm.cnw)[:x] = Dict{Int,JuMP.ConstraintRef}()
        for l in pm.data["undecided_branches"]
            con(pm, pm.cnw)[:x][l] = constraint_abs_branch_flow_ordering(pm, l)
            JuMP.set_name(con(pm, pm.cnw)[:x][l],"x[$l]")  
        end
        objective_minmax_problem(pm)
        removeRedundantConstraints(pm)
        if rm_therm_line_lim
	        removeThermalLineLimits(pm)
        end
        return pm
    else
        println(Base.stderr,"Not Supported Power Model Option. WARNING: returning 'nothing'")
	    return nothing
    end
end

function PF_FeasModel(pm_data, pm_form; rm_therm_line_lim=false, x_vals::Dict{Int64,Float64}=Dict{Int64,Float64}() )
    if pm_form in supported_pm
        pm = instantiate_model(pm_data, pm_form, Post_PF)
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
        removeRedundantConstraints(pm)
        if rm_therm_line_lim
	        removeThermalLineLimits(pm)
        end
        return pm
    else
        println(Base.stderr,"PF_FeasModel: ", pm_form," Not Supported Power Model Option. WARNING: returning 'nothing'")
	    return nothing
    end
end

function DualizeMinmaxModel(minmax_model_pm::AbstractPowerModel)
    if typeof(minmax_model_pm) in conic_supported_pm 
	    dualizable_minmax_model = ConvertModelToDualizableForm(minmax_model_pm.model)
        AssignModelDefaultConstraintNames(dualizable_minmax_model)
        dualized_minmax_model = dualize(dualizable_minmax_model;dual_names = DualNames("", ""))
        return dualized_minmax_model
    else 
        println(Base.stderr,"DualizeMinmaxModel: Not supported for ",typeof(minmax_model_pm)," PowerModels, returning 'nothing'.")
	    return nothing
    end
end

function DualizeMinmaxModel(dualizable_minmax_model::JuMP.Model)
    AssignModelDefaultConstraintNames(dualizable_minmax_model)
    dualized_minmax_model = dualize(dualizable_minmax_model;dual_names = DualNames("", ""))
end


function AssignModelDefaultConstraintNames(model::JuMP.Model)
    all_vars = JuMP.all_variables(model)
    n_vars = length(all_vars)
    var_ids = collect(1:n_vars)
    artificial_lb_var_ids=[]
    artificial_ub_var_ids=[]
    for vv in var_ids
        var_name = JuMP.name(all_vars[vv])
        if JuMP.has_upper_bound(all_vars[vv])
  		    JuMP.set_name(UpperBoundRef(all_vars[vv]), string(var_name,"_UB"))
        end
        if JuMP.has_lower_bound(all_vars[vv])
  		    JuMP.set_name(LowerBoundRef(all_vars[vv]), string(var_name,"_LB"))
        end
    end
    
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
	            up1m = var(pm,pm.cnw, :pd_br)[a,0]
	            up1p = var(pm,pm.cnw, :pd_br)[a,1]
	            JuMP.set_normalized_coefficient(pbusref,up1m,-1)
	            JuMP.set_normalized_coefficient(pbusref,up1p,1)
	        end
	    end
	    if haskey(sol(pm, pm.cnw, :bus, i),:lam_kcl_i) && !(typeof(pm) in [DCPPowerModel,DCMPPowerModel,DCPLLPowerModel,NFAPowerModel])
	        qbusref=sol(pm, pm.cnw, :bus, i)[:lam_kcl_i]
            JuMP.set_name(qbusref,string("q_bal[",i,"]"))
	        for a in ref(pm,pm.cnw, :bus_arcs, i)
	            uq1m = var(pm,pm.cnw, :qd_br)[a,0]
	            uq1p = var(pm,pm.cnw, :qd_br)[a,1]
	            JuMP.set_normalized_coefficient(qbusref,uq1m,-1)
	            JuMP.set_normalized_coefficient(qbusref,uq1p,1)
	        end
	    end
    end
    #for l in setdiff(ids(pm, :branch),pm.data["protected_branches"])
    for l in ids(pm, :branch)
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

end # "end of module"
