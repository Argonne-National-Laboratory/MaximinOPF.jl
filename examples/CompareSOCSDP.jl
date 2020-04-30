include("../../MaximinOPF/src/MaximinOPF.jl")
using PowerModels
using JuMP
using Mosek
using MosekTools
using Ipopt
PowerModels.silence()

global case_instance="30"
attack_budget=4
for aa in 1:length(ARGS)
    global case_instance
    if occursin("--case=",ARGS[aa])
            case_instance=ARGS[aa][(length("--case=")+1):length(ARGS[aa])]
            println("case being set to ",case_instance)
    else
            println(Base.stderr,"Argument ",ARGS[aa]," not recognized.")
    end
end

testcase = Dict(
    "file" => string("data/case",case_instance,".m"),
    "name" => string("case",case_instance,"K",attack_budget),
    "attack_budget" => attack_budget,
    "inactive_indices" => [],
    "protected_indices" => [],
)

start_time = time_ns()
io = open(string("case",case_instance,"_pmform_comp.txt"),"w")

pm_data = PowerModels.parse_file(testcase["file"])
branch_ids = MaximinOPF.getBranchIds(pm_data)
n_branches=length(branch_ids)
println(io,"Branches: ",branch_ids)

all_k1_attacks=Int64[]
all_k2_attacks=Tuple{Int64,Int64}[]
all_k3_attacks=Tuple{Int64,Int64,Int64}[]

#for k1 in branch_ids
for k1 in 1:n_branches
    push!(all_k1_attacks,k1)
    for k2 in (k1+1):length(branch_ids)
        push!(all_k2_attacks,(k1,k2))
        for k3 in (k2+1):length(branch_ids)
            push!(all_k3_attacks,(k1,k2,k3))
        end
    end
end
#println("All K=1 attacks: ",all_k1_attacks)
#println("All K=2 attacks: ",all_k2_attacks)
#println("All K=3 attacks: ",all_k3_attacks)


pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry
 

function carry_out_test(pm_data,aa,io)
    conic_solver=with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0)
    ip_solver=with_optimizer(Ipopt.Optimizer,
        print_level=0,
        linear_solver="ma57",
        max_iter=50000
    )
    ac_model=MaximinOPF.PF_FeasModel(pm_data, ACRPowerModel)
    JuMP.set_optimizer(ac_model.model,ip_solver)
    sdp_model=MaximinOPF.PF_FeasModel(pm_data, SparseSDPWRMPowerModel)
    JuMP.set_optimizer(sdp_model.model,conic_solver)
    qc_model=MaximinOPF.PF_FeasModel(pm_data, QCRMPowerModel)
    JuMP.set_optimizer(qc_model.model,ip_solver)
    soc_model=MaximinOPF.PF_FeasModel(pm_data, SOCWRConicPowerModel)
    JuMP.set_optimizer(soc_model.model,conic_solver)


    message_str=""
    JuMP.optimize!(ac_model.model)
    ac_val = round(JuMP.objective_value(ac_model.model);digits=3)
    ac_status = JuMP.termination_status(ac_model.model)
    if ac_val > 0
        message_str=string(message_str,"Line(s) ",aa,":\t")
        JuMP.optimize!(sdp_model.model)
        sdp_val = round(JuMP.objective_value(sdp_model.model);digits=3)
        sdp_status = JuMP.termination_status(sdp_model.model)
        JuMP.optimize!(qc_model.model)
        qc_val = round(JuMP.objective_value(qc_model.model);digits=3)
        qc_status = JuMP.termination_status(qc_model.model)
        soc_val = 0
        if sdp_val==0 && qc_val==0
            soc_val=0
            message_str=string(message_str,"\tac ",ac_val,"\t0=sdp=qc=soc","\t",ac_status)
        else
            if sdp_val > 0 && qc_val > 0
                if qc_val - sdp_val > 1e-3 
                    if ac_val - qc_val > 1e-3
                        message_str=string(message_str,"\tac=",ac_val,"\tqc=",qc_val,"\tsdp=",sdp_val)
                    else
                        message_str=string(message_str,"\tac=qc=",ac_val,"\tsdp=",sdp_val)
                    end
                elseif sdp_val - qc_val > 1e-3
                    if ac_val - sdp_val > 1e-3
                        message_str=string(message_str,"\tac=",ac_val,"\tsdp=",sdp_val,"\tqc=",qc_val)
                    else
                        message_str=string(message_str,"\tac=sdp=",ac_val,"\tqc=",qc_val)
                    end
                else
                    if ac_val - sdp_val > 1e-3
                        message_str=string(message_str,"\tac=",ac_val,"\tsdp=qc=",sdp_val)
                    else
                        message_str=string(message_str,"\tac=sdp=qc=",ac_val)
                    end
                end
                JuMP.optimize!(soc_model.model)
                soc_val = round(JuMP.objective_value(soc_model.model);digits=3)
                soc_status = JuMP.termination_status(soc_model.model)
                if soc_val > 0
                    if soc_val == min(sdp_val,qc_val)
                        message_str=string(message_str,"=soc")
                    else
                        message_str=string(message_str,"\tsoc=",soc_val)
                    end
                else
                    message_str=string(message_str,"\t0=soc")
                end
            elseif sdp_val > 0
                soc_val=0
                if ac_val - sdp_val > 1e-3
                    message_str=string(message_str,"\tac=",ac_val,"\tsdp=",sdp_val,"\t0=qc=soc","\t",ac_status)
                else
                    message_str=string(message_str,"\tac=sdp=",ac_val,"\t0=qc=soc","\t",ac_status)
                end
            elseif qc_val > 0
                soc_val=0
                if ac_val - qc_val > 1e-3
                    message_str=string(message_str,"\tac=",ac_val,"\tqc=",qc_val,"\t0=sdp=soc","\t",ac_status)
                else
                    message_str=string(message_str,"\tac=qc=",ac_val,"\t0=sdp=soc","\t",ac_status)
                end
            end
        end
        println(io,message_str)
        #println(io,"\tFor debugging: ac=",ac_val," sdp=",sdp_val," qc=",qc_val," soc=",soc_val)
    end
end


println(io,"Evaluating K=1 Attacks:")
for aa in all_k1_attacks
    pm_data["attacker_budget"] = 1
    pm_data["inactive_branches"] = [aa,]
    carry_out_test(pm_data,aa,io)
end

println(io,"Evaluating K=2 Attacks:")
for aa in all_k2_attacks
    pm_data["attacker_budget"] = 2
    pm_data["inactive_branches"] = [aa[1],aa[2]]
    carry_out_test(pm_data,aa,io)
end

println(io,"Evaluating K=3 Attacks:")
for aa in all_k3_attacks
    pm_data["attacker_budget"] = 3
    pm_data["inactive_branches"] = [aa[1],aa[2],aa[3]]
    carry_out_test(pm_data,aa,io)
end

end_time = time_ns()
runtime = (end_time-start_time)/1e9
println(io,"Runtime: ",runtime)

close(io)

