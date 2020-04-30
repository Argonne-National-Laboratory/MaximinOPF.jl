using MaximinOPF
using PowerModels
using JuMP
using SCS, Ipopt
using Test

include("testcases.jl")
PowerModels.silence()

function run_instance(testcase, PMOption, maximin = true, is_conic = true)
    pm_data = PowerModels.parse_file( testcase["file"] )
    pm_form = PMOption
    pm_data["name"]=testcase["name"]
    pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
    pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
    pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry
    pm_data["pm_form"] = pm_form

    #Create JUMP Model
    if maximin
        model = MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=false) # "SCS does not have mixed-integer support"
    else
        pm = MaximinOPF.MinimaxOPFModel(pm_data, pm_form)
        model = pm.model
    end

    if is_conic
        set_optimizer(model, optimizer_with_attributes(SCS.Optimizer, "verbose" => 0))
    else
        set_optimizer(model, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    end
    result = @elapsed JuMP.optimize!(model)
    status = JuMP.termination_status(model)
    println(" [$(testcase["name"])]: $(result) seconds, status $(status).")
    @test status in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

    # Note: We do not test writing files, because the functions are defined in MOI.

    # t_sdpa = @elapsed JuMP.write_to_file(model, "$(testcase["name"]).dat-s", format = MOI.FileFormats.FORMAT_SDPA)
    # println("  Wrote $(testcase["name"]).dat-s file: $(t_sdpa) seconds")
    # rm("$(testcase["name"]).dat-s")

    # t_cbf = @elapsed JuMP.write_to_file(model, "$(testcase["name"]).cbf", format = MOI.FileFormats.FORMAT_CBF)
    # println("  Wrote $(testcase["name"]).cbf file: $(t_cbf) seconds")
    # rm("$(testcase["name"]).cbf")
end

@testset "Maximin Model Tests" begin
    println("Maximin Models:")
    @testset "SOC Model Tests" begin
        println("Type: SOCWRConicPowerModel")
        for case in testcases
            run_instance(case, SOCWRConicPowerModel)
        end
    end
    @testset "SDP Model Tests" begin
        println("Type: SDPWRMPowerModel")
        for case in testcases
            run_instance(case, SDPWRMPowerModel)
        end
    end
    @testset "Sparse SDP Model Tests" begin
        println("Type: SparseSDPWRMPowerModel")
        for case in testcases
            run_instance(case, SparseSDPWRMPowerModel)
        end
    end
end

@testset "Minimax Model Tests" begin
    println("Minimax Models:")
    @testset "Nonconvex AC Models" begin
        @testset "ACP Model Tests" begin
            println("Type: ACPPowerModel")
            for case in testcases
                run_instance(case, ACPPowerModel, false, false)
            end
        end
        @testset "ACR Model Tests" begin
            println("Type: ACRPowerModel")
            for case in testcases
                run_instance(case, ACRPowerModel, false, false)
            end
        end
        @testset "ACT Model Tests" begin
            println("Type: ACTPowerModel")
            for case in testcases
                run_instance(case, ACTPowerModel, false, false)
            end
        end
    end
    @testset "Linear Approximation Models" begin
        @testset "DCP Model Tests" begin
            println("Type: DCPPowerModel")
            for case in testcases
                run_instance(case, DCPPowerModel, false, false)
            end
        end
        @testset "DCMP Model Tests" begin
            println("Type: DCMPPowerModel")
            for case in testcases
                run_instance(case, DCMPPowerModel, false, false)
            end
        end
        @testset "NFA Model Tests" begin
            println("Type: NFAPowerModel")
            for case in testcases
                run_instance(case, NFAPowerModel, false, false)
            end
        end
    end
    @testset "Quadratic Approximation Models" begin
        @testset "DCPLL Model Tests" begin
            println("Type: DCPLLPowerModel")
            for case in testcases
                run_instance(case, DCPLLPowerModel, false, false)
            end
        end
        @testset "LPACC Model Tests" begin
            println("Type: LPACCPowerModel")
            for case in testcases
                run_instance(case, LPACCPowerModel, false, false)
            end
        end
    end
    @testset "Quadratic Relaxation Models" begin
        @testset "SOCWR Model Tests" begin
            println("Type: SOCWRPowerModel")
            for case in testcases
                run_instance(case, SOCWRPowerModel, false, false)
            end
        end
        @testset "SOCBF Model Tests" begin
            println("Type: SOCBFPowerModel")
            for case in testcases
                run_instance(case, SOCBFPowerModel, false, false)
            end
        end
        @testset "QCRM Model Tests" begin
            println("Type: QCRMPowerModel")
            for case in testcases
                run_instance(case, QCRMPowerModel, false, false)
            end
        end
        @testset "QCLS Model Tests" begin
            println("Type: QCLSPowerModel")
            for case in testcases
                run_instance(case, QCLSPowerModel, false, false)
            end
        end
    end
    @testset "Quadratic Conic Models" begin
        @testset "SOCWR Model Tests" begin
            println("Type: SOCWRConicPowerModel")
            for case in testcases
                run_instance(case, SOCWRConicPowerModel, false)
            end
        end
        @testset "SOCBF Model Tests" begin
            println("Type: SOCBFConicPowerModel")
            for case in testcases
                run_instance(case, SOCBFConicPowerModel, false)
            end
        end
    end
    @testset "SDP Models" begin
        @testset "SDPWRM Model Tests" begin
            println("Type: SDPWRMPowerModel")
            for case in testcases
                run_instance(case, SDPWRMPowerModel, false)
            end
        end
        @testset "SparseSDPWRM Model Tests" begin
            println("Type: SparseSDPWRMPowerModel")
            for case in testcases
                run_instance(case, SparseSDPWRMPowerModel, false)
            end
        end
    end
end
