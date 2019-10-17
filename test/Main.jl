#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

#using Libdl
#libomp=Libdl.dlopen("/sandbox/schanen/petsc/wgpu/lib/libpetsc.so.3.011.3",RTLD_GLOBAL)
using Libdl
libmkl2=Libdl.dlopen("/soft/com/packages/intel/19/u2/mkl/lib/intel64/libmkl_core.so",RTLD_GLOBAL)
libmkl2=Libdl.dlopen("/soft/com/packages/intel/19/u2/mkl/lib/intel64/libmkl_intel_thread.so",RTLD_GLOBAL)
libmkl3=Libdl.dlopen("/soft/com/packages/intel/19/u2/mkl/lib/intel64/libmkl_def.so",RTLD_GLOBAL)
libmkl1=Libdl.dlopen("/soft/com/packages/intel/19/u2/mkl/lib/intel64/libmkl_avx512.so",RTLD_GLOBAL)


include("../src/opfdata.jl")
include("../src/typedefs.jl")
#CASE_NUM in {9,30,57,118,300,1354pegase,2869pegase}
MAX_TIME = 24*3600 # in seconds

#TOL = 1e-8  #feasibility tolerance for lazy constraints
useLocalCuts = false
verboseOut = false


######### READING ARGS #############
## Default arguments
CASE_NUM="30"
K=4
HEUR = 0  # 0 no heuristic, 1 heuristic1,  2 lambdaF=lambdaT and muF=muT  3 heuristic2,
FORM = 0; ECP=0; AC=1; SOCP=2; ProxPtSDP=3; DC=4; SDP=5; ACSOC=6

if length(ARGS) > 0
  CASE_NUM = ARGS[1]
  if length(ARGS) > 1
    K = parse(Int,ARGS[2])
    if length(ARGS) > 2
      HEUR = parse(Int,ARGS[3])
      if length(ARGS) > 3
        FORM = parse(Int,ARGS[4])
        if length(ARGS) > 4
          MAX_TIME = parse(Int,ARGS[4])
	end
      end
    end
  end
else
  println("Usage: julia BranchAndCut.jl CASE_NUM K HEUR FORM")
  exit()
end
############# DONE READING ARGS ###############

print("Loading data ... "); start_load = time_ns()
# Load the bus system topology
  opfdata = opf_loaddata(CASE_NUM)
print("finished loading the data after ",(time_ns()-start_load)/1e9," seconds.\n")

#include("../src/EvalDSP.jl") ### Initialize the defender subproblems with power flow balance enforced




# Import the appropriate subproblem formulation
N, L, G = 1:opfdata.nbuses, 1:opfdata.nlines, 1:opfdata.ngens 
node_data=create_node(opfdata)
#= Params
  ALG::Int
  maxNIter::Int
  maxNSG::Int
  tMin::Float64
  tMax::Float64
  tStart::Float64
  ssc1::Float64
  ssc2::Float64
  ssc3::Float64
  tol1::Float64
  tol2::Float64
  tol3::Float64
=#


if FORM == ECP 
  include("../src/lECP.jl")
  solveLECP(opfdata,K,HEUR)
elseif FORM == ProxPtSDP
  #testLevelBM(opfdata,K,HEUR,node_data)
  #include("../src/PBM-SagastizabalSolodov.jl")
  #PBM_SagastizabalSolodov(opfdata,params,K,HEUR,node_data)
  #testProxPt(opfdata,K,HEUR,node_data)
  #testProxTraj(opfdata,K,HEUR,node_data)

  #include("../src/PBM-DelfinoOliveira.jl")
  #params=Params(3,10000,100,1e-2,1000.0,1.0,0.1,0.5,0.5,1e-4,1e-2,1e-2)
  #plot_params,plot_opt,plot_data,plot_params_ss,plot_opt_ss,plot_data_ss=PBM_DelfinoOliveira(opfdata,params,K,HEUR,node_data)
  #exptype="pbm"

  include("../src/PBM-SagastizabalSolodov.jl")
  params=Params(3,10000,100,1e-4,1000.0,1.0,0.1,0.5,0.5,1e-6,1e-5,1e-5)
  PBM_SagastizabalSolodov(opfdata,params,K,HEUR,node_data)

  #include("../src/CPAlg.jl")
  #plot_data=CPAlg(opfdata,params,K,HEUR,node_data)
  #exptype="cp"

  #include("../src/scipSDP.jl")
  #createSCIPModel(opfdata)

write_data = false
if write_data
  n_data,n_data_ss=size(plot_data)[1],size(plot_data_ss)[1]
  fn_base=string("ExpOut/",exptype,"exp",CASE_NUM,"-",K,"-",params.maxNSG)
  fn_base_ss=string("ExpOut/",exptype,"exp",CASE_NUM,"-",K,"-",params.maxNSG,"ss")

  fn_params=string(fn_base,"-params.dat")
  io = open(fn_params,"w")
  for kk=1:n_data
    write(io,string(kk," ",plot_params[kk,1]," ",plot_params[kk,2]," ",plot_params[kk,3],"\n"))
  end
  close(io)
  fn_params_ss=string(fn_base_ss,"-params.dat")
  io = open(fn_params_ss,"w")
  for kk=1:n_data_ss
    write(io,string(plot_params_ss[kk,4]," ",plot_params_ss[kk,1]," ",plot_params_ss[kk,2]," ",plot_params_ss[kk,3],"\n"))
  end
  close(io)

  fn_obj=string(fn_base,"-obj.dat")
  io = open(fn_obj,"w")
  for kk=1:n_data
    write(io,string(kk," ",plot_opt[kk,1],"\n"))
  end
  close(io)

  fn_obj_ss=string(fn_base_ss,"-obj.dat")
  io = open(fn_obj_ss,"w")
  for kk=1:n_data_ss
    write(io,string(plot_opt_ss[kk,6]," ",plot_opt_ss[kk,1],"\n"))
  end
  close(io)

  fn_eta=string(fn_base,"-eta.dat")
  io = open(fn_eta,"w")
  for kk=1:n_data
    write(io,string(kk," ",plot_opt[kk,2],"\n"))
  end
  close(io)

  fn_eta_ss=string(fn_base_ss,"-eta.dat")
  io = open(fn_eta_ss,"w")
  for kk=1:n_data_ss
    write(io,string(plot_opt_ss[kk,6]," ",plot_opt_ss[kk,2],"\n"))
  end
  close(io)

  fn_dual=string(fn_base,"-dual.dat")
  io = open(fn_dual,"w")
  for kk=1:n_data
    write(io,string(kk," ",max(plot_opt[kk,4],plot_opt[kk,5]),"\n"))
  end
  close(io)
  fn_dual_ss=string(fn_base_ss,"-dual.dat")
  io = open(fn_dual_ss,"w")
  for kk=1:n_data_ss
    write(io,string(plot_opt_ss[kk,6]," ",max(plot_opt_ss[kk,4],plot_opt[kk,5]),"\n"))
  end
  close(io)

  fn_init=string(fn_base,"-init.dat")
  io = open(fn_init,"w")
  for kk=1:n_data
    write(io,string(kk," ",plot_data[kk,3],"\n"))
  end
  close(io)
  fn_init_ss=string(fn_base_ss,"-init.dat")
  io = open(fn_init_ss,"w")
  for kk=1:n_data_ss
    write(io,string(plot_data_ss[kk,8]," ",plot_data_ss[kk,3],"\n"))
  end
  close(io)

  fn_solve=string(fn_base,"-solve.dat")
  io = open(fn_solve,"w")
  for kk=1:n_data
    write(io,string(kk," ",plot_data[kk,4],"\n"))
  end
  close(io)
  fn_solve_ss=string(fn_base_ss,"-solve.dat")
  io = open(fn_solve_ss,"w")
  for kk=1:n_data_ss
    write(io,string(plot_data_ss[kk,8]," ",plot_data_ss[kk,4],"\n"))
  end
  close(io)

  fn_sg=string(fn_base,"-sg.dat")
  io = open(fn_sg,"w")
  for kk=1:n_data
    write(io,string(kk," ",plot_data[kk,5],"\n"))
  end
  close(io)
  fn_sg_ss=string(fn_base_ss,"-sg.dat")
  io = open(fn_sg_ss,"w")
  for kk=1:n_data_ss
    write(io,string(plot_data_ss[kk,8]," ",plot_data_ss[kk,5],"\n"))
  end
  close(io)

  fn_pp=string(fn_base,"-pp.dat")
  io = open(fn_pp,"w")
  for kk=1:n_data
    write(io,string(kk," ",plot_data[kk,6],"\n"))
  end
  close(io)
  fn_pp_ss=string(fn_base_ss,"-pp.dat")
  io = open(fn_pp_ss,"w")
  for kk=1:n_data_ss
    write(io,string(plot_data_ss[kk,8]," ",plot_data_ss[kk,6],"\n"))
  end
  close(io)

  fn_bundle=string(fn_base,"-bundle.dat")
  io = open(fn_bundle,"w")
  for kk=1:n_data
    write(io,string(kk," ",plot_data[kk,7],"\n"))
  end
  close(io)
  fn_bundle_ss=string(fn_base_ss,"-bundle.dat")
  io = open(fn_bundle_ss,"w")
  for kk=1:n_data_ss
    write(io,string(plot_data_ss[kk,8]," ",plot_data_ss[kk,7],"\n"))
  end
  close(io)
end

elseif FORM == AC
  include("../src/DualAC.jl")

  #testSCSonRoot(opfdata)

  finalXSoln=create_node(opfdata)
  bestUBVal,nNodes,runtime = solveBnBSDP(opfdata,finalXSoln)
  println("No. Nodes: ", nNodes)
  println("Best bound:  ", bestUBVal)
  @printf("Objective value: %.3f\n", finalXSoln.nodeBd)
  @show runtime
  printX(opfdata,finalXSoln.x_soln)
  print("\n")
  @show finalXSoln.nodeBd
elseif FORM == SOCP
  include("../src/DualSOCP.jl")
else
end

