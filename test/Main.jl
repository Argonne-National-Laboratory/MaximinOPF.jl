#=
Template for branch-and-cut method

July 2018--present 2020
Kibaek Kim
Brian Dandurand
=#
using PowerModels

#using Libdl
#libomp=Libdl.dlopen("/sandbox/schanen/petsc/wgpu/lib/libpetsc.so.3.011.3",RTLD_GLOBAL)
using Libdl
libmkl2=Libdl.dlopen("/soft/com/packages/intel/19/u2/mkl/lib/intel64/libmkl_core.so",RTLD_GLOBAL)
libmkl2=Libdl.dlopen("/soft/com/packages/intel/19/u2/mkl/lib/intel64/libmkl_intel_thread.so",RTLD_GLOBAL)
libmkl3=Libdl.dlopen("/soft/com/packages/intel/19/u2/mkl/lib/intel64/libmkl_def.so",RTLD_GLOBAL)
libmkl1=Libdl.dlopen("/soft/com/packages/intel/19/u2/mkl/lib/intel64/libmkl_avx512.so",RTLD_GLOBAL)


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
  data_path = "data/case"
  fn_base=string(data_path,CASE_NUM)
  fn_mat=string(fn_base,".m")
  pm_data = PowerModels.parse_file(fn_mat)
  pm_data["attacker_budget"]=K
  pm_data["name"]=string(fn_base,"K",K)
print("finished loading the data after ",(time_ns()-start_load)/1e9," seconds.\n")

if FORM == ECP 
  #include("../src/lECP.jl")
  #solveECP(pm_data)
  include("../src/phiECP.jl")
  solvePhiECP(pm_data)
elseif FORM == AC
  include("../src/SolveMaxminViaBnB.jl")
  solveMaxminViaBnB(pm_data,SparseSDPWRMPowerModel)
elseif FORM == SOCP
  include("../src/DualSOCP.jl")
elseif FORM == DC
  include("../src/DualDC.jl")
  testDualDC(opfdata)
else
  println("Invalid option in main(): Doing nothing.")
end

