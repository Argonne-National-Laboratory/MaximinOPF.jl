#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

using JuMP
#using CPLEX, Ipopt, SCS 
using Ipopt 
#using Mosek,MosekTools
using Arpack
using DelimitedFiles
using Printf
using SparseArrays
#using MathProgBase
using LinearAlgebra
using ReverseDiffSparse


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

include("../src/EvalDSP.jl") ### Initialize the defender subproblems with power flow balance enforced




# Import the appropriate subproblem formulation
N, L, G = opfdata.N, opfdata.L, opfdata.G 
x_lbs=zeros(opfdata.nlines)
x_ubs=ones(opfdata.nlines)
node_data=NodeInfo(x_lbs,x_ubs,1e20)
### Params
  #ALG::Int
  #maxNSG::Int
  #tMin::Float64
  #tMax::Float64
  #tVal::Float64
  #rho::Float64
  #rhoUB::Float64
  #ssc::Float64
  #ssc_cntr::Int
  #tol::Float64
###
params=Params(3,100000,0.01,20,10.0,0.0,0.0,0.25,0,1e-5,1e-3,1e-3)


if FORM == ECP 
  include("../src/lECP.jl")
  solveLECP(opfdata,K,HEUR)
elseif FORM == ProxPtSDP
  #testLevelBM(opfdata,K,HEUR,node_data)
  #include("../src/PBM-SagastizabalSolodov.jl")
  #PBM_SagastizabalSolodov(opfdata,params,K,HEUR,node_data)
  #testProxPt(opfdata,K,HEUR,node_data)
  include("../src/PBM-DelfinoOliveira.jl")
  PBM_DelfinoOliveira(opfdata,params,K,HEUR,node_data)
  #testProxTraj(opfdata,K,HEUR,node_data)
elseif FORM == AC
  include("../src/DualAC.jl")
elseif FORM == SOCP
  include("../src/DualSOCP.jl")
else
end

