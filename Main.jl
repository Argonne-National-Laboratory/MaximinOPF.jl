#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

using JuMP
using CPLEX, Ipopt, SCS 
using Mosek
#using Arpack
#using Printf
using MathProgBase

include("opfdata.jl")
#CASE_NUM in {9,30,57,118,300,1354pegase,2869pegase}
MAX_TIME = 24*3600 # in seconds


TOL = 1e-8  #feasibility tolerance for lazy constraints
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
@printf("%.2f seconds\n", (time_ns()-start_load)/1e9)

include("EvalDSP.jl") ### Initialize the defender subproblems with power flow balance enforced



function printX(opfdata,x_soln)
  for l in opfdata.L
   if x_soln[l] > 0.5 
	@printf(" %d",l)
   end
  end
end

  # Import the appropriate subproblem formulation
if FORM == ECP 
  include("lECP.jl")
  solveLECP(opfdata,K,HEUR)
elseif FORM == ProxPtSDP
  include("ProxPtSDP.jl")
  #testLevelBM(opfdata,K,HEUR)
  #testProxPt(opfdata,K,HEUR)
  testProxPt0(opfdata,K,HEUR)
elseif FORM == AC
  include("DualAC.jl")
elseif FORM == SOCP
  include("DualSOCP.jl")
else
   x_val = zeros(opfdata.nlines)
   x_val .= 0; x_val[45]=1
   @show solveFullModelSDP(opfdata,x_val,true)
   x_val .= 0; x_val[63]=1; x_val[65]=1
   @show solveFullModelSDP(opfdata,x_val,true)
   x_val .= 0; x_val[38]=1; x_val[41]=1; x_val[80]=1
   @show solveFullModelSDP(opfdata,x_val,true)
   x_val .= 0; x_val[33]=1; x_val[41]=1; x_val[48]=1; x_val[80]=1
   @show solveFullModelSDP(opfdata,x_val,true)

   x_val .= 0; x_val[45]=1
   @show solveFullModelSOCP(opfdata,x_val)
   x_val .= 0; x_val[63]=1; x_val[65]=1
   @show solveFullModelSOCP(opfdata,x_val)
   x_val .= 0; x_val[38]=1; x_val[41]=1; x_val[80]=1
   @show solveFullModelSOCP(opfdata,x_val)
   x_val .= 0; x_val[33]=1; x_val[41]=1; x_val[48]=1; x_val[80]=1
   @show solveFullModelSOCP(opfdata,x_val)

   x_val .= 0; x_val[41]=1
   @show solveFullModelSDP(opfdata,x_val,false)
   x_val .= 0; x_val[41]=1; x_val[80]=1
   @show solveFullModelSDP(opfdata,x_val,false)
   x_val .= 0; x_val[33]=1; x_val[41]=1; x_val[80]=1
   @show solveFullModelSDP(opfdata,x_val,false)
   x_val .= 0; x_val[60]=1; x_val[65]=1; x_val[66]=1; x_val[72]=1
   @show solveFullModelSDP(opfdata,x_val,false)

   x_val .= 0; x_val[41]=1
   @show solveFullModelSDP(opfdata,x_val,true)
   x_val .= 0; x_val[41]=1; x_val[80]=1
   @show solveFullModelSDP(opfdata,x_val,true)
   x_val .= 0; x_val[33]=1; x_val[41]=1; x_val[80]=1
   @show solveFullModelSDP(opfdata,x_val,true)
   x_val .= 0; x_val[60]=1; x_val[65]=1; x_val[66]=1; x_val[72]=1
   @show solveFullModelSDP(opfdata,x_val,true)
end

