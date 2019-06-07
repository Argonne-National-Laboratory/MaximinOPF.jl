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


TOL = 1e-4  #feasibility tolerance for lazy constraints
useLocalCuts = false
verboseOut = false

######### READING ARGS #############
## Default arguments
CASE_NUM="30"
K=4
HEUR = 0  # 0 no heuristic, 1 heuristic1,  2 lambdaF=lambdaT and muF=muT  3 heuristic2,
FORM = 0; AC=0; SOCP=1; SOCPBds=2; DC=3; SDP=4; ACSOC=5

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

# define data
lines, buses, generators, baseMVA = opfdata.lines, opfdata.buses, opfdata.generators, opfdata.baseMVA
nbuses, nlines, ngens = length(buses), length(lines), length(generators)

# build a dictionary between buses ids and their indexes
busIdx = mapBusIdToIdx(buses)

# set up the fromLines and toLines for each bus
fromLines, toLines = mapLinesToBuses(buses, lines, busIdx)

# generators at each bus
BusGeners = mapGenersToBuses(buses, generators, busIdx)


# collect overall running time
wtime = 0.0

# collect number of cuts applied
ncuts = 0

bestAttack = zeros(nlines+6)
bestVoltages = zeros(nbuses,2)
attacks = Dict()
noAttacks = 0  # number of attacks discovered
noNoGoodCuts = 0 # number of no-good cuts employed
attacks[0] = zeros(nlines+6) # storing the attacks with last three entries storing attack values for different relaxations
isNewAttack = false
isNewIncumbent = false


dualValUB = 1e20

acYffR = opfdata.admittancesAC.YffR; acYffI = opfdata.admittancesAC.YffI;
acYttR = opfdata.admittancesAC.YttR; acYttI = opfdata.admittancesAC.YttI;
acYftR = opfdata.admittancesAC.YftR; acYftI = opfdata.admittancesAC.YftI;
acYtfR = opfdata.admittancesAC.YtfR; acYtfI = opfdata.admittancesAC.YtfI;
acYshR = opfdata.admittancesAC.YshR; acYshI = opfdata.admittancesAC.YshI;

dcYffR = opfdata.admittancesDC.YffR; dcYffI = opfdata.admittancesDC.YffI;
dcYttR = opfdata.admittancesDC.YttR; dcYttI = opfdata.admittancesDC.YttI;
dcYftR = opfdata.admittancesDC.YftR; dcYftI = opfdata.admittancesDC.YftI;
dcYtfR = opfdata.admittancesDC.YtfR; dcYtfI = opfdata.admittancesDC.YtfI;
dcYshR = opfdata.admittancesDC.YshR; dcYshI = opfdata.admittancesDC.YshI;

if FORM == DC
  YffR = dcYffR; YffI = dcYffI;
  YttR = dcYttR; YttI = dcYttI;
  YftR = dcYftR; YftI = dcYftI;
  YtfR = dcYtfR; YtfI = dcYtfI;
  YshR = dcYshR; YshI = dcYshI;
else
  YffR = acYffR; YffI = acYffI;
  YttR = acYttR; YttI = acYttI;
  YftR = acYftR; YftI = acYftI;
  YtfR = acYtfR; YtfI = acYtfI;
  YshR = acYshR; YshI = acYshI;
end



N = 1:nbuses
L = 1:nlines
G = 1:ngens




PD = zeros(nbuses)
QD = zeros(nbuses)
Wmin = zeros(nbuses)
Wmax = zeros(nbuses)
for i in N
	PD[i] = buses[i].Pd / baseMVA
	QD[i] = buses[i].Qd / baseMVA
	Wmin[i] = (buses[i].Vmin)^2
	Wmax[i] = (buses[i].Vmax)^2
end
Pmin = zeros(ngens)
Pmax = zeros(ngens)
Qmin = zeros(ngens)
Qmax = zeros(ngens)
for g in G
	Pmin[g] = generators[g].Pmin
	Pmax[g] = generators[g].Pmax
	Qmin[g] = generators[g].Qmin
	Qmax[g] = generators[g].Qmax
end
PminI = zeros(nbuses)
PmaxI = zeros(nbuses)
QminI = zeros(nbuses)
QmaxI = zeros(nbuses)
  
for i in N
    for g in BusGeners[i]
	PminI[i] += Pmin[g]
	PmaxI[i] += Pmax[g]
	QminI[i] += Qmin[g]
	QmaxI[i] += Qmax[g]
    end
end

include("BranchAndCutDSP.jl") ### Initialize the defender subproblems with power flow balance enforced


η0Val = 0
α_val = zeros(nbuses)
β_val = zeros(nbuses)
γ_val = zeros(nbuses)
δ_val = zeros(nbuses)
λF_val = zeros(nlines)
μF_val = zeros(nlines)
λT_val = zeros(nlines)
μT_val = zeros(nlines)
x_val = zeros(nlines)

# Collect the timing results
time_Eta0SP = 0.0
time_MP = 0.0
total_TimeMP = 0.0
time_root = -1 
nMPUpdates=0
avg_TimeMP=0

# some counters
numcalls_Eta0SP = 0

# dual function value
dualValLB = 0.0
println("Done with initial setup.")

finalXSoln = zeros(nlines+1)
bestUBVal=1e20
nNodes=0
incVal=0
runtime=0

function printX(finalXSoln)
  for l in L
   if finalXSoln[l] > 0.5 
	@printf(" %d",l)
   end
  end
end

  # Import the appropriate subproblem formulation
if FORM == AC
  #include("lECP.jl")
  include("ProxPtSDP.jl")
  testECP()
elseif FORM == SDP
  #The PSD constraints of mMP 
  include("DualAC.jl")
elseif FORM == SOCP
  include("DualSOCP.jl")
end

