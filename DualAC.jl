#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#


function solveNodeAC(ndata,solninfo,xDualVals)
# The dual problem 

  nThreads=1
  mMP = Model(solver=MosekSolver(MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=nThreads))
  println("Using ",nThreads," threads.")
  @variable(mMP, -1 <= α[i=N] <= 1)
  @variable(mMP, -1 <= β[i=N] <= 1)
  @variable(mMP, γp[i=N] >= 0)
  @variable(mMP, γm[i=N] >= 0)
  @constraint(mMP, [i=N], γp[i]+γm[i] <= 1)
  @variable(mMP, ζpUB[g=G] >= 0)
  @variable(mMP, ζpLB[g=G] >= 0)
  @variable(mMP, ζqUB[g=G] >= 0)
  @variable(mMP, ζqLB[g=G] >= 0)
  @variable(mMP, λF[l=L])
  @variable(mMP, λT[l=L])
  @variable(mMP, μF[l=L])
  @variable(mMP, μT[l=L])
  @variable(mMP, 0.0 <= x[l=L] <= 1.0)
  @constraint(mMP, sum(x[l] for l in L) <= K)
  print("Fixing the following x: ")
  for l in L
    if ndata[l]!=-1
       setlowerbound(x[l],ndata[l])
       setupperbound(x[l],ndata[l])
       print(" x[$l]=",ndata[l])
    end
  end 
  print("\n")
  
  for i in N
   for g in BusGeners[i]
    @constraint(mMP, 
	-α[i] + ζpUB[g] - ζpLB[g] == 0 )
    @constraint(mMP, 
	-β[i] + ζqUB[g] - ζqLB[g] == 0 ) 
   end
  end

  @expression(mMP, C[i=1:(2*nbuses),j=i:(2*nbuses)], 0)
  for i in N
	C[i,i] += γp[i] - γm[i] + α[i]*acYshR[i] - β[i]*acYshI[i] 
        C[nbuses+i,nbuses+i] += γp[i] - γm[i] + α[i]*acYshR[i] - β[i]*acYshI[i] 
    for l in fromLines[i]
      C[i,i] += λF[l]*acYffR[l] - μF[l]*acYffI[l]  
      C[nbuses+i,nbuses+i] += λF[l]*acYffR[l] - μF[l]*acYffI[l]  
    end
    for l in toLines[i]
      C[i,i] += λT[l]*acYttR[l] - μT[l]*acYttI[l]  
      C[nbuses+i,nbuses+i] += λT[l]*acYttR[l] - μT[l]*acYttI[l]  
    end
  end
  for l in L
      from=busIdx[lines[l].from]; to=busIdx[lines[l].to]
      C[from,nbuses+to] -= 0.5*(  λF[l]*acYftI[l] + μF[l]*acYftR[l] - λT[l]*acYtfI[l] - μT[l]*acYtfR[l]  )   
      C[to,nbuses+from] += 0.5*(  λF[l]*acYftI[l] + μF[l]*acYftR[l] - λT[l]*acYtfI[l] - μT[l]*acYtfR[l]  )   
      if from < to
        C[from,to] += 0.5*(λF[l]*acYftR[l]  - μF[l]*acYftI[l] + λT[l]*acYtfR[l] - μT[l]*acYtfI[l])   
        C[nbuses+from,nbuses+to] += 0.5*(λF[l]*acYftR[l]  - μF[l]*acYftI[l] + λT[l]*acYtfR[l] - μT[l]*acYtfI[l])   
      else
        C[to,from] += 0.5*(λF[l]*acYftR[l]  - μF[l]*acYftI[l] + λT[l]*acYtfR[l] - μT[l]*acYtfI[l])   
        C[nbuses+to,nbuses+from] += 0.5*(λF[l]*acYftR[l]  - μF[l]*acYftI[l] + λT[l]*acYtfR[l] - μT[l]*acYtfI[l])   
      end
  end
  @variable(mMP, H[i=1:(2*nbuses),j=1:(2*nbuses)], SDP)
  @constraint(mMP, SetH[i=1:(2*nbuses),j=i:(2*nbuses)], C[i,j] - H[i,j] == 0)

  @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γm[i]*Wmin[i]-γp[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N))


  @constraint(mMP, [l in L], α[busIdx[lines[l].from]] - x[l] - λF[l] <= 0)
  @constraint(mMP, [l in L], α[busIdx[lines[l].from]] + x[l] - λF[l] >= 0)
  @constraint(mMP, [l in L], -(1 - x[l]) - λF[l] <= 0)
  @constraint(mMP, [l in L], (1 - x[l]) - λF[l] >= 0)
  @constraint(mMP, [l in L], α[busIdx[lines[l].to]] - x[l] - λT[l] <= 0)
  @constraint(mMP, [l in L], α[busIdx[lines[l].to]] + x[l] - λT[l] >= 0)
  @constraint(mMP, [l in L], -(1 - x[l]) - λT[l] <= 0)
  @constraint(mMP, [l in L], (1 - x[l]) - λT[l] >= 0)

  @constraint(mMP, [l in L], β[busIdx[lines[l].from]] - x[l] - μF[l] <= 0)
  @constraint(mMP, [l in L], β[busIdx[lines[l].from]] + x[l] - μF[l] >= 0)
  @constraint(mMP, [l in L], -(1 - x[l]) - μF[l] <= 0)
  @constraint(mMP, [l in L], (1 - x[l]) - μF[l] >= 0)
  @constraint(mMP, [l in L], β[busIdx[lines[l].to]] - x[l] - μT[l] <= 0)
  @constraint(mMP, [l in L], β[busIdx[lines[l].to]] + x[l] - μT[l] >= 0)
  @constraint(mMP, [l in L], -(1 - x[l]) - μT[l] <= 0)
  @constraint(mMP, [l in L], (1 - x[l]) - μT[l] >= 0)

  #These constraints are not quite valid, but their inclusion often results in much faster time to near optimal solution.

  if HEUR == 1 
  	@constraint(mMP, LambdaMuConstr1[l in L], λF[l]*YftI[l] - λT[l]*YtfI[l] + μF[l]*YftR[l] - μT[l]*YtfR[l] == 0.0)
  elseif HEUR == 2
	@constraint(mMP, LambdaFequalsT[l in L], λF[l] - λT[l] == 0) 
	@constraint(mMP, muFequalsT[l in L], μF[l] - μT[l] == 0)
  elseif HEUR == 3
	@constraint(mMP, LambdaMuConstr2[l in L], λF[l]*YtfR[l] - λT[l]*YftR[l] - μF[l]*YtfI[l] + μT[l]*YftI[l] == 0.0)
  end

  status=solve(mMP)
  if status == :Optimal || status == :Stall
      for l in L
        solninfo[l] = getvalue(x[l])
        from=busIdx[lines[l].from]; to=busIdx[lines[l].to]
	xDualVals[l] = 1 + abs(getvalue(α[from]))+abs(getvalue(α[to]))+abs(getvalue(β[from]))+abs(getvalue(β[to]))
        #print(" d[$l]=",getdual(x[l]))
      end
      #print("\n")
      solninfo[nlines+1] = getobjectivevalue(mMP)
      solninfo[nlines+2] = getsolvetime(mMP)
      #println(" with optimal value ",solninfo[nlines+1])
      if status == :Stall
	println("solveNodeAC: Return status $status")
      end	
  else
	println("solveNodeAC: Return status $status")
  end
end #end of function

function xIntTol(x_val)
    tol = 1e-6
    for l in L
      if min(abs(x_val[l]),abs(1-x_val[l])) > tol
	return false
      end
    end
    # at this point, x_val is verified to be binary within tolerance
    for l in L
	x_val[l] = round(x_val[l])
    end
    return true
end
function findBranchIdx(x_val)
  maxidx=1
  maxval=min(x_val[1], 1-x_val[1])
  for l in L
    if min(x_val[l], 1-x_val[l]) > maxval
	maxidx=l
	maxval = min(x_val[l],1-x_val[l])
    end
  end
  return maxidx,maxval
end

absCoeff=zeros(nlines)
for l in L
#=
  absCoeff[l] = 
		  abs(acYffR[l]) + abs(acYffI[l]) 
	        + abs(acYttR[l]) + abs(acYttI[l]) 
	  	+ abs(acYftR[l]) + abs(acYftI[l]) + abs(acYtfR[l]) + abs(acYtfI[l]) 
	  	+ abs(acYftI[l]) + abs(acYtfI[l]) + abs(acYftR[l]) + abs(acYtfR[l])
=#
  absCoeff[l] = 
		  sqrt(acYffR[l]^2 + acYffI[l]^2 + acYttR[l]^2 + acYttI[l]^2
	  	+ acYftR[l]^2 + acYftI[l]^2 + acYtfR[l]^2 + acYtfI[l]^2)
end
function findBranchIdx2(x_val)
  global absCoeff
  maxidx=1
  maxval=absCoeff[1]*min(x_val[1], 1-x_val[1])
  for l in L
    if absCoeff[l]*min(x_val[l], 1-x_val[l]) > maxval
	maxidx=l
	maxval = absCoeff[l]*min(x_val[l],1-x_val[l])
    end
  end
  return maxidx,maxval
end
function findBranchIdx3(x_val,xDualVals)
  maxidx=1
  maxval=xDualVals[1]*min(x_val[1], 1-x_val[1])
  for l in L
    if xDualVals[l]*min(x_val[l], 1-x_val[l]) > maxval
	maxidx=l
	maxval = xDualVals[l]*min(x_val[l],1-x_val[l])
    end
  end
  return maxidx,maxval
end
function findNextNode(E)
  global nlines
  nodekey=-1
  weakestUBVal = -1
  for (k,n) in E
    if n[2][nlines+1] > weakestUBVal
	nodekey = n[1]
        weakestUBVal = n[2][nlines+1]
    end
  end
  return nodekey,weakestUBVal
end

function testDualAC()
  x_val=zeros(nlines)
  x_val[41]=1
  x_val[80]=1
  solninfo=zeros(nlines+2)
  solveNodeAC(x_val,solninfo)
  dualObjval = solninfo[nlines+1]
  primalObjval = solveFullModelSDP(x_val)
  println("Primal value: ",primalObjval," and dual value: ",dualObjval)
end

#Used for computing pseudocosts for branching
nIp=zeros(nlines)
nIm=zeros(nlines)
pIp=ones(nlines)
pIm=ones(nlines)
function solveBnBSDP(incSoln)
  global nIp,nIm, MAX_TIME
  start_time = time_ns()
  nodedata = Dict()
  nodedata[0] = -1*ones(nlines+1)
  nodedata[0][nlines+1]=0
  nodekey=0
  feasXs = Dict()
  nXs = 0
  incVal= -1
  bestUBVal=1e20
  nNodes=1
  solninfo = zeros(nlines+2)
  xDualVals = ones(nlines)
  maxidx=1
  maxval=-1
  E=enumerate(nodedata)
  while true
    println("Best UB $bestUBVal versus incumbent value $incVal")
    currNode = pop!(nodedata,nodekey)
    objval=currNode[nlines+1]
    if objval < incVal
        println("\tFathoming due to initial testing of bound $objval < $incVal")
    else
      solveNodeAC(currNode,solninfo,xDualVals)
      x_val = solninfo[1:nlines]
      objval = solninfo[nlines+1]
      soltime = solninfo[nlines+2]
      if objval < incVal
        println("\tFathoming due to bound $objval < $incVal")
      else
	# Apply primal heuristic
        if xIntTol(x_val)
          println("\tFathoming due to optimality")
          # Apply primal heuristic
	  incVal,nXs=primHeurXInt(x_val,objval,feasXs,nXs,incSoln,incVal)
        else
          #maxidx,maxval=findBranchIdx(x_val)
          maxidx,maxval=findBranchIdx2(x_val)
          #maxidx,maxval=findBranchIdx3(x_val,xDualVals)
          println("\tBranching on index $maxidx with value $maxval")
          nodedata[nNodes]=-1*ones(nlines+1)
          for l in L
            nodedata[nNodes][l]=currNode[l]
          end
          nodedata[nNodes][maxidx]=0
          nodedata[nNodes][nlines+1]=objval
          nNodes += 1
          nodedata[nNodes]=-1*ones(nlines+1)
          for l in L
            nodedata[nNodes][l]=currNode[l]
          end
          nodedata[nNodes][maxidx]=1
          nodedata[nNodes][nlines+1]=objval
          nNodes += 1
          # Apply primal heuristic
          if (nNodes % 10) == 0
	    incVal,nXs=primHeur(x_val,feasXs,nXs,incSoln,incVal,xDualVals)
          end
        end 
      end #not fathomed due to bound after solving node
    end # no intial fathoming
    E=enumerate(nodedata)
    println("There are ",length(E)," nodes left.") 
    if length(E) > 0
      nodekey,bestUBVal=findNextNode(E)
    else
      bestUBVal = incVal
      break
    end
    if (time_ns()-start_time)/1e9 > MAX_TIME
	break
    end
  end # while
  println("Best UB $bestUBVal versus incumbent value $incVal")
  end_time = time_ns()

  runtime = (end_time-start_time)/1e9
  return bestUBVal,nNodes,incVal,runtime
end

function primHeur(x_val,feasXs,nXs,incSoln,incVal,xDualVals)
  sortIdx = zeros(Int,nlines)
  sortperm!(sortIdx,x_val[1:nlines])
  bestKIdx = sortIdx[(nlines-K+1):nlines] 
  pX = zeros(Int,nlines)
  pX[bestKIdx]= 1
  Xs = enumerate(feasXs)
  isNewX = true
  for (k,feasx) in Xs
    if sum( abs(pX[l]-feasx[2][l]) for l=1:nlines ) < 1e-6
	isNewX=false
	break
    end
  end
  if isNewX
    sinfo = zeros(nlines+2)
    nXs += 1
    feasXs[nXs]=pX
    solveNodeAC(pX,sinfo,xDualVals)
    pVal = sinfo[nlines+1]
    if pVal > incVal
      incVal = pVal
      incSoln[1:nlines]=pX[1:nlines]
      incSoln[nlines+1]=incVal
      print("New incumbent solution: ")
      printX(pX)
      println("with value: $pVal")
    end
  end
  return incVal,nXs
end

function primHeurXInt(x_val,opt_val,feasXs,nXs,incSoln,incVal)
  for l in L
    x_val[l] = round(x_val[l])
  end
  Xs = enumerate(feasXs)
  isNewX = true
  for (k,feasx) in Xs
    if sum( abs(x_val[l]-feasx[2][l]) for l=1:nlines ) < 1e-6
	isNewX=false
	break
    end
  end
  if isNewX
    nXs += 1
    pX = x_val
    pVal = opt_val
    if pVal > incVal
      incVal = pVal
      incSoln[1:nlines]=pX[1:nlines]
      print("New incumbent solution: ")
      printX(x_val)
      print("with value: $pVal")
    end
  end
  return incVal,nXs
end


finalXSoln = zeros(Int,nlines)
bestUBVal,nNodes,incVal,runtime = solveBnBSDP(finalXSoln)

@printf("\n********************FINAL RESULT FOR CASE %s WITH %d LINES CUT H%dR%d*******************\n",CASE_NUM,K, HEUR,FORM)

@printf("Summary\t")
if(HEUR == 0 && FORM == AC)
    @printf("%d  &  ",K)
else
    @printf("    &   ")
end
if(FORM == AC )
    @printf("\$\\HSet^{%d}\$  &  ",HEUR)
else
    @printf("    &   ")
end

printResults=true
if printResults
if FORM == AC
    @printf("AC &    ")
elseif FORM == SDP  
    @printf("SDP &    ")
elseif FORM == SOCP || FORM == SOCPBds 
    @printf("SOCP &    ")
elseif FORM == DC 
    @printf("DC &    ")
elseif FORM == ACSOC 
    @printf("ACSOC &    ")
else
    println("Not implemented")
    printResults=false
end

#@printf("%d & ", ncuts)
#@printf("%.2f  & ", 100*time_Eta0SP/runtime)
@printf("%d &  ", nNodes)
@printf("%.2f &  ", runtime/nNodes)
#@printf("%d & ", ncuts)
# @printf("%d &   ", round(time_root-tShift))
#@printf("%d &  ", round(time_Eta0SP))
#@printf("%d &  ", round(total_TimeMP))
@printf("%d &  ", runtime)
printX(finalXSoln)
@printf(" & %.3f &  ", incVal)
if abs(bestUBVal-incVal) < 1e-3
  @printf("0.0\\%% &")
elseif incVal > 1e-3
  percentGap= 100*(bestUBVal-incVal)/incVal
  @printf("%.1f\\%% &",percentGap)
else
  print("---    & ")
end
  #@printf("\t&\t%.3f \t&\t%.3f\t & \t%.3f", optSDP, optSOCP, optDC)
@printf(" \\\\ \n")

println("No. Nodes: ", nNodes)
println("Best bound:  ", bestUBVal)
@printf("Objective value: %.3f\n", incVal)
@show runtime
@show finalXSoln 
end
