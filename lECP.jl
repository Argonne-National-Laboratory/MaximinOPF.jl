#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

function solveLECP(opfdata,K,HEUR)

  # define data
  lines, buses, generators, baseMVA = opfdata.lines, opfdata.buses, opfdata.generators, opfdata.baseMVA
  nbuses, nlines, ngens = length(buses), length(lines), length(generators)
  N = 1:nbuses; L = 1:nlines; G = 1:ngens
  # build a dictionary between buses ids and their indexes
  busIdx = mapBusIdToIdx(buses)
  # set up the fromLines and toLines for each bus
  fromLines, toLines = mapLinesToBuses(buses, lines, busIdx)
  fromBus=zeros(Int,nlines); toBus=zeros(Int,nlines)
  for l in L
    fromBus[l] = busIdx[lines[l].from]; toBus[l] = busIdx[lines[l].to] 
  end
  # generators at each bus
  BusGeners = mapGenersToBuses(buses, generators, busIdx)

  Y = Dict()  # Admittances
  Y["ffR"] = opfdata.admittancesAC.YffR; Y["ffI"] = opfdata.admittancesAC.YffI;
  Y["ttR"] = opfdata.admittancesAC.YttR; Y["ttI"] = opfdata.admittancesAC.YttI;
  Y["ftR"] = opfdata.admittancesAC.YftR; Y["ftI"] = opfdata.admittancesAC.YftI;
  Y["tfR"] = opfdata.admittancesAC.YtfR; Y["tfI"] = opfdata.admittancesAC.YtfI;
  Y["shR"] = opfdata.admittancesAC.YshR; Y["shI"] = opfdata.admittancesAC.YshI;

  PD = zeros(nbuses); QD = zeros(nbuses)
  Wmin = zeros(nbuses); Wmax = zeros(nbuses)
  for i in N
	PD[i] = buses[i].Pd / baseMVA; QD[i] = buses[i].Qd / baseMVA
	Wmin[i] = (buses[i].Vmin)^2; Wmax[i] = (buses[i].Vmax)^2
  end
  Pmin = zeros(ngens); Pmax = zeros(ngens)
  Qmin = zeros(ngens); Qmax = zeros(ngens)
  for g in G
	Pmin[g] = generators[g].Pmin; Pmax[g] = generators[g].Pmax
	Qmin[g] = generators[g].Qmin; Qmax[g] = generators[g].Qmax
  end

  println("Done with initial setup.")

  # The master problem MP
  mMP = Model(solver=CplexSolver(CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50))
  #mMP = Model(solver=CplexSolver(CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50,CPX_PARAM_LPMETHOD=4,CPX_PARAM_SOLUTIONTYPE=2,CPX_PARAM_STARTALG=4))
  # Define the model here
  @variable(mMP, x[l=L], Bin, start=0)
  @variable(mMP, -1 <= α[i=N] <= 1, start=0); @variable(mMP, -1 <= β[i=N] <= 1, start=0)
  @variable(mMP, δ[i=N] >= 0); @variable(mMP, γ[i=N] >= 0); @constraint(mMP, [i=N], δ[i]+γ[i] <= 1)

  @variable(mMP, λF[l=L], start=0); @variable(mMP, λT[l=L], start=0); @variable(mMP, μF[l=L], start=0); @variable(mMP, μT[l=L], start=0)

  @variable(mMP, ζpUB[g=G] >=0); @variable(mMP, ζpLB[g=G] >=0); @variable(mMP, ζqUB[g=G] >=0); @variable(mMP, ζqLB[g=G] >=0)

  for i in N
   for g in BusGeners[i]
    @constraint(mMP, -α[i] + ζpUB[g] - ζpLB[g] == 0 )
    @constraint(mMP, -β[i] + ζqUB[g] - ζqLB[g] == 0 ) 
   end
  end
#=
  if HEUR == 2
    @variable(mMP, xiPLB[l=L] >= 0); @variable(mMP, xiPUB[l=L] >= 0); @variable(mMP, xiQLB[l=L] >= 0); @variable(mMP, xiQUB[l=L] >= 0)
    hBd = 0.05
    @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γ[i]*Wmin[i]-δ[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N) - hBd*sum(xiPLB[l]+xiPUB[l]+xiQLB[l]+xiQUB[l]     for l in L)   )
  else
    @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γ[i]*Wmin[i]-δ[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N))
  end
=#
    @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γ[i]*Wmin[i]-δ[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N))

  @constraint(mMP, sum(x[l] for l in L) <= K)

 # McCormick inequalities enforcing bilinear equalities 
  @constraint(mMP, AMcf1[l in L], α[fromBus[l]] - x[l] <= λF[l]); @constraint(mMP, AMcf2[l in L], α[fromBus[l]] + x[l] >= λF[l])
  @constraint(mMP, AMcf3[l in L], -(1 - x[l]) <= λF[l]); @constraint(mMP, AMcf4[l in L], (1 - x[l]) >= λF[l])
  @constraint(mMP, AMct1[l in L], α[toBus[l]] - x[l] <= λT[l]); @constraint(mMP, AMct2[l in L], α[toBus[l]] + x[l] >= λT[l])
  @constraint(mMP, AMct3[l in L], -(1 - x[l]) <= λT[l]); @constraint(mMP, AMct4[l in L], (1 - x[l]) >= λT[l])

  @constraint(mMP, BMcf1[l in L], β[fromBus[l]] - x[l] <= μF[l]); @constraint(mMP, BMcf2[l in L], β[fromBus[l]] + x[l] >= μF[l])
  @constraint(mMP, BMcf3[l in L], -(1 - x[l]) <= μF[l]); @constraint(mMP, BMcf4[l in L], (1 - x[l]) >= μF[l])
  @constraint(mMP, BMct1[l in L], β[toBus[l]] - x[l] <= μT[l]); @constraint(mMP, BMct2[l in L], β[toBus[l]] + x[l] >= μT[l])
  @constraint(mMP, BMct3[l in L], -(1 - x[l]) <= μT[l]); @constraint(mMP, BMct4[l in L], (1 - x[l]) >= μT[l])

  #These constraints are not quite valid, but their inclusion often results in much faster time to near optimal solution.

  if HEUR == 1 
  	@constraint(mMP, LambdaMuConstr1[l in L], λF[l]*Y["ftI"][l] - λT[l]*Y["tfI"][l] + μF[l]*Y["ftR"][l] - μT[l]*Y["tfR"][l] == 0.0)
  elseif HEUR == 2
#=
	@constraint(mMP, LambdaFequalsT[l in L], λF[l] - λT[l] + xiPUB[l] - xiPLB[l] == 0) 
	@constraint(mMP, muFequalsT[l in L], μF[l] - μT[l] + xiQUB[l] - xiQLB[l] == 0)
=#
	@constraint(mMP, LambdaFequalsT[l in L], λF[l] - λT[l]  == 0) 
	@constraint(mMP, muFequalsT[l in L], μF[l] - μT[l]  == 0)
  elseif HEUR == 3
	#@constraint(mMP, LambdaMuConstr2[l in L], λF[l]*Y["tfR"][l] - λT[l]*Y["ftR"][l] - μF[l]*Y["tfI"][l] + μT[l]*Y["ftI"][l] == 0.0)
	for l in L
	 if (l<34 || l>38) && l!=44 && l!=46 && l!=69 && l!=70
	  @constraint(mMP, λF[l] - λT[l]  == 0) 
	  @constraint(mMP, μF[l] - μT[l]  == 0)
	 end
	end
  end


  maxNSG = 1
  sg = Dict()
  sg["α"] = zeros(nbuses,maxNSG)
  sg["β"] = zeros(nbuses,maxNSG)
  sg["γ"] = zeros(nbuses,maxNSG) 
  sg["δ"] = zeros(nbuses,maxNSG) 
  sg["λF"] = zeros(nlines,maxNSG)
  sg["μF"] = zeros(nlines,maxNSG)
  sg["λT"] = zeros(nlines,maxNSG)
  sg["μT"] = zeros(nlines,maxNSG)

  η0Val = 0

  x_val = zeros(nlines)

  pi_val = Dict()
  pi_val["α"] = zeros(nbuses)
  pi_val["β"] = zeros(nbuses)
  pi_val["γ"] = zeros(nbuses)
  pi_val["δ"] = zeros(nbuses)
  pi_val["λF"] = zeros(nlines)
  pi_val["λT"] = zeros(nlines)
  pi_val["μF"] = zeros(nlines)
  pi_val["μT"] = zeros(nlines)
  v = Dict()
  v["R"] = zeros(nbuses)
  v["I"] = zeros(nbuses)

  function addMPCutsLazyCB(cb)

	for l in L
    	   x_val[l] = round(getvalue(x[l]))
	end

	# generate cuts
	time_Eta0Start = time_ns()
	#  generate cut(s)
        for i in N
	  pi_val["α"][i] = getvalue(α[i]); pi_val["β"][i] = getvalue(β[i]); pi_val["γ"][i] = getvalue(γ[i]); pi_val["δ"][i] = getvalue(δ[i])
        end
	for l in L
	  pi_val["λF"][l] = getvalue(λF[l])
	  pi_val["λT"][l] = getvalue(λT[l])
	  pi_val["μF"][l] = getvalue(μF[l])
	  pi_val["μT"][l] = getvalue(μT[l])
	end
        H=spzeros(2*nbuses,2*nbuses)
	updateHess(Y,opfdata,pi_val,H)
  	try
    	  η0Val = solveEta0Eigs(H,opfdata,v)
  	catch exc
    	  println("Exception caught with eigs(), solving η0Val subproblem with Ipopt as recourse.") 
    	  println(exc)
    	  η0Val = solveEta0SDP(H,opfdata,v)
  	end
  	if η0Val < -TOL 
	    computeSG(Y,opfdata,v,sg)
	    @lazyconstraint(cb, 0.0 <= sum( (sg["α"][i,1])* α[i] for i in N) + sum( (sg["β"][i,1])* β[i] for i in N )  
		+ sum( (sg["γ"][i,1])* γ[i] for i in N)  + sum( (sg["δ"][i,1])* δ[i] for i in N)
		+ sum( (sg["λF"][l,1])* λF[l] + (sg["λT"][l,1])* λT[l] for l in L) + sum( (sg["μF"][l,1])* μF[l] + (sg["μT"][l,1])* μT[l] for l in L),
		localcut=useLocalCuts)
	    ncuts += 1
  	else
    	  println("Tolerance met for not generating a new lazy cut.")
  	end
	time_Eta0SP += (time_ns()-time_Eta0Start)/1e9
	numcalls_Eta0SP += 1
  end

  addlazycallback(mMP, addMPCutsLazyCB)
  ncuts = 0   # collect number of cuts applied
  mpRootStartTime = time_ns() #For solving root node
  # Collect the timing results
  time_Eta0SP = 0.0
  # some counters
  numcalls_Eta0SP = 0
  nMPUpdates=0
  wtime = @elapsed status = solve(mMP)
  @show wtime
  @printf("LAST ATTACK && ")
  for l in L
	if round(getvalue(x[l])) == 1
		println("Cut line $l")
	end
  end 
  @printf("\n")
  @printf("\n********************FINAL RESULT FOR CASE %s WITH %d LINES CUT H%dR%d*******************\n",CASE_NUM,K, HEUR,FORM)


  ###Printing summary as one line to aid in latex writeup
  @printf("Summary\t")
  if(HEUR == 0 )
    @printf("%d  &  ",K)
  else
    @printf("    &   ")
  end
  @printf("\$\\HSet^{%d}\$  &  ",HEUR)

  @printf("ECP &    ")


  @printf("%d  & ", ncuts)
  @printf("%.2f  & ", 100*time_Eta0SP/getsolvetime(mMP))
  @printf("%d   &  ", getnodecount(mMP))
  @printf("%.2f   &  ", getsolvetime(mMP)/getnodecount(mMP))
  @printf("%d   &  ", getsolvetime(mMP))
 #@printf("%.3f   &  ", getobjectivebound(mMP))

 #@printf("%d    &   ", round(time_root-tShift))
 #@printf("%d   &  ", round(time_Eta0SP))
 #@printf("%d   &  ", round(total_TimeMP))
  for l in L
   if getvalue(x[l]) > 0.5
        @printf(" %d",l)
   end
  end
  #@printf("\t&\t%.3f \t&\t%.3f\t & \t%.3f", bestAttack[nlines+1+SDP], bestAttack[nlines+1+SOCP],bestAttack[nlines+1+DC])
 @printf("  & %.3f ", getobjectivevalue(mMP))
 if abs(getobjectivebound(mMP)-getobjectivevalue(mMP)) < 1e-3
  @printf("&  0.0\\%% ")
 elseif incVal > 1e-3
  percentGap= 100*(getobjectivebound(mMP)-getobjectivevalue(mMP))/getobjectivevalue(mMP)
  @printf("&  %.1f\\%% ",percentGap)
 else
  print("& ---     ")
 end
 @printf(" \\\\ \n")
   return mMP
 end

# Update Hessian
function updateHess(Y,opfdata,pi_val,H)
    lines, buses, generators, baseMVA = opfdata.lines, opfdata.buses, opfdata.generators, opfdata.baseMVA
    nbuses, nlines, ngens = length(buses), length(lines), length(generators)
    N = 1:nbuses; L = 1:nlines; G = 1:ngens
    busIdx = mapBusIdToIdx(buses)
    fromBus=zeros(Int,nlines); toBus=zeros(Int,nlines)
    for l in L
      fromBus[l] = busIdx[lines[l].from]; toBus[l] = busIdx[lines[l].to] 
    end
    for i in N
       H[i,i] +=  pi_val["α"][i] * Y["shR"][i] - pi_val["β"][i] * Y["shI"][i]  + pi_val["δ"][i] - pi_val["γ"][i]
       H[nbuses+i,nbuses+i] += pi_val["α"][i] * Y["shR"][i] - pi_val["β"][i] * Y["shI"][i] + pi_val["δ"][i] - pi_val["γ"][i]
    end
    for l in L
       from = fromBus[l]; to = toBus[l] 
       H[from,from] += pi_val["λF"][l] * Y["ffR"][l] - pi_val["μF"][l] * Y["ffI"][l]
       H[nbuses+from,nbuses+from] += pi_val["λF"][l] * Y["ffR"][l] - pi_val["μF"][l] * Y["ffI"][l]
       H[to,to] += pi_val["λT"][l] * Y["ttR"][l] - pi_val["μT"][l] * Y["ttI"][l]
       H[nbuses+to,nbuses+to] += pi_val["λT"][l] * Y["ttR"][l] - pi_val["μT"][l] * Y["ttI"][l]
       H[from,to] += 0.5*( pi_val["λF"][l] * Y["ftR"][l] - pi_val["μF"][l] * Y["ftI"][l] + pi_val["λT"][l] * Y["tfR"][l] - pi_val["μT"][l] * Y["tfI"][l] )
       H[to,from] += 0.5*( pi_val["λF"][l] * Y["ftR"][l] - pi_val["μF"][l] * Y["ftI"][l] + pi_val["λT"][l] * Y["tfR"][l] - pi_val["μT"][l] * Y["tfI"][l] )
       H[nbuses+from, nbuses+to] += 0.5*( pi_val["λF"][l] * Y["ftR"][l] - pi_val["μF"][l] * Y["ftI"][l] + pi_val["λT"][l] * Y["tfR"][l] - pi_val["μT"][l] * Y["tfI"][l] )
       H[nbuses+to, nbuses+from] += 0.5*( pi_val["λF"][l] * Y["ftR"][l] - pi_val["μF"][l] * Y["ftI"][l] + pi_val["λT"][l] * Y["tfR"][l] - pi_val["μT"][l] * Y["tfI"][l] )
       H[to, nbuses+from] += 0.5*( pi_val["λF"][l] * Y["ftI"][l] - pi_val["λT"][l] * Y["tfI"][l] + pi_val["μF"][l] * Y["ftR"][l] - pi_val["μT"][l] * Y["tfR"][l] )
       H[nbuses+from, to] += 0.5*( pi_val["λF"][l] * Y["ftI"][l] - pi_val["λT"][l] * Y["tfI"][l] + pi_val["μF"][l] * Y["ftR"][l] - pi_val["μT"][l] * Y["tfR"][l] )
       H[from,nbuses+to] -= 0.5*( pi_val["λF"][l] * Y["ftI"][l] - pi_val["λT"][l] * Y["tfI"][l] + pi_val["μF"][l] * Y["ftR"][l] - pi_val["μT"][l] * Y["tfR"][l] )
       H[nbuses+to,from] -= 0.5*( pi_val["λF"][l] * Y["ftI"][l] - pi_val["λT"][l] * Y["tfI"][l] + pi_val["μF"][l] * Y["ftR"][l] - pi_val["μT"][l] * Y["tfR"][l] )
    end
end

# Define callback function for generating and adding cuts
function solveEta0Eigs(H,optdata,v)
    lines, buses, generators, baseMVA = opfdata.lines, opfdata.buses, opfdata.generators, opfdata.baseMVA
    nbuses, nlines, ngens = length(buses), length(lines), length(generators)
    N = 1:nbuses; L = 1:nlines; G = 1:ngens

    E=eigs(H,nev=6,which=:SR, maxiter=100000, tol=1e-8)
    η0Val = E[1][1]
    for i in N
	v["R"][i] = E[2][i,1]; v["I"][i] = E[2][nbuses+i,1] 
    end
    return η0Val
end

function solveEta0SDP(H,opfdata,v)
    lines, buses, generators, baseMVA = opfdata.lines, opfdata.buses, opfdata.generators, opfdata.baseMVA
    nbuses, nlines, ngens = length(buses), length(lines), length(generators)
    N = 1:nbuses; L = 1:nlines; G = 1:ngens

    #The QP subproblem
    mSDP = Model(solver=IpoptSolver())
    @variable(mSDP, e[i=N], start=0); @variable(mSDP, f[i=N], start=0)
    η0Val = 0

    for i in N
      setvalue(e[i], 1); setvalue(f[i], 0)
    end

    # Adjust QP subproblem
    @NLobjective(mSDP, Min, sum( H[i,i]*(e[i]^2+f[i]^2) for i in N) 
	+ 2*sum( H[fromBus[l],toBus[l]]*(e[fromBus[l]]*e[toBus[l]]+f[fromBus[l]]*f[toBus[l]])   for l in L)
	- 2*sum( H[fromBus[l],nbuses+toBus[l]]*(f[fromBus[l]]*e[toBus[l]]-e[fromBus[l]]*f[toBus[l]])   for l in L)
    )
    status = solve(mSDP)
    if status == :Optimal || status == :UserLimit
        η0Val = getobjectivevalue(mSDP)
	for i in N
	    v["R"][i]=getvalue(e[i]); v["I"][i]=getvalue(f[i])
	end
	if(status == :UserLimit)
	    println("solveEta0SDP solve status $status") 
	end
    else
	println("solveEta0SDP: solve status $status") 
	η0Val = 0
    end
    return η0Val
end

function computeSG(Y,opfdata,v,sg)
    lines, buses, generators, baseMVA = opfdata.lines, opfdata.buses, opfdata.generators, opfdata.baseMVA
    nbuses, nlines, ngens = length(buses), length(lines), length(generators)
    N = 1:nbuses; L = 1:nlines; G = 1:ngens; busIdx = mapBusIdToIdx(buses)
    fromBus=zeros(Int,nlines); toBus=zeros(Int,nlines)
    for l in L
      fromBus[l] = busIdx[lines[l].from]; toBus[l] = busIdx[lines[l].to] 
    end
	for i in N
 	    W_val = v["R"][i]^2 + v["I"][i]^2
            sg["α"][i,1] = Y["shR"][i] * W_val; sg["β"][i,1] = -Y["shI"][i] * W_val 
	    sg["δ"][i,1] = W_val; sg["γ"][i,1] = -W_val
	end
	for l in L
	    from = fromBus[l]; to = toBus[l]
 	    e_valF = v["R"][from]; f_valF = v["I"][from]; W_valF = e_valF^2 + f_valF^2
 	    e_valT = v["R"][to]; f_valT = v["I"][to]; W_valT = e_valT^2 + f_valT^2
	    Wr_val = e_valF*e_valT + f_valF*f_valT; Wi_val = e_valT*f_valF - e_valF*f_valT
	    sg["λF"][l,1] = (Y["ffR"][l] * W_valF + Y["ftR"][l] * Wr_val + Y["ftI"][l] * Wi_val)
	    sg["λT"][l,1] = (Y["ttR"][l] * W_valT + Y["tfR"][l] * Wr_val - Y["tfI"][l] * Wi_val)
	    sg["μF"][l,1] = (-Y["ffI"][l] * W_valF - Y["ftI"][l] * Wr_val + Y["ftR"][l] * Wi_val)
	    sg["μT"][l,1] = (-Y["ttI"][l] * W_valT - Y["tfI"][l] * Wr_val - Y["tfR"][l] * Wi_val)
	end
end



