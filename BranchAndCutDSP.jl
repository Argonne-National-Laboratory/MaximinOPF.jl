





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



function solveFullModelAC(x_val)
# The full lower level model, AC formulation
  mFullAC = Model(solver=IpoptSolver())
  @variable(mFullAC, Pmin[g] <= PG[g=G] <= Pmax[g])
  @variable(mFullAC, Qmin[g] <= QG[g=G] <= Qmax[g])
  @variable(mFullAC, PSlack[i=N] >= 0, start=0)
  @variable(mFullAC, QSlack[i=N] >= 0, start=0)
  @variable(mFullAC, VMSlack[i=N] >= 0, start=0)
  @variable(mFullAC, vR[i=N], start=1)
  @variable(mFullAC, vI[i=N], start=0)

  @NLexpression(mFullAC, W[i=N], vR[i]^2 + vI[i]^2)
  @NLexpression(mFullAC, WR[l=L], vR[busIdx[lines[l].from]]*vR[busIdx[lines[l].to]] + vI[busIdx[lines[l].from]]*vI[busIdx[lines[l].to]])
  @NLexpression(mFullAC, WI[l=L], vR[busIdx[lines[l].to]]*vI[busIdx[lines[l].from]] - vR[busIdx[lines[l].from]]*vI[busIdx[lines[l].to]])
  @NLconstraint(mFullAC, VoltMagUB[i=N], W[i] - Wmax[i] - VMSlack[i] <= 0)
  @NLconstraint(mFullAC, VoltMagLB[i=N], W[i] - Wmin[i] + VMSlack[i] >= 0)

  @NLexpression(mFullAC, PFlowFrom[l=L], (W[busIdx[lines[l].from]] * Y["ffR"][l] +  Y["ftR"][l]*WR[l] + Y["ftI"][l]*(WI[l]) ) )
  @NLexpression(mFullAC, PFlowTo[l=L], (W[busIdx[lines[l].to]] * Y["ttR"][l] +  Y["tfR"][l]*WR[l] - Y["tfI"][l]*(WI[l]) ) )  
  @NLexpression(mFullAC, QFlowFrom[l=L], (-Y["ffI"][l]*W[busIdx[lines[l].from]] - Y["ftI"][l]*WR[l] + Y["ftR"][l]*(WI[l]) ) )
  @NLexpression(mFullAC, QFlowTo[l=L], (-Y["ttI"][l]*W[busIdx[lines[l].to]] - Y["tfI"][l]*WR[l] - Y["tfR"][l]*(WI[l]) ) )

   # ... Power flow balance constraints
   # real part
   @NLconstraint( mFullAC, BusPBalanceUB[i=N],
	Y["shR"][i] * W[i] - sum(PG[g] for g in BusGeners[i]) + PD[i]       # Sbus part
      + sum( (1-x_val[l])*PFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*PFlowTo[l]   for l in toLines[i] )  
      - PSlack[i] <= 0 )
   @NLconstraint( mFullAC, BusPBalanceLB[i=N],
	Y["shR"][i] * W[i] - sum(PG[g] for g in BusGeners[i]) + PD[i]       # Sbus part
      + sum( (1-x_val[l])*PFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*PFlowTo[l]   for l in toLines[i] )  
      + PSlack[i] >= 0 )

   #imaginary part
     @NLconstraint( mFullAC, BusQBalanceUB[i=N],
	-Y["shI"][i] * W[i] - sum(QG[g] for g in BusGeners[i]) + QD[i]        #Sbus part
      + sum( (1-x_val[l])*QFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*QFlowTo[l]   for l in toLines[i] )  
      - QSlack[i] <=0 )
     @NLconstraint( mFullAC, BusQBalanceLB[i=N],
	-Y["shI"][i] * W[i] - sum(QG[g] for g in BusGeners[i]) + QD[i]        #Sbus part
      + sum( (1-x_val[l])*QFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*QFlowTo[l]   for l in toLines[i] )  
      + QSlack[i] >=0 )
  @objective(mFullAC, Min, sum(PSlack[i] + QSlack[i] + VMSlack[i] for i in N))
  vR = getindex(mFullAC, :vR)
  vI = getindex(mFullAC, :vI)

  for i in N
    setvalue(vR[i],1)
    setvalue(vI[i],0)
    setvalue(PSlack[i],0)
    setvalue(QSlack[i],0)
    setvalue(VMSlack[i],0)
  end
  status = solve(mFullAC)
  if status == :Optimal || status == :UserLimit
      println("Full AC model solved with status: $status")
      return getobjectivevalue(mFullAC)
  else
      println("Full AC model solved with status: $status")
      return 0
  end
end

function solveFullModelSDP(x_val)
  # The full lower level model, SDP formulation
  #mFullSDP = Model(solver=SCSSolver(verbose=0))
  mFullSDP = Model(solver=MosekSolver(MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=4))
  @variable(mFullSDP, Pmin[g] <= PG[g=G] <= Pmax[g])
  @variable(mFullSDP, Qmin[g] <= QG[g=G] <= Qmax[g])
  @variable(mFullSDP, PSlack[i=N] >= 0)
  @variable(mFullSDP, QSlack[i=N] >= 0)
  @variable(mFullSDP, VMSlack[i=N] >= 0)
  @variable(mFullSDP, W[1:(2*nbuses), 1:(2*nbuses)], SDP)

  @constraint(mFullSDP, VoltMagUB[i=N], W[i,i] - Wmax[i] - VMSlack[i] <= 0)
  @constraint(mFullSDP, VoltMagLB[i=N], -W[i,i] + Wmin[i] - VMSlack[i] <= 0)


  @constraint(mFullSDP, Sym1[i=N,j=i:nbuses], W[i,j] - W[nbuses+i,nbuses+j] == 0) # enforcing form [WR WI; WI WR]
  @constraint(mFullSDP, Sym2[i=N,j=i:nbuses], W[i,nbuses+j] + W[j,nbuses+i] == 0)  # Enforcing WI skew symmetric

  @expression(mFullSDP, WR[l=L], W[busIdx[lines[l].from],busIdx[lines[l].to]] )
  @expression(mFullSDP, WI[l=L], W[busIdx[lines[l].from],nbuses+busIdx[lines[l].to]] ) 
  @expression(mFullSDP, PFlowFrom[l=L], 
	Y["ffR"][l]*W[busIdx[lines[l].from],busIdx[lines[l].from]]  +  Y["ftR"][l]*WR[l] + Y["ftI"][l]*WI[l] ) 
  @expression(mFullSDP, PFlowTo[l=L], 
	Y["ttR"][l]*W[busIdx[lines[l].to],busIdx[lines[l].to]] +  Y["tfR"][l]*WR[l] - Y["tfI"][l]*WI[l]  )  
  @expression(mFullSDP, QFlowFrom[l=L], -Y["ffI"][l]*W[busIdx[lines[l].from],busIdx[lines[l].from]] - Y["ftI"][l]*WR[l] + Y["ftR"][l]*WI[l] )
  @expression(mFullSDP, QFlowTo[l=L], -Y["ttI"][l]*W[busIdx[lines[l].to],busIdx[lines[l].to]] - Y["tfI"][l]*WR[l] - Y["tfR"][l]*WI[l] )

   # ... Power flow balance constraints
   # real part
   @constraint( mFullSDP, BusPBalanceUB[i=N],
	Y["shR"][i] * W[i,i] 
- sum(PG[g] for g in BusGeners[i]) 
+ PD[i]        # Sbus part
      + sum( (1-x_val[l])*PFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*PFlowTo[l]   for l in toLines[i] )  
     <= PSlack[i])
   @constraint( mFullSDP, BusPBalanceLB[i=N],
	Y["shR"][i] * W[i,i] 
- sum(PG[g] for g in BusGeners[i]) 
+ PD[i]        # Sbus part
      + sum( (1-x_val[l])*PFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*PFlowTo[l]   for l in toLines[i] )  
     >= -PSlack[i])

   #imaginary part
     @constraint( mFullSDP, BusQBalanceUB[i=N],
	-Y["shI"][i] * W[i,i] 
- sum(QG[g] for g in BusGeners[i]) 
+ QD[i]        #Sbus part
      + sum( (1-x_val[l])*QFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*QFlowTo[l]   for l in toLines[i] )  
      <= QSlack[i] )
     @constraint( mFullSDP, BusQBalanceLB[i=N],
	-Y["shI"][i] * W[i,i] 
- sum(QG[g] for g in BusGeners[i]) 
+ QD[i]        #Sbus part
      + sum( (1-x_val[l])*QFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*QFlowTo[l]   for l in toLines[i] )  
      >= -QSlack[i])
  @objective(mFullSDP, Min, sum(PSlack[i] + QSlack[i] + VMSlack[i] for i in N))

  status = solve(mFullSDP)
  if status == :Optimal || status == :Stall
#=
for i in N
@show getdual(BusPBalanceUB[i]) + getdual(BusPBalanceLB[i])
@show getdual(BusQBalanceUB[i]) + getdual(BusQBalanceLB[i])
@show getdual(VoltMagLB[i])
@show -getdual(VoltMagUB[i])
end
=#
    return getobjectivevalue(mFullSDP)
  else
    println("Full SDP model solved with status: $status")
    return 0
  end
end

function solveFullModelSOCP(x_val)
# The full lower level model, SOCP relaxation
  mFullSOCP = Model(solver=IpoptSolver())
  @variable(mFullSOCP, Pmin[g] <= socp_PG[g=G] <= Pmax[g])
  @variable(mFullSOCP, Qmin[g] <= socp_QG[g=G] <= Qmax[g])
  @variable(mFullSOCP, socp_PgSlack[i=N] >= 0, start=0)
  @variable(mFullSOCP, socp_PdSlack[i=N] >= 0, start=0)
  @variable(mFullSOCP, socp_QgSlack[i=N] >= 0, start=0)
  @variable(mFullSOCP, socp_QdSlack[i=N] >= 0, start=0)
  @variable(mFullSOCP, Wmin[i] <= WSOCP[i=N] <= Wmax[i], start=1)
  @variable(mFullSOCP, WRSOCP[l=L] >= 0, start=1)
  @variable(mFullSOCP, WISOCP[l=L], start=0)


  @NLexpression(mFullSOCP, socp_W[i=N], WSOCP[i])
  @NLexpression(mFullSOCP, socp_WR[l=L], WRSOCP[l])
  @NLexpression(mFullSOCP, socp_WI[l=L], WISOCP[l])
  @NLconstraint(mFullSOCP, [l=L],   socp_WR[l]^2 + socp_WI[l]^2 - socp_W[busIdx[lines[l].from]]*socp_W[busIdx[lines[l].to]] <= 0)
  @NLexpression(mFullSOCP, socp_PFlowFrom[l=L], (socp_W[busIdx[lines[l].from]] * Y["ffR"][l] +  Y["ftR"][l]*socp_WR[l] + Y["ftI"][l]*(socp_WI[l]) ) )
  @NLexpression(mFullSOCP, socp_PFlowTo[l=L], (socp_W[busIdx[lines[l].to]] * Y["ttR"][l] +  Y["tfR"][l]*socp_WR[l] - Y["tfI"][l]*(socp_WI[l]) ) )  
  @NLexpression(mFullSOCP, socp_QFlowFrom[l=L], (-Y["ffI"][l]*socp_W[busIdx[lines[l].from]] - Y["ftI"][l]*socp_WR[l] + Y["ftR"][l]*(socp_WI[l]) ) )
  @NLexpression(mFullSOCP, socp_QFlowTo[l=L], (-Y["ttI"][l]*socp_W[busIdx[lines[l].to]] - Y["tfI"][l]*socp_WR[l] - Y["tfR"][l]*(socp_WI[l]) ) )

   # ... Power flow balance constraints
   # real part
   @NLconstraint( mFullSOCP, socp_BusPBalance[i=N],
	Y["shR"][i] * socp_W[i] - sum(socp_PG[g] for g in BusGeners[i]) + socp_PgSlack[i] + (PD[i] - socp_PdSlack[i] )       # Sbus part
      + sum( (1-x_val[l])*socp_PFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*socp_PFlowTo[l]   for l in toLines[i] )  
     ==0)

   #imaginary part
   @NLconstraint( mFullSOCP, socp_BusQBalance[i=N],
	-Y["shI"][i] * socp_W[i] - sum(socp_QG[g] for g in BusGeners[i]) + socp_QgSlack[i] + (QD[i] - socp_QdSlack[i])       #Sbus part
      + sum( (1-x_val[l])*socp_QFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*socp_QFlowTo[l]   for l in toLines[i] )  
     ==0)
  @objective(mFullSOCP, Min, sum(socp_PgSlack[i] + socp_PdSlack[i] + socp_QgSlack[i] + socp_QdSlack[i] for i in N))

  for i in N
    setvalue(WSOCP[i],1)	
    setvalue(socp_PgSlack[i],0)
    setvalue(socp_PdSlack[i],0)
    setvalue(socp_QgSlack[i],0)
    setvalue(socp_QdSlack[i],0)
  end
  for l in L
    setvalue(WRSOCP[l],1)	
    setvalue(WISOCP[l],0)	
  end
  status = solve(mFullSOCP)
  if status == :Optimal
    return getobjectivevalue(mFullSOCP)
  else
    println("Full SOCP model solved with status: $status")
    return 0
  end
end

function solveFullModelDC(x_val)
# The full lower level models, DC approximation
  mFullDC = Model(solver=IpoptSolver())
  @variable(mFullDC, Pmin[g] <= PG[g=G] <= Pmax[g])
  @variable(mFullDC, PSlack[i=N] >= 0, start=0)
  @variable(mFullDC, THETA[i=N], start=0)

  @constraint(mFullDC, AngleDiffUB[l=L], THETA[busIdx[lines[l].from]] - THETA[busIdx[lines[l].to]] <= pi)
  @constraint(mFullDC, AngleDiffLB[l=L], -pi <= THETA[busIdx[lines[l].from]] - THETA[busIdx[lines[l].to]])

   # ... Power flow balance constraints
   # real part
   @NLconstraint( mFullDC, BusPBalanceUB[i=N],
	dcYshR[i] - sum(PG[g] for g in BusGeners[i]) + PD[i]       # Sbus part
      + sum( (1-x_val[l])*(dcYffR[l] +  dcYftR[l] + dcYftI[l]*( THETA[busIdx[lines[l].from]] - THETA[busIdx[lines[l].to]] )  )   for l in fromLines[i] )  
      + sum( (1-x_val[l])*(dcYttR[l] +  dcYtfR[l] - dcYtfI[l]*( THETA[busIdx[lines[l].from]] - THETA[busIdx[lines[l].to]] )  )   for l in toLines[i] )  
     - PSlack[i] <= 0)
   @NLconstraint( mFullDC, BusPBalanceLB[i=N],
	dcYshR[i] - sum(PG[g] for g in BusGeners[i]) + PD[i]       # Sbus part
      + sum( (1-x_val[l])*(dcYffR[l] +  dcYftR[l] + dcYftI[l]*( THETA[busIdx[lines[l].from]] - THETA[busIdx[lines[l].to]] )  )   for l in fromLines[i] )  
      + sum( (1-x_val[l])*(dcYttR[l] +  dcYtfR[l] - dcYtfI[l]*( THETA[busIdx[lines[l].from]] - THETA[busIdx[lines[l].to]] )  )   for l in toLines[i] )  
     + PSlack[i] >= 0)

  @objective(mFullDC, Min, sum(PSlack[i] for i in N))


  for i in N
    setvalue(THETA[i],0)
    setvalue(PSlack[i],0)
  end
  status = solve(mFullDC)
  if status == :Optimal
    return getobjectivevalue(mFullDC)
  else
    println("Full DC model solved with status: $status")
    return 0
  end
end
