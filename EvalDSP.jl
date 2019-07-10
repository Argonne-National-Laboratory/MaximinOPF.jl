function solveFullModelAC(opfdata,x_val)
# The full lower level model, AC formulation
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  mFullAC = Model(solver=IpoptSolver())
  @variable(mFullAC, opfdata.Pmin[g] <= PG[g=G] <= opfdata.Pmax[g])
  @variable(mFullAC, opfdata.Qmin[g] <= QG[g=G] <= opfdata.Qmax[g])
  @variable(mFullAC, PSlack[i=N] >= 0, start=0)
  @variable(mFullAC, QSlack[i=N] >= 0, start=0)
  @variable(mFullAC, VMSlack[i=N] >= 0, start=0)
  @variable(mFullAC, vR[i=N], start=1)
  @variable(mFullAC, vI[i=N], start=0)

  @NLexpression(mFullAC, W[i=N], vR[i]^2 + vI[i]^2)
  @NLexpression(mFullAC, WR[l=L], vR[fromBus[l]]*vR[toBus[l]] + vI[fromBus[l]]*vI[toBus[l]])
  @NLexpression(mFullAC, WI[l=L], vR[toBus[l]]*vI[fromBus[l]] - vR[fromBus[l]]*vI[toBus[l]])
  @NLconstraint(mFullAC, VoltMagUB[i=N], W[i] - Wmax[i] - VMSlack[i] <= 0)
  @NLconstraint(mFullAC, VoltMagLB[i=N], W[i] - Wmin[i] + VMSlack[i] >= 0)

  @NLexpression(mFullAC, PFlowFrom[l=L], (W[fromBus[l]] * Y["ffR"][l] +  Y["ftR"][l]*WR[l] + Y["ftI"][l]*(WI[l]) ) )
  @NLexpression(mFullAC, PFlowTo[l=L], (W[toBus[l]] * Y["ttR"][l] +  Y["tfR"][l]*WR[l] - Y["tfI"][l]*(WI[l]) ) )  
  @NLexpression(mFullAC, QFlowFrom[l=L], (-Y["ffI"][l]*W[fromBus[l]] - Y["ftI"][l]*WR[l] + Y["ftR"][l]*(WI[l]) ) )
  @NLexpression(mFullAC, QFlowTo[l=L], (-Y["ttI"][l]*W[toBus[l]] - Y["tfI"][l]*WR[l] - Y["tfR"][l]*(WI[l]) ) )

   # ... Power flow balance constraints
   # real part
   @NLconstraint( mFullAC, BusPBalanceUB[i=N],
	Y["shR"][i] * W[i] - sum(PG[g] for g in BusGeners[i]) + opfdata.PD[i]       # Sbus part
      + sum( (1-x_val[l])*PFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*PFlowTo[l]   for l in toLines[i] )  
      - PSlack[i] <= 0 )
   @NLconstraint( mFullAC, BusPBalanceLB[i=N],
	Y["shR"][i] * W[i] - sum(PG[g] for g in BusGeners[i]) + opfdata.PD[i]       # Sbus part
      + sum( (1-x_val[l])*PFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*PFlowTo[l]   for l in toLines[i] )  
      + PSlack[i] >= 0 )

   #imaginary part
     @NLconstraint( mFullAC, BusQBalanceUB[i=N],
	-Y["shI"][i] * W[i] - sum(QG[g] for g in BusGeners[i]) + opfdata.QD[i]        #Sbus part
      + sum( (1-x_val[l])*QFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*QFlowTo[l]   for l in toLines[i] )  
      - QSlack[i] <=0 )
     @NLconstraint( mFullAC, BusQBalanceLB[i=N],
	-Y["shI"][i] * W[i] - sum(QG[g] for g in BusGeners[i]) + opfdata.QD[i]        #Sbus part
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

function solveFullModelSDP(opfdata,x_val,relax=false)
  # The full lower level model, SDP formulation
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  #mFullSDP = Model(solver=SCSSolver(verbose=0))
  mFullSDP = Model(solver=MosekSolver(MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=4))
  @variable(mFullSDP, opfdata.Pmin[g] <= PG[g=G] <= opfdata.Pmax[g])
  @variable(mFullSDP, opfdata.Qmin[g] <= QG[g=G] <= opfdata.Qmax[g])
  @variable(mFullSDP, PSlack[i=N] >= 0)
  @variable(mFullSDP, QSlack[i=N] >= 0)
  @variable(mFullSDP, VMSlack[i=N] >= 0)
  @variable(mFullSDP, W[1:(2*nbuses), 1:(2*nbuses)], SDP)
  @variable(mFullSDP, hP[l=L]) 
  @variable(mFullSDP, hQ[l=L]) 
  if !relax
    for l in L
      setlowerbound(hP[l],0);setupperbound(hP[l],0);
      setlowerbound(hQ[l],0);setupperbound(hQ[l],0);
    end
  end

  @constraint(mFullSDP, VoltMagUB[i=N], W[i,i] - opfdata.Wmax[i] - VMSlack[i] <= 0)
  @constraint(mFullSDP, VoltMagLB[i=N], -W[i,i] + opfdata.Wmin[i] - VMSlack[i] <= 0)


  @constraint(mFullSDP, Sym1[i=N,j=i:nbuses], W[i,j] - W[nbuses+i,nbuses+j] == 0) # enforcing form [WR WI; WI WR]
  @constraint(mFullSDP, Sym2[i=N,j=i:nbuses], W[i,nbuses+j] + W[j,nbuses+i] == 0)  # Enforcing WI skew symmetric

  @expression(mFullSDP, WR[l=L], W[fromBus[l],toBus[l]] )
  @expression(mFullSDP, WI[l=L], W[fromBus[l],nbuses+toBus[l]] ) 
  @expression(mFullSDP, PFlowFrom[l=L], 
	Y["ffR"][l]*W[fromBus[l],fromBus[l]]  +  Y["ftR"][l]*WR[l] + Y["ftI"][l]*WI[l] ) 
  @expression(mFullSDP, PFlowTo[l=L], 
	Y["ttR"][l]*W[toBus[l],toBus[l]] +  Y["tfR"][l]*WR[l] - Y["tfI"][l]*WI[l]  )  
  @expression(mFullSDP, QFlowFrom[l=L], -Y["ffI"][l]*W[fromBus[l],fromBus[l]] - Y["ftI"][l]*WR[l] + Y["ftR"][l]*WI[l] )
  @expression(mFullSDP, QFlowTo[l=L], -Y["ttI"][l]*W[toBus[l],toBus[l]] - Y["tfI"][l]*WR[l] - Y["tfR"][l]*WI[l] )

   # ... Power flow balance constraints
   # real part
   @constraint( mFullSDP, BusPBalanceUB[i=N],
	Y["shR"][i] * W[i,i] 
- sum(PG[g] for g in BusGeners[i]) 
+ opfdata.PD[i]        # Sbus part
      + sum( (1-x_val[l])*(PFlowFrom[l]+hP[l])   for l in fromLines[i] )  + sum( (1-x_val[l])*(PFlowTo[l]-hP[l])   for l in toLines[i] )  
     <= PSlack[i])
   @constraint( mFullSDP, BusPBalanceLB[i=N],
	Y["shR"][i] * W[i,i] 
- sum(PG[g] for g in BusGeners[i]) 
+ opfdata.PD[i]        # Sbus part
      + sum( (1-x_val[l])*(PFlowFrom[l]+hP[l])   for l in fromLines[i] )  + sum( (1-x_val[l])*(PFlowTo[l]-hP[l])   for l in toLines[i] )  
     >= -PSlack[i])

   #imaginary part
     @constraint( mFullSDP, BusQBalanceUB[i=N],
	-Y["shI"][i] * W[i,i] 
- sum(QG[g] for g in BusGeners[i]) 
+ opfdata.QD[i]        #Sbus part
      + sum( (1-x_val[l])*(QFlowFrom[l]+hQ[l])   for l in fromLines[i] )  + sum( (1-x_val[l])*(QFlowTo[l]-hQ[l])   for l in toLines[i] )  
      <= QSlack[i] )
     @constraint( mFullSDP, BusQBalanceLB[i=N],
	-Y["shI"][i] * W[i,i] 
- sum(QG[g] for g in BusGeners[i]) 
+ opfdata.QD[i]        #Sbus part
      + sum( (1-x_val[l])*(QFlowFrom[l]+hQ[l])   for l in fromLines[i] )  + sum( (1-x_val[l])*(QFlowTo[l]-hQ[l])   for l in toLines[i] )  
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

function solveFullModelSOCP(opfdata,x_val)
# The full lower level model, SOCP relaxation
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  mFullSOCP = Model(solver=IpoptSolver())
  @variable(mFullSOCP, opfdata.Pmin[g] <= PG[g=G] <= opfdata.Pmax[g])
  @variable(mFullSOCP, opfdata.Qmin[g] <= QG[g=G] <= opfdata.Qmax[g])
  @variable(mFullSOCP, PSlack[i=N] >= 0, start=0)
  @variable(mFullSOCP, QSlack[i=N] >= 0, start=0)
  @variable(mFullSOCP, VMSlack[i=N] >= 0)
  @variable(mFullSOCP, W[i=N] >= 0, start=1)
  @constraint(mFullSOCP, VMUB[i=N], W[i] - opfdata.Wmax[i] - VMSlack[i] <= 0)
  @constraint(mFullSOCP, VMLB[i=N], -W[i] + opfdata.Wmin[i] - VMSlack[i] <= 0)
  @variable(mFullSOCP, WR[l=L] >= 0, start=1)
  @variable(mFullSOCP, WI[l=L], start=0)



  @NLconstraint(mFullSOCP, [l=L],   WR[l]^2 + WI[l]^2 - W[fromBus[l]]*W[toBus[l]] <= 0)
  @NLexpression(mFullSOCP, PFlowFrom[l=L], (W[fromBus[l]] * Y["ffR"][l] +  Y["ftR"][l]*WR[l] + Y["ftI"][l]*(WI[l]) ) )
  @NLexpression(mFullSOCP, PFlowTo[l=L], (W[toBus[l]] * Y["ttR"][l] +  Y["tfR"][l]*WR[l] - Y["tfI"][l]*(WI[l]) ) )  
  @NLexpression(mFullSOCP, QFlowFrom[l=L], (-Y["ffI"][l]*W[fromBus[l]] - Y["ftI"][l]*WR[l] + Y["ftR"][l]*(WI[l]) ) )
  @NLexpression(mFullSOCP, QFlowTo[l=L], (-Y["ttI"][l]*W[toBus[l]] - Y["tfI"][l]*WR[l] - Y["tfR"][l]*(WI[l]) ) )

   # ... Power flow balance constraints
   # real part
   @NLconstraint( mFullSOCP, BusPBalanceUB[i=N],
	Y["shR"][i] * W[i] - sum(PG[g] for g in BusGeners[i]) + opfdata.PD[i]        # Sbus part
      + sum( (1-x_val[l])*PFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*PFlowTo[l]   for l in toLines[i] )  
     - PSlack[i] <= 0)
   @NLconstraint( mFullSOCP, BusPBalanceLB[i=N],
	Y["shR"][i] * W[i] - sum(PG[g] for g in BusGeners[i]) + opfdata.PD[i]        # Sbus part
      + sum( (1-x_val[l])*PFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*PFlowTo[l]   for l in toLines[i] )  
     + PSlack[i] >= 0)

   #imaginary part
   @NLconstraint( mFullSOCP, BusQBalanceUB[i=N],
	-Y["shI"][i] * W[i] - sum(QG[g] for g in BusGeners[i]) + opfdata.QD[i]        #Sbus part
      + sum( (1-x_val[l])*QFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*QFlowTo[l]   for l in toLines[i] )  
      - QSlack[i] <= 0)
   @NLconstraint( mFullSOCP, BusQBalanceLB[i=N],
	-Y["shI"][i] * W[i] - sum(QG[g] for g in BusGeners[i]) + opfdata.QD[i]        #Sbus part
      + sum( (1-x_val[l])*QFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*QFlowTo[l]   for l in toLines[i] )  
      + QSlack[i] >= 0)
  @objective(mFullSOCP, Min, sum(PSlack[i] + QSlack[i] + VMSlack[i] for i in N))

  for i in N
    setvalue(W[i],1)	
    setvalue(PSlack[i],0)
    setvalue(QSlack[i],0)
    setvalue(VMSlack[i],0)
  end
  for l in L
    setvalue(WR[l],1)	
    setvalue(WI[l],0)	
  end
  status = solve(mFullSOCP)
  if status == :Optimal
    return getobjectivevalue(mFullSOCP)
  else
    println("Full SOCP model solved with status: $status")
    return 0
  end
end

function solveFullModelDC(opfdata,x_val)
# The full lower level models, DC approximation
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_DC
  mFullDC = Model(solver=IpoptSolver())
  @variable(mFullDC, Pmin[g] <= PG[g=G] <= Pmax[g])
  @variable(mFullDC, PSlack[i=N] >= 0, start=0)
  @variable(mFullDC, THETA[i=N], start=0)

  @constraint(mFullDC, AngleDiffUB[l=L], THETA[fromBus[l]] - THETA[toBus[l]] <= pi)
  @constraint(mFullDC, AngleDiffLB[l=L], -pi <= THETA[fromBus[l]] - THETA[toBus[l]])

   # ... Power flow balance constraints
   # real part
   @NLconstraint( mFullDC, BusPBalanceUB[i=N],
	Y["shR"][i] - sum(PG[g] for g in BusGeners[i]) + opfdata.PD[i]       # Sbus part
      + sum( (1-x_val[l])*(Y["ffR"][l] +  Y["ftR"][l] + Y["ftI"][l]*( THETA[fromBus[l]] - THETA[toBus[l]] )  )   for l in fromLines[i] )  
      + sum( (1-x_val[l])*(Y["ttR"][l] +  Y["tfR"][l] - Y["tfI"][l]*( THETA[fromBus[l]] - THETA[toBus[l]] )  )   for l in toLines[i] )  
     - PSlack[i] <= 0)
   @NLconstraint( mFullDC, BusPBalanceLB[i=N],
	Y["shR"][i] - sum(PG[g] for g in BusGeners[i]) + opfdata.PD[i]       # Sbus part
      + sum( (1-x_val[l])*(Y["ffR"][l] +  Y["ftR"][l] + Y["ftI"][l]*( THETA[fromBus[l]] - THETA[toBus[l]] )  )   for l in fromLines[i] )  
      + sum( (1-x_val[l])*(Y["ttR"][l] +  Y["tfR"][l] - Y["tfI"][l]*( THETA[fromBus[l]] - THETA[toBus[l]] )  )   for l in toLines[i] )  
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
