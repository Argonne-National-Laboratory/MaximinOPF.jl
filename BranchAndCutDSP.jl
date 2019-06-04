








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

  @NLexpression(mFullAC, PFlowFrom[l=L], (W[busIdx[lines[l].from]] * acYffR[l] +  acYftR[l]*WR[l] + acYftI[l]*(WI[l]) ) )
  @NLexpression(mFullAC, PFlowTo[l=L], (W[busIdx[lines[l].to]] * acYttR[l] +  acYtfR[l]*WR[l] - acYtfI[l]*(WI[l]) ) )  
  @NLexpression(mFullAC, QFlowFrom[l=L], (-acYffI[l]*W[busIdx[lines[l].from]] - acYftI[l]*WR[l] + acYftR[l]*(WI[l]) ) )
  @NLexpression(mFullAC, QFlowTo[l=L], (-acYttI[l]*W[busIdx[lines[l].to]] - acYtfI[l]*WR[l] - acYtfR[l]*(WI[l]) ) )

   # ... Power flow balance constraints
   # real part
   @NLconstraint( mFullAC, BusPBalanceUB[i=N],
	acYshR[i] * W[i] - sum(PG[g] for g in BusGeners[i]) + PD[i]       # Sbus part
      + sum( (1-x_val[l])*PFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*PFlowTo[l]   for l in toLines[i] )  
      - PSlack[i] <= 0 )
   @NLconstraint( mFullAC, BusPBalanceLB[i=N],
	acYshR[i] * W[i] - sum(PG[g] for g in BusGeners[i]) + PD[i]       # Sbus part
      + sum( (1-x_val[l])*PFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*PFlowTo[l]   for l in toLines[i] )  
      + PSlack[i] >= 0 )

   #imaginary part
     @NLconstraint( mFullAC, BusQBalanceUB[i=N],
	-acYshI[i] * W[i] - sum(QG[g] for g in BusGeners[i]) + QD[i]        #Sbus part
      + sum( (1-x_val[l])*QFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*QFlowTo[l]   for l in toLines[i] )  
      - QSlack[i] <=0 )
     @NLconstraint( mFullAC, BusQBalanceLB[i=N],
	-acYshI[i] * W[i] - sum(QG[g] for g in BusGeners[i]) + QD[i]        #Sbus part
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
  if status == :Optimal
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
	acYffR[l]*W[busIdx[lines[l].from],busIdx[lines[l].from]]  +  acYftR[l]*WR[l] + acYftI[l]*WI[l] ) 
  @expression(mFullSDP, PFlowTo[l=L], 
	acYttR[l]*W[busIdx[lines[l].to],busIdx[lines[l].to]] +  acYtfR[l]*WR[l] - acYtfI[l]*WI[l]  )  
  @expression(mFullSDP, QFlowFrom[l=L], -acYffI[l]*W[busIdx[lines[l].from],busIdx[lines[l].from]] - acYftI[l]*WR[l] + acYftR[l]*WI[l] )
  @expression(mFullSDP, QFlowTo[l=L], -acYttI[l]*W[busIdx[lines[l].to],busIdx[lines[l].to]] - acYtfI[l]*WR[l] - acYtfR[l]*WI[l] )

   # ... Power flow balance constraints
   # real part
   @constraint( mFullSDP, BusPBalanceUB[i=N],
	acYshR[i] * W[i,i] 
- sum(PG[g] for g in BusGeners[i]) 
+ PD[i]        # Sbus part
      + sum( (1-x_val[l])*PFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*PFlowTo[l]   for l in toLines[i] )  
     <= PSlack[i])
   @constraint( mFullSDP, BusPBalanceLB[i=N],
	acYshR[i] * W[i,i] 
- sum(PG[g] for g in BusGeners[i]) 
+ PD[i]        # Sbus part
      + sum( (1-x_val[l])*PFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*PFlowTo[l]   for l in toLines[i] )  
     >= -PSlack[i])

   #imaginary part
     @constraint( mFullSDP, BusQBalanceUB[i=N],
	-acYshI[i] * W[i,i] 
- sum(QG[g] for g in BusGeners[i]) 
+ QD[i]        #Sbus part
      + sum( (1-x_val[l])*QFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*QFlowTo[l]   for l in toLines[i] )  
      <= QSlack[i] )
     @constraint( mFullSDP, BusQBalanceLB[i=N],
	-acYshI[i] * W[i,i] 
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
  @NLexpression(mFullSOCP, socp_PFlowFrom[l=L], (socp_W[busIdx[lines[l].from]] * acYffR[l] +  acYftR[l]*socp_WR[l] + acYftI[l]*(socp_WI[l]) ) )
  @NLexpression(mFullSOCP, socp_PFlowTo[l=L], (socp_W[busIdx[lines[l].to]] * acYttR[l] +  acYtfR[l]*socp_WR[l] - acYtfI[l]*(socp_WI[l]) ) )  
  @NLexpression(mFullSOCP, socp_QFlowFrom[l=L], (-acYffI[l]*socp_W[busIdx[lines[l].from]] - acYftI[l]*socp_WR[l] + acYftR[l]*(socp_WI[l]) ) )
  @NLexpression(mFullSOCP, socp_QFlowTo[l=L], (-acYttI[l]*socp_W[busIdx[lines[l].to]] - acYtfI[l]*socp_WR[l] - acYtfR[l]*(socp_WI[l]) ) )

   # ... Power flow balance constraints
   # real part
   @NLconstraint( mFullSOCP, socp_BusPBalance[i=N],
	acYshR[i] * socp_W[i] - sum(socp_PG[g] for g in BusGeners[i]) + socp_PgSlack[i] + (PD[i] - socp_PdSlack[i] )       # Sbus part
      + sum( (1-x_val[l])*socp_PFlowFrom[l]   for l in fromLines[i] )  + sum( (1-x_val[l])*socp_PFlowTo[l]   for l in toLines[i] )  
     ==0)

   #imaginary part
   @NLconstraint( mFullSOCP, socp_BusQBalance[i=N],
	-acYshI[i] * socp_W[i] - sum(socp_QG[g] for g in BusGeners[i]) + socp_QgSlack[i] + (QD[i] - socp_QdSlack[i])       #Sbus part
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
