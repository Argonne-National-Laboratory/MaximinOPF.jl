#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

# The master problem MP
  mMP = Model(solver=CplexSolver(CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50))
 # Define the model here
  @variable(mMP, x[l=L], Bin, start=0)
  @variable(mMP, -1 <= α[i=N] <= 1, start=0)
  @variable(mMP, -1 <= β[i=N] <= 1, start=0)
  @variable(mMP, δ[i=N] >= 0)
  @variable(mMP, γ[i=N] >= 0)
  @constraint(mMP, [i=N], δ[i]+γ[i] <= 1)

  @variable(mMP, λF[l=L], start=0)
  @variable(mMP, λT[l=L], start=0)
  @variable(mMP, μF[l=L], start=0)
  @variable(mMP, μT[l=L], start=0)

  @variable(mMP, ζpUB[g=G] >=0)
  @variable(mMP, ζpLB[g=G] >=0)
  @variable(mMP, ζqUB[g=G] >=0)
  @variable(mMP, ζqLB[g=G] >=0)

  for i in N
   for g in BusGeners[i]
    @constraint(mMP, 
	-α[i] + ζpUB[g] - ζpLB[g] == 0 )
    @constraint(mMP, 
	-β[i] + ζqUB[g] - ζqLB[g] == 0 ) 
   end
  end
#=
  if HEUR == 2
    @variable(mMP, xiPLB[l=L] >= 0)
    @variable(mMP, xiPUB[l=L] >= 0)
    @variable(mMP, xiQLB[l=L] >= 0)
    @variable(mMP, xiQUB[l=L] >= 0)
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

  @constraint(mMP, AMcf1[l in L], α[busIdx[lines[l].from]] - x[l] <= λF[l])
  @constraint(mMP, AMcf2[l in L], α[busIdx[lines[l].from]] + x[l] >= λF[l])
  @constraint(mMP, AMcf3[l in L], -(1 - x[l]) <= λF[l])
  @constraint(mMP, AMcf4[l in L], (1 - x[l]) >= λF[l])
  @constraint(mMP, AMct1[l in L], α[busIdx[lines[l].to]] - x[l] <= λT[l])
  @constraint(mMP, AMct2[l in L], α[busIdx[lines[l].to]] + x[l] >= λT[l])
  @constraint(mMP, AMct3[l in L], -(1 - x[l]) <= λT[l])
  @constraint(mMP, AMct4[l in L], (1 - x[l]) >= λT[l])

  @constraint(mMP, BMcf1[l in L], β[busIdx[lines[l].from]] - x[l] <= μF[l])
  @constraint(mMP, BMcf2[l in L], β[busIdx[lines[l].from]] + x[l] >= μF[l])
  @constraint(mMP, BMcf3[l in L], -(1 - x[l]) <= μF[l])
  @constraint(mMP, BMcf4[l in L], (1 - x[l]) >= μF[l])
  @constraint(mMP, BMct1[l in L], β[busIdx[lines[l].to]] - x[l] <= μT[l])
  @constraint(mMP, BMct2[l in L], β[busIdx[lines[l].to]] + x[l] >= μT[l])
  @constraint(mMP, BMct3[l in L], -(1 - x[l]) <= μT[l])
  @constraint(mMP, BMct4[l in L], (1 - x[l]) >= μT[l])

  #These constraints are not quite valid, but their inclusion often results in much faster time to near optimal solution.

  if HEUR == 1 
  	@constraint(mMP, LambdaMuConstr1[l in L], λF[l]*YftI[l] - λT[l]*YtfI[l] + μF[l]*YftR[l] - μT[l]*YtfR[l] == 0.0)
  elseif HEUR == 2
#=
	@constraint(mMP, LambdaFequalsT[l in L], λF[l] - λT[l] + xiPUB[l] - xiPLB[l] == 0) 
	@constraint(mMP, muFequalsT[l in L], μF[l] - μT[l] + xiQUB[l] - xiQLB[l] == 0)
=#
	@constraint(mMP, LambdaFequalsT[l in L], λF[l] - λT[l]  == 0) 
	@constraint(mMP, muFequalsT[l in L], μF[l] - μT[l]  == 0)
  elseif HEUR == 3
	@constraint(mMP, LambdaMuConstr2[l in L], λF[l]*YtfR[l] - λT[l]*YftR[l] - μF[l]*YtfI[l] + μT[l]*YftI[l] == 0.0)
  end



sg_α = zeros(nbuses,2)
sg_β = zeros(nbuses,2)
sg_γ = zeros(nbuses,2) 
sg_δ = zeros(nbuses,2) 
sg_λF = zeros(nlines,2)
sg_μF = zeros(nlines,2)
sg_λT = zeros(nlines,2)
sg_μT = zeros(nlines,2)
W_val = zeros(nbuses)
Wr_val = zeros(nlines)
Wi_val = zeros(nlines)
e_val = zeros(nbuses)
f_val = zeros(nbuses)

maxNSG = 2

#The QP subproblem
  mSDP = Model(solver=IpoptSolver())
  #@variable(mSDP, -1 <= e[i=N] <= 1, start=0) #Add bounds later
  #@variable(mSDP, -1 <= f[i=N] <= 1, start=0)
  @variable(mSDP, e[i=N], start=0) #Add bounds later
  @variable(mSDP, f[i=N], start=0)
  @NLexpression(mSDP, exprW[i=N], e[i]^2 + f[i]^2)
  @NLexpression(mSDP, exprWR[l=L], e[busIdx[lines[l].from]]*e[busIdx[lines[l].to]] + f[busIdx[lines[l].from]]*f[busIdx[lines[l].to]])
  @NLexpression(mSDP, exprWI[l=L], e[busIdx[lines[l].to]]*f[busIdx[lines[l].from]] - e[busIdx[lines[l].from]]*f[busIdx[lines[l].to]])

  @NLconstraint(mSDP, vMagSumUB, sum(exprW[i] for i in N) <= 1) ### Trust-region constraint
  ### Does having a reference angle make sense in the context of generating cuts? 
  ###   Since the vectors e and f that are computed are normalized by construction?
  ###   For now, there is no reference angle.


function generateCuts(cb)
  global sg_α, sg_β, sg_γ, sg_δ, sg_λF, sg_μF, sg_λT, sg_μT
  global η0Val, ncuts
  try
    solveEta0Eigs()
  catch exc
    println("Exception caught with eigs(), solving η0Val subproblem with Ipopt as recourse.") 
    println(exc)
    solveEta0SDP()
  end
  if  η0Val <= -TOL 
	    @lazyconstraint(cb, 0.0 <= sum( (sg_α[i,1])* α[i] for i in N) + sum( (sg_β[i,1])* β[i] for i in N )  
		+ sum( (sg_γ[i,1])* γ[i] for i in N)  + sum( (sg_δ[i,1])* δ[i] for i in N)
		+ sum( (sg_λF[l,1])* λF[l] + (sg_λT[l,1])* λT[l] for l in L) + sum( (sg_μF[l,1])* μF[l] + (sg_μT[l,1])* μT[l] for l in L),
		localcut=useLocalCuts)
	    ncuts += 1
  else
    println("Tolerance met for not generating a new lazy cut.")
  end
end

# Define callback function for generating and adding cuts
function solveEta0Eigs()

	global nSG
	global η0Val
	global α_val, β_val, γ_val, δ_val
	global λF_val, μF_val, λT_val, μT_val
	global W_val, Wr_val, Wi_val
	global sg_α, sg_β, sg_γ, sg_δ, sg_λF, sg_μF, sg_λT, sg_μT
	H=spzeros(2*nbuses,2*nbuses)
	H2=spzeros(Complex,nbuses,nbuses)
        for i in N
          H[i,i] +=  α_val[i] * acYshR[i] - β_val[i] * acYshI[i]  + δ_val[i] - γ_val[i]
          H[nbuses+i,nbuses+i] += α_val[i] * acYshR[i] - β_val[i] * acYshI[i] + δ_val[i] - γ_val[i]
          H2[i,i] +=  α_val[i] * acYshR[i] - β_val[i] * acYshI[i]  + δ_val[i] - γ_val[i]
        end
        for l in L
          from = busIdx[lines[l].from];to = busIdx[lines[l].to] 
          H[from,from] += λF_val[l] * acYffR[l] - μF_val[l] * acYffI[l]
          H[nbuses+from,nbuses+from] += λF_val[l] * acYffR[l] - μF_val[l] * acYffI[l]
          H2[from,from] += λF_val[l] * acYffR[l] - μF_val[l] * acYffI[l]
          H[to,to] += λT_val[l] * acYttR[l] - μT_val[l] * acYttI[l]
          H[nbuses+to,nbuses+to] += λT_val[l] * acYttR[l] - μT_val[l] * acYttI[l]
          H2[to,to] += λT_val[l] * acYttR[l] - μT_val[l] * acYttI[l]
          H[from,to] += 0.5*( λF_val[l] * acYftR[l] - μF_val[l] * acYftI[l] + λT_val[l] * acYtfR[l] - μT_val[l] * acYtfI[l] )
          H[to,from] += 0.5*( λF_val[l] * acYftR[l] - μF_val[l] * acYftI[l] + λT_val[l] * acYtfR[l] - μT_val[l] * acYtfI[l] )
          H[nbuses+from, nbuses+to] += 0.5*( λF_val[l] * acYftR[l] - μF_val[l] * acYftI[l] + λT_val[l] * acYtfR[l] - μT_val[l] * acYtfI[l] )
          H[nbuses+to, nbuses+from] += 0.5*( λF_val[l] * acYftR[l] - μF_val[l] * acYftI[l] + λT_val[l] * acYtfR[l] - μT_val[l] * acYtfI[l] )
          H2[from,to] += 0.5*( λF_val[l] * acYftR[l] - μF_val[l] * acYftI[l] + λT_val[l] * acYtfR[l] - μT_val[l] * acYtfI[l] )
          H2[to,from] += 0.5*( λF_val[l] * acYftR[l] - μF_val[l] * acYftI[l] + λT_val[l] * acYtfR[l] - μT_val[l] * acYtfI[l] )
          H[to, nbuses+from] += 0.5*( λF_val[l] * acYftI[l] - λT_val[l] * acYtfI[l] + μF_val[l] * acYftR[l] - μT_val[l] * acYtfR[l] )
          H[nbuses+from, to] += 0.5*( λF_val[l] * acYftI[l] - λT_val[l] * acYtfI[l] + μF_val[l] * acYftR[l] - μT_val[l] * acYtfR[l] )
          H[from,nbuses+to] -= 0.5*( λF_val[l] * acYftI[l] - λT_val[l] * acYtfI[l] + μF_val[l] * acYftR[l] - μT_val[l] * acYtfR[l] )
          H[nbuses+to,from] -= 0.5*( λF_val[l] * acYftI[l] - λT_val[l] * acYtfI[l] + μF_val[l] * acYftR[l] - μT_val[l] * acYtfR[l] )
          H2[from,to] -= 0.5*( λF_val[l] * acYftI[l] - λT_val[l] * acYtfI[l] + μF_val[l] * acYftR[l] - μT_val[l] * acYtfR[l] )*im
          H2[to,from] += 0.5*( λF_val[l] * acYftI[l] - λT_val[l] * acYtfI[l] + μF_val[l] * acYftR[l] - μT_val[l] * acYtfR[l] )*im
        end
	# Reference bus angle is zero
	e_val = zeros(nbuses)
	f_val = zeros(nbuses)
        H3=Hermitian(H2)
	#E=eigs(H3,nev=1,which=:SR, maxiter=100000, tol=1e-8)
	E=eigs(H,nev=1,which=:SR, maxiter=100000, tol=1e-8)
	η0Val = real(E[1][1])
	nSG = 0
	 if real(E[1][1]) <= -TOL 
	  nSG = 1
	  for i in N
 	    e_val[i] = E[2][i,1]
	    f_val[i] = E[2][nbuses+i,1]
#=
 	    e_val[i] = real(E[2][i,1])
	    f_val[i] = imag(E[2][i,1])
=#
	  end
	  
	  for i in N
	    W_val[i] = e_val[i]^2 + f_val[i]^2
            sg_α[i,1] = acYshR[i] * W_val[i] 
            sg_β[i,1] = -acYshI[i] * W_val[i] 
	    sg_δ[i,1] = W_val[i] 
	    sg_γ[i,1] = -W_val[i]
	  end
	  #sgnorm = sum(sg_α[i,1]^2 + sg_β[i,1]^2 + sg_δ[i,1]^2 + sg_γ[i,1]^2 for i in N)
	  for l in L
	    from = busIdx[lines[l].from]; to = busIdx[lines[l].to]
	    Wr_val[l] = e_val[from]*e_val[to] + f_val[from]*f_val[to]
	    Wi_val[l] = e_val[to]*f_val[from] - e_val[from]*f_val[to]
	    sg_λF[l,1] = (acYffR[l] * W_val[from] + acYftR[l] * Wr_val[l] + acYftI[l] * Wi_val[l])
	    sg_λT[l,1] = (acYttR[l] * W_val[to] + acYtfR[l] * Wr_val[l] - acYtfI[l] * Wi_val[l])
	    sg_μF[l,1] = (-acYffI[l] * W_val[from] - acYftI[l] * Wr_val[l] + acYftR[l]* Wi_val[l])
	    sg_μT[l,1] = (-acYttI[l] * W_val[to] - acYtfI[l] * Wr_val[l] - acYtfR[l] * Wi_val[l])
	  end
          #sgnorm += sum(sg_λF[l,1]^2 + sg_λT[l,1]^2 + sg_μF[l,1]^2 + sg_μT[l,1]^2 for l in L)
	  #sgnorm = sqrt(sgnorm)
#print(" ",sgnorm)
#=
	  if sgnorm > 1e-6
	    for i in N
		sg_α[i,1] /= sgnorm
		sg_β[i,1] /= sgnorm
		sg_δ[i,1] /= sgnorm
		sg_γ[i,1] /= sgnorm
	    end
	    for l in L
		sg_λF[l,1] /= sgnorm
		sg_λT[l,1] /= sgnorm
		sg_μF[l,1] /= sgnorm
		sg_μT[l,1] /= sgnorm
	    end
	  end
=#
	end
#=
   if nSG == 1
     for i in N
	print(" (",E[2][i,1],",",E[2][i+nbuses,1],")")
     end
	print("\n")
   elseif nSG == 2
     for i in N
	print(" (",E[2][i,1],",",E[2][i+nbuses,1],",")
	print(E[2][i,2],",",E[2][i+nbuses,2],")")
     end
	print("\n")
   end
=#
#print("\n")
end


###################OLD VERSION################
#=
	global nSG
	global η0Val
	global α_val, β_val, γ_val, δ_val
	global λF_val, μF_val, λT_val, μT_val
	global W_val, Wr_val, Wi_val
        global sg_α, sg_β, sg_γ, sg_δ, sg_λF, sg_μF, sg_λT, sg_μT
        H=spzeros(2*nbuses,2*nbuses)
        for i in N
          H[i,i] +=  α_val[i] * acYshR[i] - β_val[i] * acYshI[i]  + δ_val[i] - γ_val[i]
          H[nbuses+i,nbuses+i] += α_val[i] * acYshR[i] - β_val[i] * acYshI[i] + δ_val[i] - γ_val[i]
        end
        for l in L
          from = busIdx[lines[l].from];to = busIdx[lines[l].to] 
          H[from,from] += λF_val[l] * acYffR[l] - μF_val[l] * acYffI[l]
          H[nbuses+from,nbuses+from] += λF_val[l] * acYffR[l] - μF_val[l] * acYffI[l]
          H[to,to] += λT_val[l] * acYttR[l] - μT_val[l] * acYttI[l]
          H[nbuses+to,nbuses+to] += λT_val[l] * acYttR[l] - μT_val[l] * acYttI[l]
          H[from,to] += 0.5*( λF_val[l] * acYftR[l] - μF_val[l] * acYftI[l] + λT_val[l] * acYtfR[l] - μT_val[l] * acYtfI[l] )
          H[to,from] += 0.5*( λF_val[l] * acYftR[l] - μF_val[l] * acYftI[l] + λT_val[l] * acYtfR[l] - μT_val[l] * acYtfI[l] )
          H[nbuses+from, nbuses+to] += 0.5*( λF_val[l] * acYftR[l] - μF_val[l] * acYftI[l] + λT_val[l] * acYtfR[l] - μT_val[l] * acYtfI[l] )
          H[nbuses+to, nbuses+from] += 0.5*( λF_val[l] * acYftR[l] - μF_val[l] * acYftI[l] + λT_val[l] * acYtfR[l] - μT_val[l] * acYtfI[l] )
          H[to, nbuses+from] += 0.5*( λF_val[l] * acYftI[l] - λT_val[l] * acYtfI[l] + μF_val[l] * acYftR[l] - μT_val[l] * acYtfR[l] )
          H[nbuses+from, to] += 0.5*( λF_val[l] * acYftI[l] - λT_val[l] * acYtfI[l] + μF_val[l] * acYftR[l] - μT_val[l] * acYtfR[l] )
          H[from,nbuses+to] -= 0.5*( λF_val[l] * acYftI[l] - λT_val[l] * acYtfI[l] + μF_val[l] * acYftR[l] - μT_val[l] * acYtfR[l] )
          H[nbuses+to,from] -= 0.5*( λF_val[l] * acYftI[l] - λT_val[l] * acYtfI[l] + μF_val[l] * acYftR[l] - μT_val[l] * acYtfR[l] )
        end
        #H = Symmetric(H)
        e_val = zeros(nbuses)
        f_val = zeros(nbuses)
        E=eigs(H,nev=maxNEigs,which=:SR, maxiter=3000, tol=1e-6)
        #E=eigs(H,nev=1,which=:SR, maxiter=1000, tol=1e-5, v0=ones(2*nbuses))
        #E=eigs(H,nev=1,which=:SR, maxiter=1000)

        D2_val = E[1][1]
        D2AC_val = D2_val
        nSG = 0
        for s=1:maxNEigs
         if E[1][s] <= -TOL
          nSG += 1
          for i in N
            e_val[i] = E[2][i,s]
            f_val[i] = E[2][nbuses+i,s]
            W_val[i] = e_val[i]^2 + f_val[i]^2
            sg_α[i,s] = acYshR[i] * W_val[i]
            sg_β[i,s] = -acYshI[i] * W_val[i]
            sg_δ[i,s] = W_val[i]
            sg_γ[i,s] = -W_val[i]
          end
          for l in L
            from = busIdx[lines[l].from]; to = busIdx[lines[l].to]
            Wr_val[l] = e_val[from]*e_val[to] + f_val[from]*f_val[to]
            Wi_val[l] = e_val[to]*f_val[from] - e_val[from]*f_val[to]
            sg_λF[l,s] = acYffR[l] * W_val[from] + acYftR[l] * Wr_val[l] + acYftI[l] * Wi_val[l]
            sg_λT[l,s] = acYttR[l] * W_val[to] + acYtfR[l] * Wr_val[l] - acYtfI[l] * Wi_val[l]
            sg_μF[l,s] = -acYffI[l] * W_val[from] - acYftI[l] * Wr_val[l] + acYftR[l]* Wi_val[l]
            sg_μT[l,s] = -acYttI[l] * W_val[to] - acYtfI[l] * Wr_val[l] - acYtfR[l] * Wi_val[l]
	  end
	else
	  break
	end
       end # s=1:maxNSG
end
=#
                
#############################################


function solveEta0SDP()
	global nSG, η0Val
	global W_val, Wr_val, Wi_val, e_val, f_val
	global α_val, β_val, γ_val, δ_val
	global λF_val, μF_val, λT_val, μT_val
	global sg_α, sg_β, sg_γ, sg_δ, sg_λF, sg_μF, sg_λT, sg_μT

	e_ = getindex(mSDP, :e)
	f_ = getindex(mSDP, :f)
	η0Val = 0

	for i in N
		setvalue(e_[i], 1)
		setvalue(f_[i], 0)
		e_val[i]=1
		f_val[i]=0
	end

	# Adjust QP subproblem
	@NLobjective(mSDP, Min,
	  sum( ( α_val[i] * acYshR[i] - β_val[i] * acYshI[i] + δ_val[i] - γ_val[i]) * exprW[i] for i in N)
	  + sum(
		  ( λF_val[l] * acYffR[l] - μF_val[l] * acYffI[l]) 						    * exprW[busIdx[lines[l].from]]
	        + ( 						λT_val[l] * acYttR[l] - μT_val[l] * acYttI[l] ) * exprW[busIdx[lines[l].to]]
	  	+ ( λF_val[l] * acYftR[l] - μF_val[l] * acYftI[l] + λT_val[l] * acYtfR[l] - μT_val[l] * acYtfI[l] ) * exprWR[l]
	  	+ ( λF_val[l] * acYftI[l] - λT_val[l] * acYtfI[l] + μF_val[l] * acYftR[l] - μT_val[l] * acYtfR[l] ) * exprWI[l]
		for l in L))
	status = solve(mSDP)
	if status == :Optimal || status == :UserLimit
	    nSG = 1
	    η0Val = getobjectivevalue(mSDP)
print(" versus opt val: ",η0Val,"\n")
	    for i in N
 	        e_val[i] = getvalue(e_[i])
	        f_val[i] = getvalue(f_[i])
	    end
	    if(status == :UserLimit)
	       println("solveEta0SDP solve status $status") 
	    end
	else
	   println("solveEta0SDP: solve status $status") 
	   η0Val = 0
	   nSG = 0
	   return
	end
	for i in N
	    W_val[i] = e_val[i]^2 + f_val[i]^2
            sg_α[i,1] = acYshR[i] * W_val[i] 
            sg_β[i,1] = -acYshI[i] * W_val[i] 
	    sg_δ[i,1] = W_val[i] 
	    sg_γ[i,1] = -W_val[i]
	end
	for l in L
	    from = busIdx[lines[l].from]; to = busIdx[lines[l].to]
	    Wr_val[l] = e_val[from]*e_val[to] + f_val[from]*f_val[to]
	    Wi_val[l] = e_val[to]*f_val[from] - e_val[from]*f_val[to]
	    sg_λF[l,1] = (acYffR[l] * W_val[from] + acYftR[l] * Wr_val[l] + acYftI[l] * Wi_val[l])
	    sg_λT[l,1] = (acYttR[l] * W_val[to] + acYtfR[l] * Wr_val[l] - acYtfI[l] * Wi_val[l])
	    sg_μF[l,1] = (-acYffI[l] * W_val[from] - acYftI[l] * Wr_val[l] + acYftR[l]* Wi_val[l])
	    sg_μT[l,1] = (-acYttI[l] * W_val[to] - acYtfI[l] * Wr_val[l] - acYtfR[l] * Wi_val[l])
	end
end

function addMPCutsLazyCB(cb)
	global time_Eta0SP, time_MP, mpstart, avg_TimeMP, total_TimeMP, nMPUpdates, numcalls_Eta0SP, ncuts, attacks, noAttacks, bestAttack, noNoGoodCuts
  	global η0Val
  	global W_val, Wr_val, Wi_val
  	global α_val, β_val, γ_val, δ_val 
  	global λF_val, μF_val, λT_val, μT_val, x_val
  	global sg_α, sg_β, sg_γ, sg_δ, sg_λF, sg_μF, sg_λT, sg_μT

	time_MP= (time_ns() - mpstart)/1e9
	total_TimeMP += time_MP
	#nMPUpdates += 1
	#avg_TimeMP = total_TimeMP/nMPUpdates
	# get the variable values
	for i in N
		α_val[i] = getvalue(α[i])
		β_val[i] = getvalue(β[i])
		γ_val[i] = getvalue(γ[i])
		δ_val[i] = getvalue(δ[i]) 
	end
	for l in L
    	   x_val[l] = round(getvalue(x[l]))
		  λF_val[l] = getvalue(λF[l])
		  μF_val[l] = getvalue(μF[l]) 
		  λT_val[l] = getvalue(λT[l])
		  μT_val[l] = getvalue(μT[l]) 
	end

	# generate cuts
	time_Eta0Start = time_ns()
	#  generate cut(s)
	generateCuts(cb) ### the implementation depends on which eta0 SP code is included above
	time_Eta0SP += (time_ns()-time_Eta0Start)/1e9
	numcalls_Eta0SP += 1
	#mpstart = time_ns() #For updating MP
end

addlazycallback(mMP, addMPCutsLazyCB)
mpRootStartTime = time_ns() #For solving root node
mpstart = time_ns() #For updating MP
wtime = @elapsed status = solve(mMP)
@show wtime

@printf("\n********************FINAL RESULT FOR CASE %s WITH %d LINES CUT H%dR%d*******************\n",CASE_NUM,K, HEUR,FORM)

@printf("LAST ATTACK && ")
for l in L
	if round(getvalue(x[l])) == 1
		println("Cut line $l")
	end
end
@printf("\n")
@printf("MPs took: avg: %f, total %f.\n",avg_TimeMP, total_TimeMP)

###Printing summary as one line to aid in latex writeup
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

if FORM == AC
    @printf("AC &    ")
elseif FORM == SDP  
    @printf("SDP &    ")
elseif FORM == SOCP || FORM == SOCPBds 
    @printf("SOCP &    ")
elseif FORM == DC 
    @printf("DC &    ")
end

tShift = wtime - getsolvetime(mMP)

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

