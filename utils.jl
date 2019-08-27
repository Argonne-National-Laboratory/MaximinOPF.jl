using LightGraphs

#USEFUL SUBROUTINES
function get_graph(mat::SparseMatrixCSC{Float64,Int64})
    rows = rowvals(mat)
    vals = nonzeros(mat)
    m, n = size(mat)
    g = SimpleGraph(m)
    for i = 1:n
        for j in nzrange(mat, i)
            # if abs(vals[j]) > 0
                add_edge!(g, rows[j], i)
            # end
        end
    end
    return g
end

# Get the n-by-n chordal extension
function get_chordal_extension(opfdata)
    N, L, fromBus, toBus = opfdata.N, opfdata.L, opfdata.fromBus, opfdata.toBus  
    # Laplacian graph
    I = Int64[]; J = Int64[]; V = Float64[]
    for l in L
        fid = fromBus[l]
        tid = toBus[l]
        push!(I,fid); push!(J,tid); push!(V,-1)
        push!(I,tid); push!(J,fid); push!(V,-1)
    end
    A = sparse(I,J,V)
    for i in N
        push!(I,i); push!(J,i)
        push!(V,-sum(A[i,:])+1)
    end
    A = sparse(I,J,V)
    C = sparse(cholfact(A))
    return C
end

# Get the 2n-by-2n chordal extension
function get_chordal_extension_complex(opfdata)
    N, L, fromBus, toBus = opfdata.N, opfdata.L, opfdata.fromBus, opfdata.toBus  
    # Laplacian graph
    num = length(N)
    I = Int64[]; J = Int64[]; V = Float64[]
    for l in L
        fid = fromBus[l]
        tid = toBus[l]
        push!(I,fid); push!(J,tid); push!(V,-1)
        push!(I,tid); push!(J,fid); push!(V,-1)
        push!(I,fid); push!(J,num+tid); push!(V,-1)
        push!(I,num+tid); push!(J,fid); push!(V,-1)
        push!(I,tid); push!(J,num+fid); push!(V,-1)
        push!(I,num+fid); push!(J,tid); push!(V,-1)
        push!(I,num+fid); push!(J,num+tid); push!(V,-1)
        push!(I,num+tid); push!(J,num+fid); push!(V,-1)
    end
    A = sparse(I,J,V)
    for i in 1:(2*num)
        push!(I,i); push!(J,i)
        push!(V,-sum(A[i,:])+1)
    end
    A = sparse(I,J,V,2*num,2*num)
    # λ, ϕ = eigs(A)
    # @show λ
    C = sparse(cholfact(A))
    return C
end

function KiwielRhoUpdate(opfdata,mpsoln,sscval,vval,agg_norm,epshat,rho,tVal,v_est,ssc_cntr)
  global ssc
  nextTVal=tVal
  oldTVal=tVal
  if sscval >= ssc
    if sscval >= 0.5 && ssc_cntr > 0
      nextTVal = 2*tVal*(1-sscval)
    elseif ssc_cntr > 3
      nextTVal = tVal/2
    end
    tVal = max(nextTVal, tVal/10, tMin)
    v_est = max(v_est,2*vval)
    ssc_cntr = max(ssc_cntr+1,1)
    if oldTVal != tVal
      ssc_cntr = 1
    end
  else
    v_est=min(v_est,agg_norm+epshat)
#@show rho*mpsoln.linerr,max(v_est,10*vval),ssc_cntr
    if rho*mpsoln.linerr > max(v_est,10*vval) && ssc_cntr < -3
      nextTVal = 2*tVal*(1-sscval)
    end
    tVal = min(nextTVal,10*tVal)
    ssc_cntr = min(ssc_cntr-1,-1)
    if oldTVal != tVal
      ssc_cntr = -1	
    end
  end
  return tVal,v_est,ssc_cntr
end

function update_agg(opfdata,bundles,ctr,mpsoln,sg_agg,ctr_bundles,agg_bundles,agg_bundle)
  nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
  N, L, G = opfdata.N, opfdata.L, opfdata.G 
  comp_agg(opfdata,ctr.soln,mpsoln.soln,sg_agg)
  agg_norm=comp_norm(opfdata,sg_agg)
  rho = 0 #ctr.cut_dual
  ncuts,naggcuts=length(bundles),length(agg_bundles)
  if ncuts > 0
    #@show sum(bundles[n].cut_dual for n in 1:ncuts) 
    rho += sum(bundles[n].cut_dual for n in 1:ncuts) 
  end
  if naggcuts > 0
    #@show sum(agg_bundles[n].cut_dual for n in 1:naggcuts)
    rho += sum(agg_bundles[n].cut_dual for n in 1:naggcuts)
  end
  epshat = mpsoln.linobjval - (ctr.linobjval - rho*ctr.eta)
  - (dot(sg_agg.α[N],(mpsoln.soln.α[N]-ctr.soln.α[N])) - dot(sg_agg.β[N],(mpsoln.soln.β[N]-ctr.soln.β[N])) 
  + dot(sg_agg.γ[N],(mpsoln.soln.γ[N]-ctr.soln.γ[N])) - dot(sg_agg.δ[N],(mpsoln.soln.δ[N]-ctr.soln.δ[N])) 
  + dot(sg_agg.λF[L],(mpsoln.soln.λF[L] - ctr.soln.λF[L])) - dot(sg_agg.λT[L],(mpsoln.soln.λT[L] - ctr.soln.λT[L]))
  + dot(sg_agg.x[L], (mpsoln.soln.x[L] - ctr.soln.x[L]))
  + dot(sg_agg.μF[L],(mpsoln.soln.μF[L] - ctr.soln.μF[L])) - dot(sg_agg.μT[L],(mpsoln.soln.μT[L] - ctr.soln.μT[L])) )
  return agg_norm,epshat,rho
end

function updateCenter(opfdata,mpsoln,ctr,trl_bundles,ctr_bundles,agg_bundles)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC
  cpy_bundle(opfdata,mpsoln,ctr)
  for n=1:length(trl_bundles)
    updateLinErr(opfdata,ctr,trl_bundles[n])
  end
  for n=1:length(ctr_bundles)
    updateLinErr(opfdata,ctr,ctr_bundles[n])
  end
  for n=1:length(agg_bundles)
    updateLinErr(opfdata,ctr,agg_bundles[n])
  end
end

function updateLinErr(opfdata,ctr,bundle)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  bundle.linerr = ctr.eta - bundle.eta 
    - dot(bundle.eta_sg.α[N],(bundle.soln.α[N] - ctr.soln.α[N]))
    - dot(bundle.eta_sg.β[N],(bundle.soln.β[N] - ctr.soln.β[N]))
    - dot(bundle.eta_sg.γ[N],(bundle.soln.γ[N] - ctr.soln.γ[N]))
    - dot(bundle.eta_sg.δ[N],(bundle.soln.δ[N] - ctr.soln.δ[N]))
    - dot(bundle.eta_sg.λF[L],(bundle.soln.λF[L] - ctr.soln.λF[L]))
    - dot(bundle.eta_sg.λT[L],(bundle.soln.λT[L] - ctr.soln.λT[L]))
    - dot(bundle.eta_sg.μF[L],(bundle.soln.μF[L] - ctr.soln.μF[L]))
    - dot(bundle.eta_sg.μT[L],(bundle.soln.μT[L] - ctr.soln.μT[L]))
end

function computeSG(opfdata,mpsoln)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC
      vR = zeros(nbuses)
      vI = zeros(nbuses)
      try
        mpsoln.eta = -solveEta0Eigs(opfdata,mpsoln.soln,vR,vI)
      catch exc
        println("Exception caught with eigs(), solving η0Val subproblem with Ipopt as recourse.")
        println(exc)
        mpsoln.eta = -solveEta0SDP(opfdata,mpsoln.soln,vR,vI)
      end
      for i in N
        W_val = vR[i]^2 + vI[i]^2
        mpsoln.eta_sg.α[i] = Y["shR"][i] * W_val
	mpsoln.eta_sg.β[i] = -Y["shI"][i] * W_val
        mpsoln.eta_sg.δ[i] = W_val
	mpsoln.eta_sg.γ[i] = -W_val
      end
      for l in L
        from = fromBus[l]; to = toBus[l]
        e_valF = vR[from]; f_valF = vI[from]; W_valF = e_valF^2 + f_valF^2
        e_valT = vR[to]; f_valT = vI[to]; W_valT = e_valT^2 + f_valT^2
        Wr_val = e_valF*e_valT + f_valF*f_valT; Wi_val = e_valT*f_valF - e_valF*f_valT
        mpsoln.eta_sg.λF[l] = (Y["ffR"][l] * W_valF + Y["ftR"][l] * Wr_val + Y["ftI"][l] * Wi_val)
        mpsoln.eta_sg.λT[l] = (Y["ttR"][l] * W_valT + Y["tfR"][l] * Wr_val - Y["tfI"][l] * Wi_val)
        mpsoln.eta_sg.μF[l] = (-Y["ffI"][l] * W_valF - Y["ftI"][l] * Wr_val + Y["ftR"][l] * Wi_val)
        mpsoln.eta_sg.μT[l] = (-Y["ttI"][l] * W_valT - Y["tfI"][l] * Wr_val - Y["tfR"][l] * Wi_val)
      end
      return mpsoln.eta
end
# SUBROUTINE FOR COMPUTING THE MINIMUM EIGENVALUE OF H WITH A CORRESPONDING EIGENVECTOR
function solveEta0Eigs(opfdata,soln,vR,vI)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus = opfdata.fromBus, opfdata.toBus
      H=spzeros(2*nbuses,2*nbuses)
      updateHess(opfdata,soln,H)
      E=eigs(H,nev=1,which=:SR, maxiter=100000, tol=1e-6)
      η0Val = E[1][1]
      for i in N
        vR[i] = E[2][i,1]; vI[i] = E[2][nbuses+i,1]
      end
      return η0Val
end
# Update Hessian
function updateHess(opfdata,pi_val,H)
      #lines, buses, generators, baseMVA = opfdata.lines, opfdata.buses, opfdata.generators, opfdata.baseMVA
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus = opfdata.fromBus, opfdata.toBus  
      Y = opfdata.Y_AC
      for i in N
        H[i,i] +=  pi_val.α[i] * Y["shR"][i] - pi_val.β[i] * Y["shI"][i]  + pi_val.δ[i] - pi_val.γ[i]
        H[nbuses+i,nbuses+i] += pi_val.α[i] * Y["shR"][i] - pi_val.β[i] * Y["shI"][i] + pi_val.δ[i] - pi_val.γ[i]
      end
      for l in L
        from = fromBus[l]; to = toBus[l]
        H[from,from] += pi_val.λF[l] * Y["ffR"][l] - pi_val.μF[l] * Y["ffI"][l]
        H[nbuses+from,nbuses+from] += pi_val.λF[l] * Y["ffR"][l] - pi_val.μF[l] * Y["ffI"][l]
        H[to,to] += pi_val.λT[l] * Y["ttR"][l] - pi_val.μT[l] * Y["ttI"][l]
        H[nbuses+to,nbuses+to] += pi_val.λT[l] * Y["ttR"][l] - pi_val.μT[l] * Y["ttI"][l]
        H[from,to] += 0.5*( pi_val.λF[l] * Y["ftR"][l] - pi_val.μF[l] * Y["ftI"][l] + pi_val.λT[l] * Y["tfR"][l] - pi_val.μT[l] * Y["tfI"][l] )
        H[to,from] += 0.5*( pi_val.λF[l] * Y["ftR"][l] - pi_val.μF[l] * Y["ftI"][l] + pi_val.λT[l] * Y["tfR"][l] - pi_val.μT[l] * Y["tfI"][l] )
        H[nbuses+from, nbuses+to] += 0.5*( pi_val.λF[l] * Y["ftR"][l] - pi_val.μF[l] * Y["ftI"][l] + pi_val.λT[l] * Y["tfR"][l] - pi_val.μT[l] * Y["tfI"][l] )
        H[nbuses+to, nbuses+from] += 0.5*( pi_val.λF[l] * Y["ftR"][l] - pi_val.μF[l] * Y["ftI"][l] + pi_val.λT[l] * Y["tfR"][l] - pi_val.μT[l] * Y["tfI"][l] )
        H[to, nbuses+from] += 0.5*( pi_val.λF[l] * Y["ftI"][l] - pi_val.λT[l] * Y["tfI"][l] + pi_val.μF[l] * Y["ftR"][l] - pi_val.μT[l] * Y["tfR"][l] )
        H[nbuses+from, to] += 0.5*( pi_val.λF[l] * Y["ftI"][l] - pi_val.λT[l] * Y["tfI"][l] + pi_val.μF[l] * Y["ftR"][l] - pi_val.μT[l] * Y["tfR"][l] )
        H[from,nbuses+to] -= 0.5*( pi_val.λF[l] * Y["ftI"][l] - pi_val.λT[l] * Y["tfI"][l] + pi_val.μF[l] * Y["ftR"][l] - pi_val.μT[l] * Y["tfR"][l] )
        H[nbuses+to,from] -= 0.5*( pi_val.λF[l] * Y["ftI"][l] - pi_val.λT[l] * Y["tfI"][l] + pi_val.μF[l] * Y["ftR"][l] - pi_val.μT[l] * Y["tfR"][l] )
      end
end


# SUBROUTINE FOR COMPUTING THE MINIMUM EIGENVALUE OF H WITH A CORRESPONDING EIGENVECTOR
  # VIA AN OPTIMIZATION PROBLEM
function solveEta0SDP(opfdata,soln,vR,vI)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus = opfdata.fromBus, opfdata.toBus

      #The QP subproblem
      mSDP = Model(solver=IpoptSolver())
      @variable(mSDP, e[i=N], start=0); @variable(mSDP, f[i=N], start=0)
      η0Val = 0

      for i in N
        setvalue(e[i], 1); setvalue(f[i], 0)
      end

      H=spzeros(2*nbuses,2*nbuses)
      updateHess(opfdata,soln,H)

      # Adjust QP subproblem
      @NLobjective(mSDP, Min, sum( H[i,i]*(e[i]^2+f[i]^2) for i in N)
        + 2*sum( H[fromBus[l],toBus[l]]*(e[fromBus[l]]*e[toBus[l]]+f[fromBus[l]]*f[toBus[l]])   for l in L)
        - 2*sum( H[fromBus[l],nbuses+toBus[l]]*(f[fromBus[l]]*e[toBus[l]]-e[fromBus[l]]*f[toBus[l]])   for l in L)
      )
      status = solve(mSDP)
      if status == :Optimal || status == :UserLimit
        η0Val = getobjectivevalue(mSDP)
        for i in N
          vR[i],vI[i]=getvalue(e[i]),getvalue(f[i])
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

# SUBROUTINE FOR COMPUTING A SUBGRADIENT OF ETA(PI), WHICH IS THE FUNCTION TAKING THE VALUE OF THE MINIMUM EIGENVALUE OF H(PI)
function purgeSG(opfdata,bundle,minAge=20,maxAge=80)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC

      orig_ncuts = length(bundle)
      ncuts = orig_ncuts
      for n=orig_ncuts:-1:1
        if (abs(bundle[n].cut_dual) < 1e-8 && bundle[n].age > minAge) || bundle[n].age > maxAge
	  bundle[n]=bundle[ncuts]
	  delete!(bundle,ncuts)
	  ncuts -= 1
	end
      end
      for n=1:length(bundle)
	bundle[n].age += 1
      end
      return length(bundle)
end

function aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC
      agg_bundle=create_bundle(opfdata)

      cpy_soln(opfdata,mpsoln.soln,agg_bundle.soln)
      sumDuals,ncuts,nctrcuts,naggcuts = 0.0,length(trl_bundles),length(ctr_bundles),length(agg_bundles)
      if ncuts > 0
        sumDuals += sum(trl_bundles[n].cut_dual for n in 1:ncuts) 
      end
      if nctrcuts > 0
        sumDuals += sum(ctr_bundles[n].cut_dual for n in 1:nctrcuts) 
      end
      if naggcuts > 0
        sumDuals += sum(agg_bundles[n].cut_dual for n in 1:naggcuts) 
      end

      agg_bundle.eta_sg.α[N] .= 0
      agg_bundle.eta_sg.β[N] .= 0
      agg_bundle.eta_sg.γ[N] .= 0
      agg_bundle.eta_sg.δ[N] .= 0
      agg_bundle.eta_sg.ζpLB[G] .= 0
      agg_bundle.eta_sg.ζpUB[G] .= 0
      agg_bundle.eta_sg.ζqLB[G] .= 0
      agg_bundle.eta_sg.ζqUB[G] .= 0
      agg_bundle.eta_sg.x[L] .= 0
      agg_bundle.eta_sg.λF[L] .= 0
      agg_bundle.eta_sg.λT[L] .= 0
      agg_bundle.eta_sg.μF[L] .= 0
      agg_bundle.eta_sg.μT[L] .= 0

      agg_bundle.objval = 0
      agg_bundle.linobjval = 0
      agg_bundle.penval = 0
      agg_bundle.psival = 0

      if ncuts > 0
        agg_bundle.eta_sg.α[N] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.α[N] for n in 1:ncuts) )
        agg_bundle.eta_sg.β[N] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.β[N] for n in 1:ncuts) )
        agg_bundle.eta_sg.γ[N] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.γ[N] for n in 1:ncuts) )
        agg_bundle.eta_sg.δ[N] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.δ[N] for n in 1:ncuts) )
        agg_bundle.eta_sg.ζpLB[G] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.ζpLB[G] for n in 1:ncuts) )
        agg_bundle.eta_sg.ζpUB[G] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.ζpUB[G] for n in 1:ncuts) )
        agg_bundle.eta_sg.ζqLB[G] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.ζqLB[G] for n in 1:ncuts) )
        agg_bundle.eta_sg.ζqUB[G] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.ζqUB[G] for n in 1:ncuts) )
        agg_bundle.eta_sg.x[L]  += (1/sumDuals)* ( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.x[L] for n in 1:ncuts) )
        agg_bundle.eta_sg.λF[L] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.λF[L] for n in 1:ncuts) )
        agg_bundle.eta_sg.λT[L] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.λT[L] for n in 1:ncuts) )
        agg_bundle.eta_sg.μF[L] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.μF[L] for n in 1:ncuts) )
        agg_bundle.eta_sg.μT[L] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.μT[L] for n in 1:ncuts) )

        agg_bundle.objval += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].objval for n in 1:ncuts) )
        agg_bundle.linobjval += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].linobjval for n in 1:ncuts) )
        agg_bundle.penval += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].penval for n in 1:ncuts) )
        agg_bundle.psival += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].psival for n in 1:ncuts) )

        agg_bundle.etahat -= (1/sumDuals)*sum( trl_bundles[n].cut_dual*(trl_bundles[n].eta 
	+ dot(trl_bundles[n].eta_sg.α[N],trl_bundles[n].soln.α[N]) + dot(trl_bundles[n].eta_sg.β[N],trl_bundles[n].soln.β[N])
	  + dot(trl_bundles[n].eta_sg.γ[N],trl_bundles[n].soln.γ[N]) + dot(trl_bundles[n].eta_sg.δ[N],trl_bundles[n].soln.δ[N]) 
        + dot(trl_bundles[n].eta_sg.λF[L],trl_bundles[n].soln.λF[L]) + dot(trl_bundles[n].eta_sg.λT[L],trl_bundles[n].soln.λT[L]) 
	  + dot(trl_bundles[n].eta_sg.μF[L],trl_bundles[n].soln.μF[L]) + dot(trl_bundles[n].eta_sg.μT[L],trl_bundles[n].soln.μT[L]) ) for n in 1:ncuts)
      end
      if nctrcuts > 0
        agg_bundle.eta_sg.α[N] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.α[N] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.β[N] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.β[N] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.γ[N] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.γ[N] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.δ[N] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.δ[N] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.ζpLB[G] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.ζpLB[G] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.ζpUB[G] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.ζpUB[G] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.ζqLB[G] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.ζqLB[G] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.ζqUB[G] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.ζqUB[G] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.x[L]  += (1/sumDuals)* ( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.x[L] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.λF[L] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.λF[L] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.λT[L] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.λT[L] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.μF[L] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.μF[L] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.μT[L] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.μT[L] for n in 1:nctrcuts) )

        agg_bundle.objval += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].objval for n in 1:nctrcuts) )
        agg_bundle.linobjval += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].linobjval for n in 1:nctrcuts) )
        agg_bundle.penval += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].penval for n in 1:nctrcuts) )
        agg_bundle.psival += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].psival for n in 1:nctrcuts) )

        agg_bundle.etahat -= (1/sumDuals)*sum( ctr_bundles[n].cut_dual*(ctr_bundles[n].eta 
	+ dot(ctr_bundles[n].eta_sg.α[N],ctr_bundles[n].soln.α[N]) + dot(ctr_bundles[n].eta_sg.β[N],ctr_bundles[n].soln.β[N])
	  + dot(ctr_bundles[n].eta_sg.γ[N],ctr_bundles[n].soln.γ[N]) + dot(ctr_bundles[n].eta_sg.δ[N],ctr_bundles[n].soln.δ[N]) 
        + dot(ctr_bundles[n].eta_sg.λF[L],ctr_bundles[n].soln.λF[L]) + dot(ctr_bundles[n].eta_sg.λT[L],ctr_bundles[n].soln.λT[L]) 
	  + dot(ctr_bundles[n].eta_sg.μF[L],ctr_bundles[n].soln.μF[L]) + dot(ctr_bundles[n].eta_sg.μT[L],ctr_bundles[n].soln.μT[L]) ) for n in 1:nctrcuts)
      end
      if naggcuts > 0
        agg_bundle.eta_sg.α[N] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.α[N] for n in 1:naggcuts) )
        agg_bundle.eta_sg.β[N] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.β[N] for n in 1:naggcuts) )
        agg_bundle.eta_sg.γ[N] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.γ[N] for n in 1:naggcuts) )
        agg_bundle.eta_sg.δ[N] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.δ[N] for n in 1:naggcuts) )
        agg_bundle.eta_sg.ζpLB[G] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.ζpLB[G] for n in 1:naggcuts) )
        agg_bundle.eta_sg.ζpUB[G] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.ζpUB[G] for n in 1:naggcuts) )
        agg_bundle.eta_sg.ζqLB[G] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.ζqLB[G] for n in 1:naggcuts) )
        agg_bundle.eta_sg.ζqUB[G] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.ζqUB[G] for n in 1:naggcuts) )
        agg_bundle.eta_sg.x[L]  += (1/sumDuals)* ( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.x[L] for n in 1:naggcuts) )
        agg_bundle.eta_sg.λF[L] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.λF[L] for n in 1:naggcuts) )
        agg_bundle.eta_sg.λT[L] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.λT[L] for n in 1:naggcuts) )
        agg_bundle.eta_sg.μF[L] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.μF[L] for n in 1:naggcuts) )
        agg_bundle.eta_sg.μT[L] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.μT[L] for n in 1:naggcuts) )

        agg_bundle.objval += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].objval for n in 1:naggcuts) )
        agg_bundle.linobjval += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].linobjval for n in 1:naggcuts) )
        agg_bundle.penval += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].penval for n in 1:naggcuts) )
        agg_bundle.psival += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].psival for n in 1:naggcuts) )

        agg_bundle.etahat -= (1/sumDuals)*sum( agg_bundles[n].cut_dual*(agg_bundles[n].eta 
	+ dot(agg_bundles[n].eta_sg.α[N],agg_bundles[n].soln.α[N]) + dot(agg_bundles[n].eta_sg.β[N],agg_bundles[n].soln.β[N])
	  + dot(agg_bundles[n].eta_sg.γ[N],agg_bundles[n].soln.γ[N]) + dot(agg_bundles[n].eta_sg.δ[N],agg_bundles[n].soln.δ[N]) 
        + dot(agg_bundles[n].eta_sg.λF[L],agg_bundles[n].soln.λF[L]) + dot(agg_bundles[n].eta_sg.λT[L],agg_bundles[n].soln.λT[L]) 
	  + dot(agg_bundles[n].eta_sg.μF[L],agg_bundles[n].soln.μF[L]) + dot(agg_bundles[n].eta_sg.μT[L],agg_bundles[n].soln.μT[L]) ) for n in 1:naggcuts)
      end
      agg_bundle.age = 1
      agg_bundle.eta = mpsoln.eta
      return agg_bundle
end

