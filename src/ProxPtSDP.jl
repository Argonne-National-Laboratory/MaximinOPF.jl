#=
Kibaek Kim
Brian Dandurand
=#

include("utils.jl")
include("MP.jl")

function testProxTraj(opfdata,params,K,HEUR)
    TOL = 1e-5
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # INITIAL ITERATION
    bundles=Dict()
    ncuts=0
    params.rho = 0
    mpsoln=create_bundle(opfdata)
    ctr=create_bundle(opfdata)
    agg_bundles=Dict()
    agg_bundle=create_bundle(opfdata)
    sg_agg=create_soln(opfdata)


    v_est,ssc_cntr,agg_norm = 1e20,0,0.0
  # MAIN LOOP
    for kk=1:params.maxNSG
     # PHASE I
      tMax = max(tMin,0.5*(comp_norm(opfdata,ctr.eta_sg)^2)/(1+abs(ctr.eta)))
      println("Penalty will be in [",tMin,",",tMax,"].")
      params.tVal = (tMax+tMin)/2
     # PHASE II
      while true
     # STEP 1
        mpsoln=create_bundle(opfdata)
        while mpsoln.status != MOI.OPTIMAL && mpsoln.status != MOI.LOCALLY_SOLVED
	  params.tVal /= 2
          println("Resolving with reduced prox parameter value: ",params.tVal)
          status = solveNodeMP(opfdata,node_data,params,bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
        end	
        if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
	  ncuts=length(bundles)
	  naggcuts=length(bundles)
	  if ncuts > 0
	    params.rho = sum(bundles[n].cut_dual for n in 1:ncuts) + sum(agg_bundles[n].cut_dual for n in 1:naggcuts) 
           # COMPUTING LINEARIZATION ERRORS
	    agg_norm,epshat,params.rho=update_agg(opfdata,bundles,ctr,mpsoln,sg_agg,ctr_bundles,agg_bundles,agg_bundle)
            if params.rho*mpsoln.eta <= params.ssc*(1/params.tVal)*agg_norm^2 
	      break
	    end
	  end
          agg_bundles[1]=agg_bundle
	  ncuts=purgeSG(opfdata,bundles)
          bundles[ncuts+1]=mpsoln
        else
	  println("Solver returned: $status")
        end
      end
      updateCenter(opfdata,mpsoln,ctr,bundles)
      @show "Traj",kk,mpsoln.linobjval,mpsoln.eta,ncuts,agg_norm,params
      if agg_norm < 1e-2
        println("Convergence to within tolerance: obj, feas ",mpsoln.linobjval," ",mpsoln.eta)
	break
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

function testProxPt0(opfdata,params,K,HEUR,node_data)
    println("Applying the algorithm of Delfino and de Oliveira 2018....")
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # INITIAL ITERATION
    trl_bundles=Dict()
    ctr_bundles=Dict()
    agg_bundles=Dict()
    ncuts=0
    params.rho,params.rhoUB = 0,0
    mpsoln=create_bundle(opfdata)
    ctr=create_bundle(opfdata)
    status = solveNodeMP(opfdata,node_data,params,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
    while mpsoln.status != MOI.OPTIMAL && mpsoln.status != MOI.LOCALLY_SOLVED
      params.tVal /= 2
      println("Status was: ",mpsoln.status,". Resolving with reduced prox parameter value: ",params.tVal)
      status = solveNodeMP(opfdata,node_data,params,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
    end	
    updateCenter(opfdata,mpsoln,ctr,trl_bundles,ctr_bundles,agg_bundles)
    ctr_bundles[1]=mpsoln

    agg_bundle=create_bundle(opfdata)
    sg_agg=create_soln(opfdata)

    v_est,ssc_cntr = 1e20,0
    tL,tU=params.tMin,params.tMax
    TOL = 1e-5
  # MAIN LOOP
    for kk=1:params.maxNSG
     # STEP 1
      mpsoln=create_bundle(opfdata)
      status = solveNodeMP(opfdata,node_data,params,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
      while mpsoln.status != MOI.OPTIMAL && mpsoln.status != MOI.LOCALLY_SOLVED
	params.tVal /= 2
        println("Status was: ",mpsoln.status,". Resolving with reduced prox parameter value: ",params.tVal)
        status = solveNodeMP(opfdata,node_data,params,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
      end	
      if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
       # STEP 2
        ncuts,nctrcuts,naggcuts = length(trl_bundles),length(ctr_bundles),length(agg_bundles)
	 # UPDATING RHO AS NECESSARY TO CORRESPOND TO EXACT PENALTY
         # COMPUTING AGGREGATION INFORMATION
          agg_norm = update_agg(opfdata,params,ctr,mpsoln,sg_agg)
          update_rho(params,trl_bundles,ctr_bundles,agg_bundles)
          epshat = compute_epshat(opfdata,params,mpsoln,ctr,sg_agg)
	  if params.rhoUB < params.rho
	    params.rhoUB = params.rho + 1
	  end
          agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
          if ctr.eta < TOL && agg_norm < 1e-3 && epshat < 1e-3 
	    println("Convergence to within tolerance: ")
	    @show kk,ncuts,ssc_cntr,params.tVal,params.rho
	    @show ctr.linobjval,ctr.eta,agg_norm,epshat
	    break
          end
         # STEP 3
	  sscval = ((mpsoln.linobjval - params.rhoUB*mpsoln.eta)-(ctr.linobjval - params.rhoUB*ctr.eta))/(mpsoln.linobjval-(ctr.linobjval - params.rhoUB*ctr.eta)) 
	  vval = (mpsoln.linobjval-(ctr.linobjval - params.rhoUB*ctr.eta)) 
          if sscval >= params.ssc 
         # UPDATE CENTER VALUES
	    nctrcuts=purgeSG(opfdata,ctr_bundles)
	    ctr_bundles[nctrcuts+1]=mpsoln
	    updateCenter(opfdata,mpsoln,ctr,trl_bundles,ctr_bundles,agg_bundles)
	    if mpsoln.eta < 1e-5
	      params.tVal *= 0.8
	    end
	    @show kk,ncuts,ssc_cntr,params.tVal,params.rho
	    @show mpsoln.linobjval,mpsoln.eta,agg_norm,epshat
	  else
	    ncuts=purgeSG(opfdata,trl_bundles)
            trl_bundles[ncuts+1]=mpsoln
          end
      else
	println("Solver returned: $status")
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

# Implements the approach of Sagastizabal and Solodov 2005
function testProxPt(opfdata,params,K,HEUR,node_data)
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # INITIAL ITERATION
    bundles=Dict()
    mpsoln=create_bundle(opfdata)
    mpsoln.soln.x[L] = x_val[L]
    sg_agg=create_soln(opfdata)
    ctr=mpsoln

    v_est,ssc_cntr = 1e20,0
  # MAIN LOOP
    for kk=1:params.maxNSG
     # STEP 1
      mpsoln=create_bundle(opfdata)
      status = solveNodeMP(opfdata,node_data,params,bundles,K,HEUR,ctr,PROX,mpsoln)
      while status != :Optimal
	params.tVal /= 2
        println("Resolving with reduced prox parameter value: ",params.tVal)
        status = solveNodeMP(opfdata,node_data,params,bundles,K,HEUR,ctr,PROX,mpsoln)
      end	
      if status == :Optimal 
	comp_agg(opfdata,ctr.soln,mpsoln.soln,sg_agg)
	agg_norm=comp_norm(opfdata,sg_agg)
	epshat=max(0,ctr.eta)-mpsoln.psival-(1.0/params.tVal)*agg_norm^2
	#del=epshat+0.5*(1.0/params.tVal)*agg_norm^2
        del=ctr.eta-mpsoln.objval
       # STEP 2
        if del < 1e-6 
	  println("Convergence to within tolerance: del = ",del," val: ",mpsoln.linobjval," and feas: ",mpsoln.eta)
	  break
        end
       # STEP 3
	hk = max(ctr.linobjval-mpsoln.linobjval,mpsoln.eta)
	sscval = max(0,ctr.eta) - hk
        if hk <= max(0,ctr.eta) - params.ssc*del || kk==1
          # UPDATE CENTER VALUES
	    updateCenter(opfdata,mpsoln,ctr,bundles)
	    #purgeSG(opfdata,bundles)
	    oldncuts = length(bundles)
	    ncuts = oldncuts
	    for n=oldncuts:-1:1
	      if ctr.linobjval-bundles[n].linobjval > bundles[n].psival
		bundles[n]=bundles[ncuts]
	        delete!(bundles,ncuts)
		ncuts -= 1
	      end
	    end
	    if ctr.linobjval-mpsoln.linobjval > mpsoln.psival
	      ncuts = length(bundles)
              bundles[ncuts+1]=mpsoln
	    end
            @show "Sagadov",kk,ctr.linobjval,ctr.eta,epshat,del,hk,length(bundles),params
	else
	  ncuts = length(bundles)
          bundles[ncuts+1]=mpsoln
        end
       # STEP 4
        #tVal,v_est,ssc_cntr = KiwielRhoUpdate(opfdata,params,mpsoln,sscval,del,agg_norm,epshat,v_est,ssc_cntr)
      else
	#println("Solve status: $status")
	#break
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

#Implements the approach of De Oliveira 2016
function testLevelBM(opfdata,params,K,HEUR)
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # INITIAL ITERATION
    optUB=1e20
    # INITIAL SOLN (ALL ZEROS)
      mpsoln=create_bundle(opfdata)
    bestsoln = mpsoln
    ctr=mpsoln 
    bundles = Dict()
    prox_term=LVLINF

    cpsoln=create_bundle(opfdata)
    solveNodeMP(opfdata,node_data,params,bundles,K,HEUR,ctr,CP,cpsoln)
    if cpsoln.status == :Optimal 
	optUB = cpsoln.objval
        hkctr = max(optUB,0)
    end

  # MAIN LOOP
    hkctr_updated = false
    optub_updated = false
    for kk=1:params.maxNSG
     # STEP 1
      hkval = max(optUB - ctr.linobjval,ctr.eta)
      bestsoln = ctr
      bestIdx = 0
      ncuts=length(bundles)
      for pp=ncuts:-1:1
	if max(optUB - bundles[pp].linobjval,bundles[pp].eta) < hkval
	  hkval = max(optUB - bundles[pp].linobjval,bundles[pp].eta) 
          bestsoln = bundles[pp]
	  bestIdx=pp
	end
      end
      node_data.nodeBd = optUB - params.ssc*hkval
      if hkval < 1e-6
	println("Convergence to within tolerance.")
	break
      end

     # STEP 2
        if hkval <= (1-params.ssc)*hkctr 
          hkctr = hkval
          hkctr_updated = true
          # UPDATE CENTER VALUES
	    if bestIdx > 0
	      updateCenter(opfdata,bestsoln,ctr,bundles)
	    end
	  @show "hk  update",optUB,ctr.linobjval,ctr.eta,hkval,ncuts
        end

     # STEP 3
      cpstatus = solveNodeMP(opfdata,node_data,params,bundles,K,HEUR,ctr,CP,cpsoln)
      if cpsoln.status != :Optimal
	println("FLAG! cpsoln.status = ",cpsoln.status)
      end
      if node_data.nodeBd-cpsoln.linobjval <= 0.0
	params.tVal = 0.1
        mpsoln=create_bundle(opfdata)
        status=solveNodeMP(opfdata,node_data,params,bundles,K,HEUR,ctr,prox_term,mpsoln)
        while status != :Optimal
	  params.tVal /= 2
          println("Resolving with reduced prox parameter value: ",params.tVal)
          status=solveNodeMP(opfdata,node_data,params,bundles,K,HEUR,ctr,prox_term,mpsoln)
        end	
	if mpsoln.status != :Optimal
	  println("Taking recourse since mpsoln does not have an optimal solution for MP. Feas: ",node_data.nodeBd - mpsoln.linobjval," val: ",node_data.nodeBd - cpsoln.linobjval)
	  cpy_soln(opfdata,cpsoln,mpsoln)
	else
	  if hkctr_updated
	    oldncuts = ncuts
	    ncuts=purgeSG(opfdata,bundles)
	    #println("There were ",oldncuts," cuts; after purging, there are ",ncuts)
	    hkctr_updated = false
	  end
        end

        if mpsoln.eta > 0
          bundles[ncuts+1] = mpsoln
	  ncuts += 1
        else
            println("Tolerance met for not generating a new lazy cut, obj,eta=",mpsoln.linobjval,mpsoln.eta,".")
	    if mpsoln.linobjval < optUB
	      optUB = mpsoln.linobjval
              optub_updated = true
	    end
	    @show "bnd opt update",optUB,mpsoln.linobjval,-mpsoln.eta,hkval,ncuts
	end
      else
	  optUB = node_data.nodeBd 
          optub_updated = true
          hkctr = max(optUB - bestsoln.linobjval,bestsoln.eta)
          hkctr_updated = true
	  if bestIdx > 0
	    updateCenter(opfdata,bestsoln,ctr,bundles)
	  end
	  @show "bnd update",kk,optUB,bestsoln.linobjval,bestsoln.eta,hkval,ncuts
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end


