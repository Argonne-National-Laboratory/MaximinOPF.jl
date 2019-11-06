#=
Kibaek Kim
Brian Dandurand
=#

include("utils.jl")
include("MP.jl")

#Implements the level bundle method approach of De Oliveira 2016
function PLBM_Oliveira(opfdata,params,K,HEUR)
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


