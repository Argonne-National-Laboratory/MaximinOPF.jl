#=
Kibaek Kim
Brian Dandurand
=#

include("utils.jl")
include("MP.jl")

#Based on Fuduli and Gaudioso 2006
function PBM_FuduliGaudioso(opfdata,params,K,HEUR)
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

