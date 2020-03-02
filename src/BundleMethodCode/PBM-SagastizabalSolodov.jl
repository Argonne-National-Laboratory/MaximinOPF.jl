#=
Kibaek Kim
Brian Dandurand
=#

include("utils.jl")
include("MP.jl")

# Implements the approach of Sagastizabal and Solodov 2005
function PBM_SagastizabalSolodov(opfdata,params,K,HEUR,node_data)
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = 1:nbuses,1:nlines,1:ngens 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # INITIAL ITERATION
    ctr_bundles=Dict()
    trl_bundles=Dict()
    agg_bundles=Dict()
    mpsoln=create_bundle(opfdata)
    sg_agg=create_soln(opfdata)
    ctr=mpsoln
    #mpsoln=computeMPSoln(opfdata,node_data,K,PROX,ctr,trl_bundles,ctr_bundles,agg_bundles)
    #ctr=mpsoln
    node_data.tVal = params.tStart
    tLow,tHigh=params.tMin,params.tMax
    last_kk=params.maxNIter
  # MAIN LOOP
    for kk=1:params.maxNIter
     # STEP 1
      mpsoln=computeMPSoln(opfdata,node_data,K,PROX,ctr,trl_bundles,ctr_bundles,agg_bundles)
       # STEP 2
        #if node_data.descent_est < 1e-4 
        if mpsoln.eta < params.tol1 && node_data.agg_sg_norm < params.tol2 && node_data.epshat < params.tol3 
	  println("Convergence to within tolerance: del = ",node_data.descent_est," val: ",mpsoln.linobjval," and feas: ",mpsoln.eta)
          last_kk=kk
	  break
        end
       # STEP 3
#@show node_data.descent,node_data.descent_est
        if node_data.descent >= params.ssc1*node_data.descent_est 
          if testSchrammZoweSSII(opfdata,params,node_data,mpsoln,ctr) || tHigh-tLow < 1e-4 
          # UPDATE CENTER VALUES
            agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
	    ctr=mpsoln
	    ntrlcuts,nnzcuts=purgeSG(opfdata,trl_bundles,params.maxNSG)
            trl_bundles[ntrlcuts+1]=mpsoln
            tLow,tHigh=params.tMin,params.tMax
	    @printf("ss: %d\t(objval,eta)=(%.4f,%.2e)\t(t,rho,ncuts)=(%.4f,%.3f,%d)\t(err,||s||,epshat,desc)=(%.2e,%.2e,%.2e,%.5e)\n",
	      kk,mpsoln.linobjval,mpsoln.eta,node_data.tVal,node_data.rho,node_data.ncuts,node_data.linerr,node_data.agg_sg_norm,node_data.epshat,node_data.descent_est)
            node_data.tVal = max(0.95*node_data.tVal,params.tMin)
	  else
            tHigh=node_data.tVal
	    node_data.tVal=2*tLow*tHigh/(tLow+tHigh)
	  end
	else
	  if testSchrammZoweNSII(opfdata,params,ctr,node_data,mpsoln) || tHigh-tLow < 1e-4
            agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
	    ntrlcuts,nnzcuts=purgeSG(opfdata,trl_bundles,params.maxNSG)
            trl_bundles[ntrlcuts+1]=mpsoln
            tLow,tHigh=params.tMin,params.tMax
	    @printf("ns: %d\t(objval,eta)=(%.4f,%.2e)\t(t,rho,ncuts)=(%.4f,%.3f,%d)\t(err,||s||,epshat,desc)=(%.2e,%.2e,%.2e,%.5e)\n",
	      kk,mpsoln.linobjval,mpsoln.eta,node_data.tVal,node_data.rho,node_data.ncuts,node_data.linerr,node_data.agg_sg_norm,node_data.epshat,node_data.descent_est)
	  else 
            tLow=node_data.tVal
	    node_data.tVal=2*tLow*tHigh/(tLow+tHigh)
	  end
        end
       # STEP 4
    end
    @printf("iter: %d\t(objval,eta)=(%.4f,%.2e)\t(t,rho,ncuts)=(%.4f,%.3f,%d)\t(err,||s||,epshat,desc)=(%.2e,%.2e,%.2e,%.5e)\n",
      last_kk,ctr.linobjval,ctr.eta,node_data.tVal,node_data.rho,node_data.ncuts,node_data.linerr,node_data.agg_sg_norm,node_data.epshat,node_data.descent_est)
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

