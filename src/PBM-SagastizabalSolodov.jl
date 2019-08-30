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

