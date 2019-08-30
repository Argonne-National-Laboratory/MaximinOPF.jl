
# Defining datatypes
mutable struct NodeInfo
  x_lbs::Array{Float64}
  x_ubs::Array{Float64}
  nodeBd::Float64
end
mutable struct Params
  ALG::Int
  maxNSG::Int
  tMin::Float64
  tMax::Float64
  tVal::Float64
  rho::Float64
  rhoUB::Float64
  ssc::Float64
  ssc_cntr::Int
  tol::Float64
end
mutable struct Soln
  α::Array{Float64}
  β::Array{Float64}
  γ::Array{Float64}
  δ::Array{Float64}
  ζpLB::Array{Float64}
  ζpUB::Array{Float64}
  ζqLB::Array{Float64}
  ζqUB::Array{Float64}
  x::Array{Float64}
  λF::Array{Float64}
  λT::Array{Float64}
  μF::Array{Float64}
  μT::Array{Float64}
end
function create_soln(opfdata)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  return Soln(zeros(nbuses),zeros(nbuses),zeros(nbuses),zeros(nbuses),
    zeros(ngens),zeros(ngens),zeros(ngens),zeros(ngens),
    zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines))
end
function cpy_soln(opfdata,fromSoln,toSoln)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  toSoln.α[N] = fromSoln.α[N]
  toSoln.β[N] = fromSoln.β[N]
  toSoln.γ[N] = fromSoln.γ[N]
  toSoln.δ[N] = fromSoln.δ[N]
  toSoln.ζpLB[G] = fromSoln.ζpLB[G]
  toSoln.ζpUB[G] = fromSoln.ζpUB[G]
  toSoln.ζqLB[G] = fromSoln.ζqLB[G]
  toSoln.ζqUB[G] = fromSoln.ζqUB[G]
  toSoln.x[L] = fromSoln.x[L]
  toSoln.λF[L] = fromSoln.λF[L]
  toSoln.λT[L] = fromSoln.λT[L]
  toSoln.μF[L] = fromSoln.μF[L]
  toSoln.μT[L] = fromSoln.μT[L]
end

function comp_agg(opfdata,params,ctr,trl,agg)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  agg.α[N] = params.tVal*(ctr.α[N]-trl.α[N])
  agg.β[N] = params.tVal*(ctr.β[N]-trl.β[N])
  agg.γ[N] = params.tVal*(ctr.γ[N]-trl.γ[N])
  agg.δ[N] = params.tVal*(ctr.δ[N]-trl.δ[N])
  agg.ζpLB[G] = params.tVal*(ctr.ζpLB[G]-trl.ζpLB[G])
  agg.ζpUB[G] = params.tVal*(ctr.ζpUB[G]-trl.ζpUB[G])
  agg.ζqLB[G] = params.tVal*(ctr.ζqLB[G]-trl.ζqLB[G])
  agg.ζqUB[G] = params.tVal*(ctr.ζqUB[G]-trl.ζqUB[G])
  agg.x[L] = params.tVal*(ctr.x[L]-trl.x[L])
  agg.λF[L] = params.tVal*(ctr.λF[L]-trl.λF[L])
  agg.λT[L] = params.tVal*(ctr.λT[L]-trl.λT[L])
  agg.μF[L] = params.tVal*(ctr.μF[L]-trl.μF[L])
  agg.μT[L] = params.tVal*(ctr.μT[L]-trl.μT[L])
end
function comp_norm(opfdata,soln)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  return sqrt(norm(soln.α[N])^2 + norm(soln.β[N])^2 + norm(soln.γ[N])^2 + norm(soln.δ[N])^2 
	+ norm(soln.ζpLB[G])^2 + norm(soln.ζpUB[G])^2 + norm(soln.ζqLB[G])^2 + norm(soln.ζqUB[G])^2
  	+ norm(soln.x[L])^2 
	+ norm(soln.λF[L])^2 + norm(soln.λT[L])^2 + norm(soln.μF[L])^2 + norm(soln.μT[L])^2)
end


mutable struct Bundle
  soln::Soln
  eta_sg::Soln
  objval::Float64
  linobjval::Float64
  penval::Float64
  psival::Float64
  eta::Float64
  etahat::Float64
  linerr::Float64
  cut_dual::Float64
  lvl_dual::Float64
  age::Float64
  solvetime::Float64
  status
end
function create_bundle(opfdata)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  return Bundle(create_soln(opfdata),create_soln(opfdata),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0)
end
function cpy_bundle(opfdata,fromBundle,toBundle)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  cpy_soln(opfdata,fromBundle.soln,toBundle.soln)
  cpy_soln(opfdata,fromBundle.eta_sg,toBundle.eta_sg)
  toBundle.objval = fromBundle.objval
  toBundle.linobjval = fromBundle.linobjval
  toBundle.penval = fromBundle.penval
  toBundle.psival = fromBundle.psival
  toBundle.eta = fromBundle.eta
  toBundle.etahat = fromBundle.etahat
  toBundle.linerr = fromBundle.linerr
  toBundle.cut_dual = fromBundle.cut_dual
  toBundle.lvl_dual = fromBundle.lvl_dual
  toBundle.solvetime = fromBundle.solvetime
  toBundle.status = fromBundle.status
end

mutable struct ConstrDuals
  var_bds::Soln
  zeta_eq::Array{Float64,2}
  budget::Float64
  AMcf1::Array{Float64} 
  AMcf2::Array{Float64} 
  AMcf3::Array{Float64} 
  AMcf4::Array{Float64} 
  AMct1::Array{Float64} 
  AMct2::Array{Float64} 
  AMct3::Array{Float64} 
  AMct4::Array{Float64} 
  BMcf1::Array{Float64} 
  BMcf2::Array{Float64} 
  BMcf3::Array{Float64} 
  BMcf4::Array{Float64} 
  BMct1::Array{Float64} 
  BMct2::Array{Float64} 
  BMct3::Array{Float64} 
  BMct4::Array{Float64} 
end

function create_constr_duals(opfdata)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  return ConstrDuals(create_soln(opfdata),zeros(nbuses,ngens),0,
    zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),
    zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines))
end
