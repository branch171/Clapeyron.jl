"""
    GFEMTPFlash(;   max_trust_region_iters::Int = 0
    				tol::Float64 = 0.01*GFEM_tol
    				GFEM_tol::Float64 = sqrt(eps)
    				verbose = true
					flashin)

"""
Base.@kwdef struct GFEMTPFlash <: TPFlashMethod
    max_trust_region_iters::Int = 0
    tol::Float64 = 0.01*sqrt(eps(Float64))
    GFEM_tol::Float64 = sqrt(eps(Float64))
    verbose::Bool = false
	flashin::FlashResult
end

function tp_flash_impl(model::EoSModel, p, T, n, method::GFEMTPFlash)
	z₀ = n
	sumz₀ = sum(z₀)
	z₀ ./= sumz₀
	if method.max_trust_region_iters == 0
		max_trust_region_iters = 2000*length(n)
	else
		max_trust_region_iters = method.max_trust_region_iters
	end
	tol = method.tol
	GFEM_tol = method.GFEM_tol
	verbose = method.verbose
	flashin = method.flashin
	if verbose == true
		println("GFEM  - Setup:")
		println("GFEM  - trust region tolerence = $(tol)")
		println("GFEM  - GFEM tolerence = $(GFEM_tol)")
	end
	beta,xp,vp,Gsol = GFEM_impl(model,p,T,z₀,max_trust_region_iters,tol,GFEM_tol,verbose,flashin)
	return FlashResult(xp,beta,vp,FlashData(p,T,Gsol))
end

function GFEMConstraints_Gibbs(x,np,n₀,lb,ub,lbrho,ubrho,s)
    n = size(x)[1]
	xs = Vector{Base.promote_eltype(x)}(undef,n)
	
    outside = false

    for i=1:n
        xs[i] = x[i] + s[i]
    end

	nc = length(n₀)
	beta = Vector{eltype(x)}(undef,0)
	xp = Vector{Vector{eltype(x)}}(undef,0)
    rhop = Vector{eltype(x)}(undef,0)
    ix = 0
    for ip=1:np-1
    	ix += 1
    	push!(beta,xs[ix])
    	xc = Vector{eltype(x)}(undef,0)
    	for ic = 1:nc-1
    		ix += 1
    		push!(xc,xs[ix])
    	end
    	push!(xc,1.0 - sum(xc[1:nc-1]))
    	push!(xp,xc)
    	ix += 1
    	push!(rhop,xs[ix])
    end
    ix += 1
    push!(rhop,xs[ix])
    push!(beta,1.0 - sum(beta[1:np-1]))
    xc = Vector{eltype(x)}(undef,0)
    for ic = 1:nc
    	sum_moles = 0.0
    	for ip = 1:np-1
    		sum_moles += beta[ip]*xp[ip][ic]
    	end
    	push!(xc,(n₀[ic] - sum_moles)/beta[np])
    end
	push!(xp,xc)

	for ip = 1:np
		if beta[ip] < lb
	    	outside = true
		end
		if beta[ip] > ub
	    	outside = true
		end
		for ic = 1:length(xp[ip])
			if xp[ip][ic] < lb
	    		outside = true
			end
			if xp[ip][ic] > ub
	    		outside = true
			end
		end
		if rhop[ip] < lbrho
	    	outside = true
		end
		if rhop[ip] > ubrho
	    	outside = true
		end
	end

    return outside

end

function GFEMProjection_Gibbs(x,np,n₀,lb,ub,lbrho,ubrho)
	n = size(x)[1]
    p = Vector{Base.promote_eltype(x)}(undef,0)
	nc = length(n₀)
	beta = Vector{eltype(x)}(undef,0)
	xp = Vector{Vector{eltype(x)}}(undef,0)
    rhop = Vector{eltype(x)}(undef,0)

#	println("GFEMProjection_Gibbs x  = $(x)")

    ix = 0
    for ip=1:np-1
    	ix += 1
    	push!(beta,x[ix])
    	xc = Vector{eltype(x)}(undef,0)
    	for ic = 1:nc-1
    		ix += 1
    		push!(xc,x[ix])
    	end
    	push!(xc,1.0 - sum(xc[1:nc-1]))
    	push!(xp,xc)
    	ix += 1
    	push!(rhop,x[ix])
    end
    ix += 1
    push!(rhop,x[ix])
    push!(beta,1.0 - sum(beta[1:np-1]))
    xc = Vector{eltype(x)}(undef,0)
    for ic = 1:nc
    	sum_moles = 0.0
    	for ip = 1:np-1
    		sum_moles += beta[ip]*xp[ip][ic]
    	end
    	push!(xc,(n₀[ic] - sum_moles)/beta[np])
    end
	push!(xp,xc)

	outside = false
	# check bounds
	for ip = 1:np
		if beta[ip] < lb
	    	beta[ip] = lb
			outside = true
		end
		if beta[ip] > ub
	    	beta[ip] = ub
			outside = true
		end
		for ic = 1:nc
			if xp[ip][ic] < lb
	    		xp[ip][ic] = lb
				outside = true
			end
			if xp[ip][ic] > ub
	    		xp[ip][ic] = ub
				outside = true
			end
		end
		if rhop[ip] < lbrho
	    	rhop[ip] = lbrho
			outside = true
		end
		if rhop[ip] > ubrho
	    	rhop[ip] = ubrho
			outside = true
		end
	end

#	println("GFEMProjection_Gibbs outside  = $(outside)")
#	println("GFEMProjection_Gibbs beta  = $(beta)")
#	println("GFEMProjection_Gibbs xp  = $(xp)")
#	println("GFEMProjection_Gibbs rhop  = $(rhop)")

	if outside
		# normalise after bounds check
		sumbeta = sum(beta)
		beta ./= sumbeta
		for ip = 1:np
			sumxp = sum(xp[ip])
			xp[ip] ./= sumxp
		end

		# normalise overall for consisitent solution
		phasemoles = Vector{Vector{Float64}}(undef,0)
		for ic = 1:nc-1			
			summoles = 0.0
			for ip = 1:np
				summoles += xp[ip][ic] * beta[ip]
			end
			mole = Vector{Float64}(undef,0)
			for ip = 1:np
				push!(mole, xp[ip][ic] * beta[ip] * n₀[ic] / summoles)
			end
			push!(phasemoles, mole)
		end
					
		summoles = 0.0
		for ip = 1:np
			x_nc = 1.0 - sum(xp[ip][1:nc-1])
			summoles += x_nc * beta[ip]
		end
		mole = Vector{Float64}(undef,0)
		for ip = 1:np
			x_nc = 1.0 - sum(xp[ip][1:nc-1])
			push!(mole, x_nc * beta[ip] * n₀[nc] / summoles)
		end
		push!(phasemoles, mole)

		for ip = 1:np
			beta[ip] = 0.0
			for ic = 1:nc
				beta[ip] += phasemoles[ic][ip]
			end
		end

		for ip = 1:np
			for ic = 1:nc-1
				xp[ip][ic] = phasemoles[ic][ip] / beta[ip]
			end
		end

		for ip = 1:np
			xp[ip][nc] = rhop[ip]
		end

		for ip = 1:length(xp)-1
			push!(p,beta[ip])
			for ic = 1:nc
				push!(p,xp[ip][ic])
			end
		end
		push!(p,xp[np][nc])

#		println("GFEMProjection_Gibbs p  = $(p)")

		return p
	else
		return x
	end
end

function GFEM_Gibbs_func(model,p,T,n₀,v₀,np,x)
	nc = length(n₀)
	beta = Vector{eltype(x)}(undef,0)
	xp = Vector{Vector{eltype(x)}}(undef,0)
    vp = Vector{eltype(x)}(undef,0)
    ix = 0
    for ip=1:np-1
    	ix += 1
    	push!(beta,deepcopy(x[ix]))
    	xc = Vector{eltype(x)}(undef,0)
    	for ic = 1:nc-1
    		ix += 1
    		push!(xc,deepcopy(x[ix]))
    	end
    	push!(xc,1.0 - sum(xc[1:nc-1]))
    	push!(xp,xc)
    	ix += 1
    	push!(vp,deepcopy(x[ix]))
    end
    ix += 1
    push!(vp,deepcopy(x[ix]))
    push!(beta,1.0 - sum(beta[1:np-1]))
    xc = Vector{eltype(x)}(undef,0)
    for ic = 1:nc
    	sum_moles = 0.0
    	for ip = 1:np-1
    		sum_moles += beta[ip]*xp[ip][ic]
    	end
    	push!(xc,(n₀[ic] - sum_moles)/beta[np])
    end
	push!(xp,xc)
    f = 0.0
    for ip=1:np
    	v = v₀/vp[ip]
	   	f += beta[ip]*(eos(model,v,T,xp[ip]) + p*v)/R̄/T
    end
    return f
end

function GFEM_impl(model,p,T,z₀,
	max_trust_region_iters,
	tol,
	GFEM_tol,
	verbose,
	flashin)
	
	# z₀ must sum to one, i.e. it is a mole fraction vector
    nc = length(z₀)
  	# calculate reference volume based on Kays rule and vc[i] and scale to give water a vref/v ~ 1
	
	# this is expensive for SAFT eos
    pure = split_pure_model(model)
    crit = crit_pure.(pure)
    vref = 0.0
    for i= 1:nc
    	Tc,pc,vc = crit[i]
    	vref += z₀[i]*vc
    end

	beta = flashin.fractions
	np = length(beta)
	xp = flashin.compositions
	vp = flashin.volumes

	xGFEM = Vector{Float64}(undef,0)
	for ip = 1:np-1
    	push!(xGFEM,beta[ip])
    	for ic = 1:nc-1
    		push!(xGFEM,xp[ip][ic])
    	end
		push!(xGFEM,vref/vp[ip])
    end
    push!(xGFEM,vref/vp[np])

	if verbose == true
		println("GFEM Step 1 - Start Gibbs Energy Minimisation:")
	end
			
	# xsol contains all the phases xgibbs works on np-1 phases and completes the missing phase via the mole balance.
	# to be sure we must make xgibbs feasible, no negatives and no greater than 1 values.
			
	Gibbs(x) = GFEM_Gibbs_func(model,p,T,z₀,vref,np,x)
	Gibbs_g(x) = Solvers.gradient(Gibbs, x)
	Gibbs_h(x) = Solvers.hessian(Gibbs, x)

	lb = eps(Float64)*1e2
	ub = 1.0 - lb
	lbrho = 1.0e-6
	ubrho = 15.0
			
	projGibbs(x) = GFEMProjection_Gibbs(x,np,z₀,lb,ub,lbrho,ubrho)
	cnstGibbs(x,s) = GFEMConstraints_Gibbs(x,np,z₀,lb,ub,lbrho,ubrho,s)
			
	xsol,Gsol,iter,error,check = Solvers.trustregion_Dennis_Schnabel(Gibbs, Gibbs_g, Gibbs_h, projGibbs, cnstGibbs, xGFEM, max_trust_region_iters, tol, false)
			
	if verbose == true
		println("GFEM Step 2 - Gibbs Energy Minimisation: iterations taken = $(iter)")
		println("GFEM Step 2 - Gibbs Energy Minimisation: error = $(error) - tol = $(tol)")
		if check
			println("GFEM Step 2 - Gibbs Energy Minimisation: did not converge to required tolerance")
			if error < GFEM_tol && iter < max_trust_region_iters
				println("GFEM Step 2 - Gibbs Energy Minimisation: solution accepted as this is less than GFEM tolerence $(GFEM_tol)")
			else
				println("GFEM Step 2 - Gibbs Energy Minimisation: not converged")
			end
		else
			println("GFEM Step 2 - Gibbs Energy Minimisation: solution found")
		end
	end
			
	# unpack xsol using mole balance
	beta = Vector{Float64}(undef,0)
	xp = Vector{Vector{Float64}}(undef,0)
	vp = Vector{Float64}(undef,0)
	ix = 0
	for ip=1:np-1
		ix += 1
		push!(beta,xsol[ix])
		xc = Vector{Float64}(undef,0)
		for ic = 1:nc-1
			ix += 1
			push!(xc,xsol[ix])
		end
		push!(xc,1.0 - sum(xc[1:nc-1]))
		push!(xp,xc)
		ix += 1
		push!(vp,vref/xsol[ix])
	end
	ix += 1
	# volume is now returned
	push!(vp,vref/xsol[ix])
	push!(beta,1.0 - sum(beta[1:np-1]))
	xc = Vector{Float64}(undef,0)
	for ic = 1:nc
		sum_moles = 0.0
		for ip = 1:np-1
			sum_moles += beta[ip]*xp[ip][ic]
		end
		push!(xc,(z₀[ic] - sum_moles)/beta[np])
	end
	push!(xp,xc)
			
	if verbose == true
		println("GFEM Step 3 - Gibbs Energy Minimisation: Normalise final solution so it mole balances")
	end
	# normalise the solution before we finish.
	phasemoles = Vector{Vector{Float64}}(undef,0)
	for ic = 1:nc			
		summoles = 0.0
		for ip = 1:np
			summoles += xp[ip][ic] * beta[ip]
		end
		mole = Vector{Float64}(undef,0)
		for ip = 1:np
			push!(mole, xp[ip][ic] * beta[ip] * z₀[ic] / summoles)
		end
		push!(phasemoles, mole)
	end

	for ip = 1:np
		beta[ip] = 0.0
		for ic = 1:nc
			beta[ip] += phasemoles[ic][ip]
		end
	end

	for ip = 1:np
		for ic = 1:nc-1
			xp[ip][ic] = phasemoles[ic][ip] / beta[ip];
		end
		ρp = RVS_density(model,p,T,xp[ip],vref)
		vp[ip] = vref/ρp
	end
	# end of 
			
	ivpsort = sortperm(vp, rev = true)

	betas = Vector{Float64}(undef,0)
	for ip in eachindex(beta)
		push!(betas,beta[ivpsort[ip]])
	end

	xps = Vector{Vector{Float64}}(undef,0)
	for ip in eachindex(xp)
		# undo sort
		xs = Vector{Float64}(undef,length(xp[ivpsort[ip]]))
		for i in eachindex(xp[ivpsort[ip]])
			xs[i] = xp[ivpsort[ip]][i]
		end 
		push!(xps,xs)
	end

	vps = Vector{Float64}(undef,0)
	for ip in eachindex(vp)
		push!(vps,vp[ivpsort[ip]])
	end
			
	if verbose == true
		println("GFEM Step 4 - Complete")
		println("GFEM Step 5 - Phases found $(length(betas))")
		println("GFEM Step 5 - Phase moles:")
		for ip = 1:length(betas)
			println("GFEM Step 5 - Phase beta[$(ip)] = $(betas[ip])")
		end
		println("GFEM Step 5 - Phase mole fraction:")
		for ip = 1:length(betas)
			println("GFEM Step 5 - Phase x[$(ip)] = $(xps[ip])")
		end
		println("GFEM Step 5 - Phase volumes:")
		for ip = 1:length(betas)
			println("GFEM Step 5 - Phase volume[$(ip)] = $(vps[ip])")
		end
		println("GFEM Step 5 - Minimum Gibbs Energy = $(Gsol)")
	end
			
	return betas,xps,vps,Gsol

end

export GFEMTPFlash
