using HiGHS

"""
    HELDTPFlash(;   max_HELD_iters::Int = 0
    				max_trust_region_iters::Int = 0
    				tol::Float64 = HELD_tol/10
    				HELD_tol::Float64 = sqrt(eps)
					add_pure_guess = true
    				add_anti_pure_guess = true
    				add_pure_component = [0]
    				add_random_guess = false
    				add_all_guess = false
    				verbose = true)

"""
Base.@kwdef struct HELDTPFlash <: TPFlashMethod
    max_HELD_iters::Int = 0
    max_trust_region_iters::Int = 0
    tol::Float64 = 0.01*sqrt(eps(Float64))
	# tol clean was sqrt(tol) which was too big near phase boundaries and it removed
	# phases guesses that were valid. its reduced to help with this issue
	tol_clean::Float64 = sqrt(eps(Float64))
    HELD_tol::Float64 = sqrt(eps(Float64))
	add_pure_guess::Bool = true
    add_anti_pure_guess::Bool = true
    add_pure_component::Vector{Bool} = Vector{Bool}(undef,0)
    add_random_guess::Bool = false
    add_all_guess::Bool = false
    verbose::Bool = false
end

function tp_flash_impl(model::EoSModel, p, T, n, method::HELDTPFlash)
	z₀ = n
	sumz₀ = sum(z₀)
	z₀ ./= sumz₀
	if method.max_HELD_iters == 0
		max_HELD_iters = 100*length(n)
	else
		max_HELD_iters = method.max_HELD_iters
	end
	if method.max_trust_region_iters == 0
		max_trust_region_iters = 2000*length(n)
	else
		max_trust_region_iters = method.max_trust_region_iters
	end
	tol = method.tol
	tol_clean = method.tol_clean
	HELD_tol = method.HELD_tol
	if length(method.add_pure_component) == 0 || length(method.add_pure_component) !== length(n)
		add_pure_component = fill(true,length(n))
	else
		add_pure_component = method.add_pure_component
	end
	add_pure_guess = method.add_pure_guess
	add_anti_pure_guess = method.add_anti_pure_guess
	add_random_guess = method.add_random_guess
	add_all_guess = method.add_all_guess
	verbose = method.verbose
	if verbose == true
		println("HELD  - Setup:")
		println("HELD  - trust region tolerence = $(tol)")
		println("HELD  - Cleaning tolerence = $(tol_clean)")
		println("HELD  - HELD tolerence = $(HELD_tol)")
		println("HELD  - add_pure_guess = $(add_pure_guess)")
		println("HELD  - add_anti_pure_guess = $(add_anti_pure_guess)")
		println("HELD  - add_pure_component = $(add_pure_component)")
		println("HELD  - add_random_guess = $(add_random_guess)")
		println("HELD  - add_all_guess = $(add_all_guess)")
	end
	beta,xp,vp,Gsol = HELD_impl(model,
								p,
								T,
								z₀,
								max_HELD_iters,
								max_trust_region_iters,
								tol,
								tol_clean,
								HELD_tol,
								add_pure_guess,
								add_anti_pure_guess,
								add_pure_component,
								add_random_guess,
								add_all_guess,verbose)
 
	return FlashResult(xp,beta,vp,FlashData(p,T,Gsol))
end

# new HELD

function HELDConstraints(x,lb,ub,lbrho,ubrho,s)
    n = size(x)[1]
	xp = Vector{Base.promote_eltype(x,lb,ub,s)}(undef,n)
	
    outside = false

    for i=1:n
        xp[i] = x[i] + s[i]
    end
    for i = 1:n-1
       if xp[i] > ub
	    	outside = true
	    end
		if xp[i] < lb
	    	outside = true
		end
    end
	if xp[n] > ubrho
	    outside = true
	end
	if xp[n] < lbrho
	    outside = true
	end

    return outside

end

function HELDConstraints_Gibbs(x,np,n₀,lb,ub,lbrho,ubrho,s)
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

function HELDProjection(x,lb,ub,lbrho,ubrho)
    n = size(x)[1]
    p = Vector{Base.promote_eltype(x)}(undef,n)
    xp = Vector{Base.promote_eltype(x)}(undef,n)
    for i = 1:n-1
		p[i] = x[i]
		if p[i] < lb
	    	p[i] = lb
		end
		if p[i] > ub
	    	p[i] = ub
		end
		xp[i] = p[i]
    end
    sumxp = sum(xp[1:n-1])
    xp[n] = 1.0-sumxp
    if xp[n] < lb
    	xp[n] = lb
    end
    if xp[n] > ub
    	xp[n] = ub
    end
    sumxp = sum(xp[1:n])
    xp ./= sumxp
    for i = 1:n-1
    	p[i] = xp[i]
    end
    p[n] = x[n]
    if p[n] < lbrho
		p[n] = lbrho
	end
	if p[n] > ubrho
	    p[n] = ubrho
	end
    return p
end

function HELDProjection_base(x,lb,ub)
    n = size(x)[1]
    p = Vector{Base.promote_eltype(x)}(undef,n)
    for i = 1:n
		p[i] = x[i]
		if p[i] < lb
	    	p[i] = lb
		end
		if p[i] > ub
	    	p[i] = ub
		end 
    end
    return p
end

function HELDProjection_Gibbs(x,np,n₀,lb,ub,lbrho,ubrho)
	n = size(x)[1]
    p = Vector{Base.promote_eltype(x)}(undef,0)
	nc = length(n₀)
	beta = Vector{eltype(x)}(undef,0)
	xp = Vector{Vector{eltype(x)}}(undef,0)
    rhop = Vector{eltype(x)}(undef,0)

#	println("HELDProjection_Gibbs x  = $(x)")

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

#	println("HELDProjection_Gibbs outside  = $(outside)")
#	println("HELDProjection_Gibbs beta  = $(beta)")
#	println("HELDProjection_Gibbs xp  = $(xp)")
#	println("HELDProjection_Gibbs rhop  = $(rhop)")

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

#		println("HELDProjection_Gibbs p  = $(p)")

		return p
	else
		return x
	end
end

function HELD_func(model,p,T,n₀,v₀,x,λ)
    nc = length(n₀)
    xₙ = append!(deepcopy(x[1:nc-1]),1.0 - sum(x[1:nc-1]))
    v = v₀/x[end]
    A = eos(model,v,T,xₙ)
    f = (A + p*v)/R̄/T + ∑(λ.*(n₀[1:nc-1] .- xₙ[1:nc-1]))
    return f
end

function HELD_Gibbs_func(model,p,T,n₀,v₀,np,x)
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

function HELD_initial_compositions(model,p,T,z,add_pure_guess,add_anti_pure_guess,add_pure_component,add_random_guess,add_all_guess)
	n = length(z)

	lb = eps(Float64)*1e2
	ub = 1.0 - lb
	
    xp = Vector{Vector{Float64}}(undef,0)
           
    # Wilson k-values
    K = wilson_k_values(model,p,T)
    # vapour liquid like estimates
    xvap = fill(1.,n)
    xvap = z.*K
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    xvap = HELDProjection_base(xvap,lb,ub)
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    push!(xp,xvap)
    xliq = fill(1.,n)
    xliq = z./K
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    xliq = HELDProjection_base(xliq,lb,ub)
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    push!(xp,xliq)

   	# Wilson k-values to power 1/3 for close to critical point
    Kn = K.^(1/3)
    xvap = z.*Kn
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    xvap = HELDProjection_base(xvap,lb,ub)
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    push!(xp,xvap)
    xliq = z./Kn
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    xliq = HELDProjection_base(xliq,lb,ub)
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    push!(xp,xliq)

   	# Wilson k-values to power 3 for leaner vapour and richer liquid
    Kn = K.^3
    xvap = z.*Kn
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    xvap = HELDProjection_base(xvap,lb,ub)
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    push!(xp,xvap)
    xliq = z./Kn
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    xliq = HELDProjection_base(xliq,lb,ub)
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    push!(xp,xliq)

	# Wilson k-values to power 1/2 for close to critical point
    Kn = K.^(1/2)
    xvap = z.*Kn
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    xvap = HELDProjection_base(xvap,lb,ub)
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    push!(xp,xvap)
    xliq = z./Kn
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    xliq = HELDProjection_base(xliq,lb,ub)
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    push!(xp,xliq)

   	# Wilson k-values to power 2 for leaner vapour and richer liquid
    Kn = K.^2
    xvap = z.*Kn
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    xvap = HELDProjection_base(xvap,lb,ub)
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    push!(xp,xvap)
    xliq = z./Kn
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    xliq = HELDProjection_base(xliq,lb,ub)
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    push!(xp,xliq)

	# Wilson k-values to power 2/3 for close to critical point
    Kn = K.^(2/3)
    xvap = z.*Kn
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    xvap = HELDProjection_base(xvap,lb,ub)
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    push!(xp,xvap)
    xliq = z./Kn
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    xliq = HELDProjection_base(xliq,lb,ub)
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    push!(xp,xliq)

   	# Wilson k-values to power 3/2 for leaner vapour and richer liquid
    Kn = K.^(3/2)
    xvap = z.*Kn
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    xvap = HELDProjection_base(xvap,lb,ub)
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    push!(xp,xvap)
    xliq = z./Kn
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    xliq = HELDProjection_base(xliq,lb,ub)
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    push!(xp,xliq)

	if add_pure_guess || add_all_guess
  		# pure generation
    	k = 10.0
    	for i in 1:n
    		if add_pure_component[i] == true
        		xi = fill(1.,n)
        		xi[i] = k
        		sumxi = sum(xi)
        		xi ./= sumxi
        		push!(xp,xi)
        	end
    	end
    end

	if add_pure_guess || add_all_guess
  		# pure generation
    	k = 100.0
    	for i in 1:n
    		if add_pure_component[i] == true
        		xi = fill(1.,n)
        		xi[i] = k
        		sumxi = sum(xi)
        		xi ./= sumxi
        		push!(xp,xi)
        	end
    	end
    end

    if add_pure_guess || add_all_guess
  		# pure generation
    	k = 1000.0
    	for i in 1:n
    		if add_pure_component[i] == true
        		xi = fill(1.,n)
        		xi[i] = k
        		sumxi = sum(xi)
        		xi ./= sumxi
        		push!(xp,xi)
        	end
    	end
    end

	if add_pure_guess || add_all_guess
  		# pure generation
    	k = 10000.0
    	for i in 1:n
    		if add_pure_component[i] == true
        		xi = fill(1.,n)
        		xi[i] = k
        		sumxi = sum(xi)
        		xi ./= sumxi
        		push!(xp,xi)
        	end
    	end
    end

	if add_anti_pure_guess || add_all_guess
    	# anti pure generation
    	k = 10.0
    	for i in 1:n
    		if add_pure_component[i] == true
        		xi = z
        		xi[i] = z[i]/k
        		sumxi = sum(xi)
        		xi ./= sumxi
        		xi = HELDProjection_base(xi,lb,ub)
    			sumxi = sum(xi)
    			xi ./= sumxi
        		push!(xp,xi)
        		# anti pure Wilson
        		K = wilson_k_values(model,p,T)
    			# vapour liquid like estimates
    			xvap = fill(1.,n)
    			xvap = xi.*K
    			sumxvap = sum(xvap)
    			xvap ./= sumxvap
    			xvap = HELDProjection_base(xvap,lb,ub)
    			sumxvap = sum(xvap)
    			xvap ./= sumxvap
    			push!(xp,xvap)
    			xliq = fill(1.,n)
    			xliq = xi./K
    			sumxliq = sum(xliq)
    			xliq ./= sumxliq
    			xliq = HELDProjection_base(xliq,lb,ub)
    			sumxliq = sum(xliq)
    			xliq ./= sumxliq
    			push!(xp,xliq)
    		end
    	end
    end

	if add_anti_pure_guess || add_all_guess
    	# anti pure generation
    	k = 100.0
    	for i in 1:n
    		if add_pure_component[i] == true
        		xi = z
        		xi[i] = z[i]/k
        		sumxi = sum(xi)
        		xi ./= sumxi
        		xi = HELDProjection_base(xi,lb,ub)
    			sumxi = sum(xi)
    			xi ./= sumxi
        		push!(xp,xi)
        		# anti pure Wilson
        		K = wilson_k_values(model,p,T)
    			# vapour liquid like estimates
    			xvap = fill(1.,n)
    			xvap = xi.*K
    			sumxvap = sum(xvap)
    			xvap ./= sumxvap
    			xvap = HELDProjection_base(xvap,lb,ub)
    			sumxvap = sum(xvap)
    			xvap ./= sumxvap
    			push!(xp,xvap)
    			xliq = fill(1.,n)
    			xliq = xi./K
    			sumxliq = sum(xliq)
    			xliq ./= sumxliq
    			xliq = HELDProjection_base(xliq,lb,ub)
    			sumxliq = sum(xliq)
    			xliq ./= sumxliq
    			push!(xp,xliq)
    		end
    	end
    end
    
    if add_anti_pure_guess || add_all_guess
    	# anti pure generation
    	k = 1000.0
    	for i in 1:n
    		if add_pure_component[i] == true
        		xi = z
        		xi[i] = z[i]/k
        		sumxi = sum(xi)
        		xi ./= sumxi
        		xi = HELDProjection_base(xi,lb,ub)
    			sumxi = sum(xi)
    			xi ./= sumxi
        		push!(xp,xi)
        		# anti pure Wilson
        		K = wilson_k_values(model,p,T)
    			# vapour liquid like estimates
    			xvap = fill(1.,n)
    			xvap = xi.*K
    			sumxvap = sum(xvap)
    			xvap ./= sumxvap
    			xvap = HELDProjection_base(xvap,lb,ub)
    			sumxvap = sum(xvap)
    			xvap ./= sumxvap
    			push!(xp,xvap)
    			xliq = fill(1.,n)
    			xliq = xi./K
    			sumxliq = sum(xliq)
    			xliq ./= sumxliq
    			xliq = HELDProjection_base(xliq,lb,ub)
    			sumxliq = sum(xliq)
    			xliq ./= sumxliq
    			push!(xp,xliq)
    		end
    	end
    end

	if add_anti_pure_guess || add_all_guess
    	# anti pure generation
    	k = 10000.0
    	for i in 1:n
    		if add_pure_component[i] == true
        		xi = z
        		xi[i] = z[i]/k
        		sumxi = sum(xi)
        		xi ./= sumxi
        		xi = HELDProjection_base(xi,lb,ub)
    			sumxi = sum(xi)
    			xi ./= sumxi
        		push!(xp,xi)
        		# anti pure Wilson
        		K = wilson_k_values(model,p,T)
    			# vapour liquid like estimates
    			xvap = fill(1.,n)
    			xvap = xi.*K
    			sumxvap = sum(xvap)
    			xvap ./= sumxvap
    			xvap = HELDProjection_base(xvap,lb,ub)
    			sumxvap = sum(xvap)
    			xvap ./= sumxvap
    			push!(xp,xvap)
    			xliq = fill(1.,n)
    			xliq = xi./K
    			sumxliq = sum(xliq)
    			xliq ./= sumxliq
    			xliq = HELDProjection_base(xliq,lb,ub)
    			sumxliq = sum(xliq)
    			xliq ./= sumxliq
    			push!(xp,xliq)
    		end
    	end
    end
    
    if add_random_guess || add_all_guess
    	# some random guesses
    	for ir = 1:n
        	xr = fill(0.,n)
    		for i = 1:n
    			xr[i] = 1.0e-6 + (1.0 - 2e-6)*rand()
    		end
    		sumxr = sum(xr)
    		xr ./= sumxr
    		push!(xp,xr)
   		end
    end
    
    return xp
    
end

function HELD_Pereira_compositions(model,p,T,z)
    
	n = length(z)

	lb = eps(Float64)*1e2
	ub = 1.0 - lb
	
    xp = Vector{Vector{Float64}}(undef,0)
	d = Vector{Float64}(undef,0)

	use_full_dset = false
	if use_full_dset
		nd = 2.0^(n-1)
		Dd = 1.0/nd
		dd = 0.0
		for id = 1:nd-1
			dd += Dd
			push!(d,dd)
		end
	else
		push!(d,0.5)
	end
    
	for id in eachindex(d)
		for i = 1:n-1
			x̂ = fill(0.0,n)
			x̄ = fill(0.0,n)
			for j = 1:n
				if (j == i)
					x̂[j] = d[id]*z[i]
					x̄[j] = z[i] + d[id]*(1.0 - z[i])
				else
					x̂[j] = (1.0 - d[id]*z[i])/(n-1)
					x̄[j] = (1.0 - (z[i] + d[id]*(1.0 - z[i])))/(n-1)
				end
			end
			x̂ = HELDProjection_base(x̂,lb,ub)
			sumx̂ = sum(x̂)
			x̂ ./= sumx̂
			push!(xp,x̂)
			x̄ = HELDProjection_base(x̄,lb,ub)
			sumx̄ = sum(x̄)
			x̄ ./= sumx̄
			push!(xp,x̄)
		end
	end

    return xp
    
end

function HELD_impl(model,p,T,z₀,
	max_HELD_iters,
	max_trust_region_iters,
	tol,
	tol_clean,
	HELD_tol,
	add_pure_guess,
	add_anti_pure_guess,
	add_pure_component,
	add_random_guess,
	add_all_guess,
	verbose)

	 if verbose == true
        	println("HELD Fliud In")
			println("p = $(p/1.0e5)")
			println("T = $(T -273.15)")
			println("z₀ = $(z₀)")
    end

	rr1_iter = 0
	rr1_error = 0.0
	rr2_iter = 0
	rr2_error = 0.0
	
	# z₀ must sum to one, i.e. it is a mole fraction vector
    nc = length(z₀)

	vref = RVS_getvref(model,z₀)

 	ρ₀ = RVS_density(model,p,T,z₀,vref)
	v₀ = vref/ρ₀
    μ₀ = VT_chemical_potential(model,v₀,T,z₀)
    λ₀ = (μ₀[1:nc-1] .- μ₀[nc])/R̄/T
 
    G(x) = HELD_func(model,p,T,z₀,vref,x,λ₀)
    G_g(x) = Solvers.gradient(G,x)
    G_h(x) = Solvers.hessian(G,x)

	lb = eps(Float64)*1e2
	ub = 1.0 - lb
	lbrho = 1.0e-6
	ubrho = 15.0
    projHELD(x) = HELDProjection(x,lb,ub,lbrho,ubrho)
	cnstHELD(x,s) = HELDConstraints(x,lb,ub,lbrho,ubrho,s)
 	x₀ = append!(deepcopy(z₀[1:nc-1]),ρ₀)
    G₀ = G(x₀)

    if verbose == true
    		println("HELD Step 1 - Initialisation:")
    		println("HELD Step 1 - UBDⱽ = $(G₀)")
    		println("HELD Step 1 - λ₀ = $(λ₀)")
    end
    xi = HELD_initial_compositions(model,p,T,z₀,add_pure_guess,add_anti_pure_guess,add_pure_component,add_random_guess,add_all_guess)
    
	fmins = Vector{Float64}(undef,0)
    xmins = Vector{Vector{Float64}}(undef,0)
	
	for ix = 1:length(xi)
		ρi = RVS_density(model,p,T,xi[ix],vref)
		xρi = append!(deepcopy(xi[ix][1:nc-1]),ρi)
		xmin,fmin,iter,error,check = Solvers.trustregion_Dennis_Schnabel(G, G_g, G_h, projHELD,cnstHELD, xρi, max_trust_region_iters, tol, false)

		xmin = projHELD(xmin)
		fmin = G(xmin)
		
		# we add fmin < G₀ as we are searching for instability
		if fmin < G₀ && check == false
			push!(fmins,fmin)
			push!(xmins,xmin)
		end

	end
    
    fmins_unique, xmins_unique, stable = HELD_clean_local_solutions(G₀, x₀, fmins, xmins, tol_clean, verbose)

    if verbose == true
    		println("HELD Step 1 - Phase stability check completed: $(length(fmins_unique)) unique solutions found")
			if length(fmins_unique) > 0
				println("HELD Step 1 - Phase stability, fmin unique = $(fmins_unique)")
				println("HELD Step 1 - Phase stability, xmin unique = $(xmins_unique)")
			end
    end
    if stable
        # return starting solution as it is stable
        if verbose == true
        	println("HELD Step 1 - Fluid is stable")
    	end
    	if verbose == true
			println("HELD Step 6 - Complete")
			println("HELD Step 6 - Phase found 1")
			println("HELD Step 6 - Phase moles:")
			println("HELD Step 6 - Phase beta(1) = 1")
			println("HELD Step 6 - Phase mole fraction:")
			println("HELD Step 6 - Phase x(1) = $(z₀)")
			println("HELD Step 6 - Phase volumes:")
			println("HELD Step 6 - Phase volume(1) = $(v₀)")
			println("HELD Step 6 - Minimum Gibbs Energy = $(G₀)")
    	end
    	return  [1.0], [z₀], [v₀], G₀
    else
        # return unique solutions, these need to be added to M and are good initial guesses for phases
        if verbose == true
            println("HELD Step 1 - Fluid is unstable, search for phases begins")
			println("HELD Step 1 - xmins unique = $(xmins_unique)")
            println("HELD Step 1 - Initialise set ℳ used for OPₓᵥ , and ℳguess used for local minimisations in IPₓᵥ")
    	end
    	
        # for cutting plane we are working on the Gibbs surface so λ is zero
        λi = fill(0.,nc-1)
        Gi(x) = HELD_func(model,p,T,z₀,vref,x,λi)
        
        # set up inital ℳ used for OPₓᵥ
        # ℳi is [xi[1:nc-1], Vref/Vi, Gi]
    	ℳ = Vector{Vector{Float64}}(undef,0)
		xm = HELD_Pereira_compositions(model,p,T,z₀)
		# set up initial ℳ set
		# add initial guesses and the newly found minimums from first iteration stability check
		for im = 1:length(xm)
			ρm = RVS_density(model,p,T,xm[im],vref)
		#	xρm = append!(deepcopy(xm[im][1:nc-1]),ρm)
    	#	xρGim = append!(deepcopy(xρm),Gi(xρm))
			xρm = append!(xm[im][1:nc-1],ρm)
    		xρGim = append!(xρm,Gi(xρm))
    		push!(ℳ,xρGim)
    	end
		for ii = 1:length(xi)
			ρi = RVS_density(model,p,T,xi[ii],vref)
		#	xρi = append!(deepcopy(xi[ii][1:nc-1]),ρi)
    	#	xρGii = append!(deepcopy(xρi),Gi(xρi))
			xρi = append!(xi[ii][1:nc-1],ρi)
    		xρGii = append!(xρi,Gi(xρi))
    		push!(ℳ,xρGii)
    	end
		for i = 1:length(fmins_unique)
		#	xminsGi_unique = append!(deepcopy(xmins_unique[i]),Gi(xmins_unique[i]))
			xminsGi_unique = append!(xmins_unique[i],Gi(xmins_unique[i]))
			push!(ℳ,xminsGi_unique)
		end

        # set up inital ℳguess used for local minimisations in IPₓᵥ
        # ℳguessi is [xi[1:nc-1], Vref/Vi, Gi]
		ℳguess = Vector{Vector{Float64}}(undef,0)
		for ii = 1:length(xi)
			ρi = RVS_density(model,p,T,xi[ii],vref)
		#	xρi = append!(deepcopy(xi[ii][1:nc-1]),ρi)
		#	xρGii = append!(deepcopy(xρi),Gi(xρi))
			xρi = append!(xi[ii][1:nc-1],ρi)
			xρGii = append!(xρi,Gi(xρi))
    		push!(ℳguess,xρGii)
    	end
		for i = 1:length(fmins_unique)
		#	xminsGi_unique = append!(deepcopy(xmins_unique[i]),Gi(xmins_unique[i]))
			xminsGi_unique = append!(xmins_unique[i],Gi(xmins_unique[i]))
			push!(ℳguess,xminsGi_unique)
		end
		
		n_unique_previous = length(fmins_unique)

		# now we have our first ℳ we can solve the cutting plane problem to get new λ
		# with new λ we can solve the local minimisation to get new unique minimums to add to ℳ
    	
    	UBDⱽ =  G₀
    	LBDⱽ = -Inf

		λnorm = norm(λ₀,Inf)
	#	λnorm_last = λnorm
		λmax = 1.05*λnorm
	#    λmax = Inf
	#	λᴸ = fill(-λmax,nc-1)
	#	λᵁ = fill( λmax,nc-1)

		limit_λs_by_bounds = true
    	
    	use_global_solution = false

    	λStalling_count = 0

    	np = 1
    	HELD_complete = false
    	xHELD = Vector{Float64}(undef,0)

	#	beta_best = Vector{Float64}(undef,0)
	#	fmins_best = Vector{Float64}(undef,0)
    #	xmins_best = Vector{Vector{Float64}}(undef,0)
    	
    	for k = 1:max_HELD_iters

        	if verbose == true
        		println("HELD Step 2 - iteration $(k)")
            	println("HELD Step 2 - Solve OPₓᵥ for new λˢ")
    		end
    		
			#=
			# JuMP interface
    		OPₓᵥ = Model(HiGHS.Optimizer)
    		set_optimizer_attribute(OPₓᵥ, "log_to_console", false)
    		set_optimizer_attribute(OPₓᵥ, "output_flag", false)
    
    		@variable(OPₓᵥ, v)
    		@variable(OPₓᵥ, λ[1:nc-1])
    
    		@constraint(OPₓᵥ,v <= UBDⱽ)
    		#v <= Gi + ∑(λi*(n-xi))
    		@constraint(OPₓᵥ,[i ∈ 1:length(ℳ)],v <= ℳ[i][nc+1]+sum(λ.*(z₀[1:nc-1] .- ℳ[i][1:nc-1])))

			if limit_λs_by_bounds
    			# λᴸ <= λ <= λᵁ 
    			@constraint(OPₓᵥ,[i ∈ 1:nc-1],λᴸ[i] <= λ[i] <= λᵁ[i])
			end

    		@objective(OPₓᵥ, Max, v)
			use_ipm = false
			if use_ipm
    			ipm_tol = max(tol, 1.0e-10)
    			set_attribute(OPₓᵥ, "solver", "ipm")    		
    			set_attribute(OPₓᵥ, "ipm_optimality_tolerance", ipm_tol)
			else
				simplex_tol = max(tol, 1.0e-10)
	   			set_attribute(OPₓᵥ, "solver", "simplex")
				set_attribute(OPₓᵥ, "dual_feasibility_tolerance", simplex_tol)
				set_attribute(OPₓᵥ, "primal_feasibility_tolerance", simplex_tol)
				set_attribute(OPₓᵥ, "dual_residual_tolerance", simplex_tol)
				set_attribute(OPₓᵥ, "primal_residual_tolerance", simplex_tol)
			end
    		optimize!(OPₓᵥ)
    		λˢ = JuMP.value.(λ)

    		UBDⱽ  = JuMP.value.(v)
			=#

			# C API interface to HiGHs - JuMP not required
			OPₓᵥ = Highs_create()
			Highs_setBoolOptionValue(OPₓᵥ, "log_to_console", false)
			Highs_setBoolOptionValue(OPₓᵥ, "output_flag", false)
			use_ipm = false
			if use_ipm
				Highs_setStringOptionValue(OPₓᵥ, "solver", "ipm")
				ipm_tol = max(tol, 1.0e-12)
				Highs_setDoubleOptionValue(OPₓᵥ, "ipm_optimality_tolerance", ipm_tol)
			else
				Highs_setStringOptionValue(OPₓᵥ, "solver", "simplex")
				simplex_tol = max(tol, 1.0e-10)
				Highs_setDoubleOptionValue(OPₓᵥ, "dual_feasibility_tolerance", simplex_tol)
				Highs_setDoubleOptionValue(OPₓᵥ, "primal_feasibility_tolerance", simplex_tol)
				Highs_setDoubleOptionValue(OPₓᵥ, "dual_residual_tolerance", simplex_tol)
				Highs_setDoubleOptionValue(OPₓᵥ, "primal_residual_tolerance", simplex_tol)
			end
			Highs_addCol(OPₓᵥ, 1.0, -Inf, Inf, 0, C_NULL, C_NULL)
			for iλ = 1:nc-1
				if limit_λs_by_bounds
					Highs_addCol(OPₓᵥ, 0.0, -λmax, λmax, 0, C_NULL, C_NULL)
				else
					Highs_addCol(OPₓᵥ, 0.0, -Inf, Inf, 0, C_NULL, C_NULL)
				end
			end
			Highs_changeObjectiveSense(OPₓᵥ, kHighsObjSenseMaximize)
			for iℳ = 1:length(ℳ)
				row_index = Vector{Cint}(undef,0)
				row_value = Vector{Cdouble}(undef,0)
				push!(row_index,0)
				push!(row_value, 1.0)
				for irow = 1:nc-1
					push!(row_index,irow)
					push!(row_value, ℳ[iℳ][irow] - z₀[irow] )
				end
				Highs_addRow(OPₓᵥ, -Inf, ℳ[iℳ][nc+1], nc, row_index, row_value)
			end
			Highs_run(OPₓᵥ)
			sol = zeros(Cdouble, nc)
			Highs_getSolution(OPₓᵥ, sol, C_NULL, C_NULL, C_NULL)
			UBDⱽ = sol[1]
			λˢ = sol[2:nc]
			Highs_destroy(OPₓᵥ)

			filter = 0.3
		#	nstep = Int(ceil(1.0/(0.2754*filter + 0.0008)))
			if k == 1
				if verbose == true
					println("HELD Step 3 - λ bounds added to OPₓᵥ")
				end
			end

			λnorm = norm(λˢ,Inf)
			if limit_λs_by_bounds
				λmax  = filter*1.05*λnorm + (1.0 - filter)*λmax
			end

			if verbose == true
    			println("HELD Step 2 - Update UBDⱽ and λˢ from OPₓᵥ: UBDⱽ = $(UBDⱽ)")
        		println("HELD Step 2 - λˢ = $(λˢ)")
				println("HELD Step 2 - λnorm = $(λnorm) λmax = $(λmax)")
    		end
    		
    		Dλ = λˢ - λ₀
			λStalling = false
    		if norm(Dλ,Inf) < HELD_tol
    		    if verbose == true
        			println("HELD Step 2 - λˢ has stalled use random inital guesses to provide a chance to converge")
    			end
				λStalling = true
				λStalling_count += 1
    		else
    			λStalling = false
    			λ₀ = λˢ 
    		end

    	    if verbose == true
        		println("HELD Step 3 - IPₓᵥ solve, generate cutting plane with λˢ")
    		end

			Gˢ(x) = HELD_func(model,p,T,z₀,vref,x,λˢ)
			Gˢ_g(x) = Solvers.gradient(Gˢ, x)
    		Gˢ_h(x) = Solvers.hessian(Gˢ, x)
    		
    		fmins = Vector{Float64}(undef,0)
    		xmins = Vector{Vector{Float64}}(undef,0)
    			
    		for ix = 1:length(ℳguess)

    			xmin,fmin,iter,error,check = Solvers.trustregion_Dennis_Schnabel(Gˢ, Gˢ_g, Gˢ_h, projHELD, cnstHELD, ℳguess[ix][1:nc], max_trust_region_iters, tol, false)
    			
				if fmin < UBDⱽ && check == false
					push!(fmins,fmin)
					push!(xmins,xmin)
				end

    		end

		#	if verbose == true
		#		println("HELD Step 3 - IPₓᵥ solve xmins set $(xmins)")
    	#	end
		
			fmins_unique, xmins_unique, stable = HELD_clean_local_solutions(UBDⱽ, x₀, fmins, xmins, tol_clean, verbose)
			
			if length(fmins_unique) < 1 || λStalling
				use_global_solution = true
			else
				use_global_solution = false
			end
			
			if verbose == true
        		println("HELD Step 3 - IPₓᵥ solve unique, $(length(fmins_unique)) unique solutions found")
			#	println("HELD Step 3 - IPₓᵥ solve unique, xmins_unique set $(xmins_unique)")
    		end
    		
    		ℒ = Vector{Vector{Float64}}(undef,0)
    		if length(fmins_unique) > 0
				
				# find lowest minimum of returned set.
				LBDⱽ = fmins_unique[1]
				iLBDⱽ = 1
				for i = 2:length(fmins_unique)
					if fmins_unique[i] < LBDⱽ
				    	LBDⱽ = fmins_unique[i]
				    	iLBDⱽ = i
					end 
				end
				
				# sometimes the are more than one solution that can be added
				push!(ℒ,xmins_unique[iLBDⱽ])
				for i = 1:length(fmins_unique)
					if i != iLBDⱽ && abs(fmins_unique[i] - fmins_unique[iLBDⱽ]) < tol
						push!(ℒ,xmins_unique[i])
					end
				end
				
				if verbose == true
        			println("HELD Step 3 - ℒ set contains $(length(ℒ)) items")
				#	println("HELD Step 3 - ℒ set $(ℒ)")
  				end
				
			end
			
			solution_found = true
			if use_global_solution
    			if verbose == true
        			println("HELD Step 3 - IPₓᵥ Global solution required")
  				end
				iter_rand_max = 3*nc
				for iter_rand = 1:iter_rand_max
					xr = fill(0.0,nc)
    				for i = 1:nc
    						xr[i] = 1.0e-9 + (1.0 - 2e-9)*rand()
    				end
    				sumxr = sum(xr)
    				xr ./= sumxr
    				xr = HELDProjection_base(xr,lb,ub)
    				sumxr = sum(xr)
    				xr ./= sumxr
					ρr = RVS_density(model,p,T,xr,vref)
					xρr = append!(deepcopy(xr[1:nc-1]),ρr)
					xρGr = append!(deepcopy(xρr),Gi(xρr))
    				xmin,fmin,iter,error,check = Solvers.trustregion_Dennis_Schnabel(Gˢ, Gˢ_g, Gˢ_h, projHELD, cnstHELD,  xρGr[1:nc], max_trust_region_iters, tol, false)

					xmin = projHELD(xmin)
					fmin = Gˢ(xmin)

					if fmin < UBDⱽ && check == false
    					push!(fmins,fmin)
    					push!(xmins,xmin)
    				end
				end

			#	if verbose == true
			#			println("HELD Step 3 - Global solution - xmins set $(xmins)")
  			#	end
				
				fmins_unique, xmins_unique, stable = HELD_clean_local_solutions(UBDⱽ, x₀, fmins, xmins, tol_clean, verbose)

			#	if verbose == true
			#			println("HELD Step 3 - Global solution - xmins_unique set $(xmins_unique)")
  			#	end

				ℒ = Vector{Vector{Float64}}(undef,0)
				if length(fmins_unique) > 0
					
					# find lowest minimum of returned set.
					LBDⱽ = fmins_unique[1]
					iLBDⱽ = 1
					for i = 2:length(fmins_unique)
						if fmins_unique[i] < LBDⱽ
				    		LBDⱽ = fmins_unique[i]
				    		iLBDⱽ = i
						end 
					end
					
					# sometimes the are more than one solution that can be added
					push!(ℒ,xmins_unique[iLBDⱽ])
					for i = 1:length(fmins_unique)
						if i != iLBDⱽ && abs(fmins_unique[i] - fmins_unique[iLBDⱽ]) < tol
							push!(ℒ,xmins_unique[i])
						end
					end
				
					if verbose == true
        				println("HELD Step 3 - Global solution - ℒ set contains $(length(ℒ)) items")
					#	println("HELD Step 3 - Global solution - ℒ set $(ℒ)")
  					end
  					
				else
					solution_found = false
				end	
    		end
    		
    		if !solution_found || λStalling_count > 3
    			# return starting solution as we have no phases to add to the solution so this is the best we can do
        		if verbose == true
        			println("HELD Step 1 - Global solution failed, its wise to check this solution")
    			end
				
    			if verbose == true
					println("HELD Step 6 - Complete")
					println("HELD Step 6 - Phase found 1")
					println("HELD Step 6 - Phase moles:")
					println("HELD Step 6 - Phase beta(1) = 1")
					println("HELD Step 6 - Phase mole fraction:")
					println("HELD Step 6 - Phase x(1) = $(z₀)")
					println("HELD Step 6 - Phase volumes:")
					println("HELD Step 6 - Phase volume(1) = $(v₀)")
					println("HELD Step 6 - Minimum Gibbs Energy = $(G₀)")
    			end
    			return  [1.0], [z₀], [v₀], G₀

    		end
				
			error = UBDⱽ - LBDⱽ 
			
			if verbose == true
        		println("HELD Step 3 - Update LBDⱽ from IPₓᵥ: LBDⱽ = $(LBDⱽ)")
        		println("HELD Step 3 - UBDⱽ - LBDⱽ = $(error) and tol = $(HELD_tol)")
    		end
			
			betaerror = 1.0
			
			bphase = Vector{Float64}(undef,length(xmins_unique[1]))
			beta = Vector{Float64}(undef,length(xmins_unique))
			aphase = Matrix{Float64}(undef, length(xmins_unique[1]), length(xmins_unique))
			if length(xmins_unique) > 1 
				for ib = 1:length(xmins_unique)
					sumx = 0.0
					for ia = 1:length(xmins_unique[1])-1
						aphase[ia,ib] = xmins_unique[ib][ia]
						sumx += aphase[ia,ib]
					end
					aphase[length(xmins_unique[1]),ib] = 1.0 - sumx
				end
				for ia = 1:length(xmins_unique[1])
					bphase[ia] = z₀[ia]
				end
				beta = aphase\bphase
				sumbeta = 0.0
				for ib = 1:length(xmins_unique)
					sumbeta += beta[ib]
				end
				betaerror = abs(1.0 - sumbeta)
				if verbose == true
					println("HELD Step 3 - Test phase mole balance: error = $(betaerror) and tol =  $(sqrt(HELD_tol))")
					println("HELD Step 3 - Phases found: np = $(length(beta))")
				end
			end
			
			if verbose == true
				println("HELD Step 3 - Test overall convergence:")
			end
			
			np = length(beta)
			for ip = 1:np-1
				if beta[ip] < 0.0
					beta[ip] = lb
				end
				if beta[ip] > 1.0
					beta[ip] = ub
				end
			end
			sumbeta = sum(beta[1:np-1])
			if sumbeta > 1.0
				for ip = 1:np-1
					beta[ip] /= sumbeta
					beta[ip] *= ub
				end
			end
			sumbeta = sum(beta[1:np-1])
			beta[np] = 1.0 - sumbeta

			if error < sqrt(HELD_tol) && betaerror < sqrt(HELD_tol)

				feasible = true
				rr1_iter, rr1_error = RadfordRice_BetaOnly_Solver!(z₀,beta,xmins_unique,feasible,verbose)

			#	fmins_best = deepcopy(fmins_unique)
			#	xmins_best = deepcopy(xmins_unique)
			#	beta_best  = deepcopy(beta)

			#	if verbose == true
			#		println("HELD Step 5 - fmins_best: $(fmins_best)")
			#		println("HELD Step 5 - xmins_best: $(xmins_best)")
			#		println("HELD Step 5 - beta_best: $(beta_best)")
			#	end

			#	if feasible
			#		rr2_iter, rr2_error = RadfordRice_SS_Solver!(model,p,T,z₀,vref,beta,xmins_unique,feasible,verbose)
			#		if verbose == true
			#			println("HELD Step 5 - RadfordRice_SS_Solver!: iterations taken = $(iter)")
			#			println("HELD Step 5 - RadfordRice_SS_Solver!: error = $(error) - tol = $(tol)")
			#		end
			#	end

			else
				feasible = false
			end

			if verbose == true
				if feasible
					println("HELD Step 3 - RadfordRice_BetaOnly_Solver! - feasible, $(rr1_iter), $(rr1_error)")
					println("HELD Step 3 - RadfordRice_BetaOnly_Solver! beta = $(beta)")
				else
					println("HELD Step 3 - RadfordRice_BetaOnly_Solver! - unfeasible")
					println("HELD Step 3 - RadfordRice_BetaOnly_Solver! beta = $(beta)")
				end
			end

			if (error < HELD_tol && feasible)
		#	if error < HELD_tol && betaerror < sqrt(HELD_tol)
		#	if error < HELD_tol

				if verbose == true
				    println("HELD Step 4 - Error within tolerences on UBDⱽ - LBDⱽ and phase mole balance")
        			println("HELD Step 4 - solution accepted")
    			end

	   			# normalise the solution before we do the Gibbs minimisation step.
	   			# its essential that xHELD moles balances and is a feasible solution

		#		rr2_iter, rr2_error = RadfordRice_SS_Solver!(model,p,T,z₀,vref,beta,xmins_unique,feasible,verbose)
	
	   			phasemoles = Vector{Vector{Float64}}(undef,0)
	   			for ic = 1:nc-1			
					summoles = 0.0
					for ip = 1:np
						summoles += xmins_unique[ip][ic] * beta[ip]
					end
					mole = Vector{Float64}(undef,0)
					for ip = 1:np
						push!(mole, xmins_unique[ip][ic] * beta[ip] * z₀[ic] / summoles)
					end
					push!(phasemoles, mole)
				end
				
				summoles = 0.0
				for ip = 1:np
					x_nc = 1.0 - sum(xmins_unique[ip][1:nc-1])
					summoles += x_nc * beta[ip]
				end
				mole = Vector{Float64}(undef,0)
				for ip = 1:np
				    x_nc = 1.0 - sum(xmins_unique[ip][1:nc-1])
					push!(mole, x_nc * beta[ip] * z₀[nc] / summoles)
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
						xmins_unique[ip][ic] = phasemoles[ic][ip] / beta[ip]
					end
				end
				# end of normalisation

				if verbose == true
					println("HELD Step 4 - xmins_unique: $(xmins_unique)")
				end

	   			for ip = 1:np-1
    				push!(xHELD,beta[ip])
    				for ic = 1:nc
    					push!(xHELD,xmins_unique[ip][ic])
    				end
    			end
    			push!(xHELD,xmins_unique[np][nc])

    			HELD_complete = true
				break

			end
    		
    		if verbose == true
				println("HELD Step 3 - Add new (x,V)s: ℒs to the ℳ set and all current minimums to the ℳguess set")
    		end

			use_only_lowest_min = false
		#	if error > 0.001 && !use_only_lowest_min
			if !use_only_lowest_min
				# add latest minimums to ℳguess
				for i = 1:length(fmins_unique)
					xminsGi_unique = append!(deepcopy(xmins_unique[i]),Gi(xmins_unique[i]))
					push!(ℳ,xminsGi_unique)
				end
			else
				# add lowest minimum to ℳ set
				for i = 1:length(ℒ)
					ℒGi = append!(deepcopy(ℒ[i]),Gi(ℒ[i]))
					push!(ℳ,ℒGi)
				end
			end

			# add minimums to ℳguess set
			# remove previous minimums from ℳguess only if number found is less than or equal to previous
			# this keeps the maximum number found in the set
		#	if length(fmins_unique) <= n_unique_previous
    	#		deleteat!(ℳguess, (length(ℳguess) - (n_unique_previous-1)):length(ℳguess))
    	#	end
    		# add latest minimums to ℳguess
    		for i = 1:length(fmins_unique)
				xminsGi_unique = append!(deepcopy(xmins_unique[i]),Gi(xmins_unique[i]))
				push!(ℳguess,xminsGi_unique)
			end

			n_unique_previous = length(fmins_unique)
			
			if verbose == true
        		println("HELD Step 3 - Overall convergence not satisfied return to step 2")
    		end
			   	
    	end
    	
    	if !HELD_complete
    	# we made it here without a HELD solution just return single phase and warn 	
    		if verbose == true
				println("HELD Step 5 - No solutions found")
				println("HELD Step 5 - Warning: fluid is being flagged as stable: try increasing max HELD iterations if this is a large component set problem")
				println("HELD Step 5 - Phase found 1")
				println("HELD Step 5 - Phase moles:")
				println("HELD Step 5 - Phase beta(1) = 1")
				println("HELD Step 5 - Phase mole fraction:")
				println("HELD Step 5 - Phase x(1) = $(z₀)")
				println("HELD Step 5 - Phase volumes:")
				println("HELD Step 5 - Phase volume(1) = $(v₀)")
				println("HELD Step 5 - Minimum Gibbs Energy = $(G₀)")
    		end
    		return  [1.0], [z₀], [v₀], G₀
    	else
    		# we made it here with a HELD solution use Gibbs Energy Minimisation to polish the solution this normally takes 2 to 3 iterations as we are close to the solution.
			if verbose == true
				println("HELD Step 5 - Start Gibbs Energy Minimisation:")
				println("HELD Step 5 - xHELD: $(xHELD)")
			end
			
			# xsol contains all the phases xgibbs works on np-1 phases and completes the missing phase via the mole balance.
			# to be sure we must make xgibbs feasible, no negatives and no greater than 1 values.
			
			Gibbs(x) = HELD_Gibbs_func(model,p,T,z₀,vref,np,x)
			Gibbs_g(x) = Solvers.gradient(Gibbs, x)
			Gibbs_h(x) = Solvers.hessian(Gibbs, x)

			#=
			lb = Vector{Float64}(undef,0)
		   	for ip = 1:np-1
				push!(lb,eps(Float64)*100.0)
				for ic = 1:nc-1
					push!(lb,eps(Float64)*100.0)
				end
				push!(lb,eps(Float64)*100.0)
			end
			push!(lb,eps(Float64)*100.0)
			
			ub = Vector{Float64}(undef,0)
		   	for ip = 1:np-1
				push!(ub,1.0 - eps(Float64)*100.0)
				for ic = 1:nc-1
					push!(ub,1.0 - eps(Float64)*100.0)
				end
				push!(ub,15.0)
			end
			push!(ub,15.0)
			projGibbs(x) = HELDProjection_base(x,lb,ub)
			cnstGibbs(x,s) = HELDConstraints(x,lb,ub,s)
			=#

			lb = eps(Float64)*1e2
			ub = 1.0 - lb
			lbrho = 1.0e-6
			ubrho = 15.0

			projGibbs(x) = HELDProjection_Gibbs(x,np,z₀,lb,ub,lbrho,ubrho)
			cnstGibbs(x,s) = HELDConstraints_Gibbs(x,np,z₀,lb,ub,lbrho,ubrho,s)

			UseGibbsMin = true
			if UseGibbsMin
				xsol,Gsol,iter,error,check = Solvers.trustregion_Dennis_Schnabel(Gibbs, Gibbs_g, Gibbs_h, projGibbs, cnstGibbs, xHELD, max_trust_region_iters, tol, false)
				if check
					xsol = xHELD
					check = false
				end
				xsol = projGibbs(xsol)
				Gsol = Gibbs(xsol)
				Gibbs_accepted = !check
				if verbose == true
					println("HELD Step 5 - Gibbs Energy Minimisation: iterations taken = $(iter)")
					println("HELD Step 5 - Gibbs Energy Minimisation: error = $(error) - tol = $(tol)")
					if check
						println("HELD Step 5 - Gibbs Energy Minimisation: did not converge to required tolerance")
						if error < HELD_tol && iter < max_trust_region_iters
							println("HELD Step 5 - Gibbs Energy Minimisation: solution accepted as this is less than HELD tolerence $(HELD_tol)")
						else
							println("HELD Step 5 - Gibbs Energy Minimisation: not converged")
						end
					else
						println("HELD Step 5 - Gibbs Energy Minimisation: solution found")
					end
				end
			else
				check = false
				feasible = true
				rr2_iter, rr2_error = RadfordRice_SS_Solver!(model,p,T,z₀,vref,beta,xmins_unique,feasible,verbose)
				for ip = 1:np-1
    				push!(xHELD,beta[ip])
    				for ic = 1:nc
    					push!(xHELD,xmins_unique[ip][ic])
    				end
    			end
    			push!(xHELD,xmins_unique[np][nc])
				xsol = xHELD
				xsol = projGibbs(xsol)
				Gsol = Gibbs(xsol)
				iter = rr2_iter
				error = rr2_error
				Gibbs_accepted = !check
				if verbose == true
					println("HELD Step 5 - RadfordRice_SS_Solver!: iterations taken = $(iter)")
					println("HELD Step 5 - RadfordRice_SS_Solver!: error = $(error) - tol = $(tol)")
					if check
						println("HELD Step 5 - RadfordRice_SS_Solver!: did not converge to required tolerance")
						if error < HELD_tol && iter < max_trust_region_iters
							println("HELD Step 5 - RadfordRice_SS_Solver!: solution accepted as this is less than HELD tolerence $(HELD_tol)")
						else
							println("HELD Step 5 - RadfordRice_SS_Solver!: not converged")
						end
					else
						println("HELD Step 5 - RadfordRice_SS_Solver!: solution found")
					end
				end

			end
			
			if Gibbs_accepted
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
					println("HELD Step 5 - Normalise final solution so it mole balances")
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
						xp[ip][ic] = phasemoles[ic][ip] / beta[ip]
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
					println("HELD Step 5 - Complete")
					println("HELD Step 6 - Phases found $(length(betas))")
					println("HELD Step 6 - Phase moles:")
					for ip = 1:length(betas)
						println("HELD Step 6 - Phase beta[$(ip)] = $(betas[ip])")
					end
					println("HELD Step 6 - Phase mole fraction:")
					for ip = 1:length(betas)
						println("HELD Step 6 - Phase x[$(ip)] = $(xps[ip])")
					end
					println("HELD Step 6 - Phase volumes:")
					for ip = 1:length(betas)
						println("HELD Step 6 - Phase volume[$(ip)] = $(vps[ip])")
					end
					println("HELD Step 6 - Minimum Gibbs Energy = $(Gsol)")
				end
			
				return betas,xps,vps,Gsol

			else

				# best we can do
				if verbose == true
					println("HELD Step 6 - Phase found 1")
					println("HELD Step 6 - Phase moles:")
					println("HELD Step 6 - Phase beta(1) = 1")
					println("HELD Step 6 - Phase mole fraction:")
					println("HELD Step 6 - Phase x(1) = $(z₀)")
					println("HELD Step 6 - Phase volumes:")
					println("HELD Step 6 - Phase volume(1) = $(v₀)")
					println("HELD Step 6 - Minimum Gibbs Energy = $(G₀)")
				end
				return  [1.0], [z₀], [v₀], G₀

			end
			
		end # Gibss Minimisation
		
    end # Search for phases
    
end

function HELD_clean_local_solutions(G₀, x₀, fmins, xmins, tol_clean, verbose)
    iremove = fill(false,length(fmins))
    iminfound = fill(false,length(fmins))
#	if verbose
#		println("HELD_clean_local_solutions xmins: $(xmins)")
#	end
    for ir = 1:length(fmins)    	
    	if !iremove[ir] && !iminfound[ir]
    	    imin = ir
    		fmin = fmins[ir]
    		for imins = 1:length(fmins)
    			if !iremove[imins] && !iminfound[imins]
        			if fmins[imins] < fmin
        				imin = imins
        			    fmin = fmins[imins]
    				end
    			end
    		end
			iminfound[imin] = true
    		for imins = 1:length(fmins)
    			if !iremove[imins] && !iminfound[imins]
    				distances = xmins[imins] .- xmins[imin]
    				distance = norm(distances, Inf)
    				if distance < tol_clean
    					iremove[imins] = true
    				end
    			end
    		end
    	end
    end
    # remove trival solutions
	#=
    for imins = 1:length(fmins)
    	if iminfound[imins]
    		distances = xmins[imins] .- x₀
    		distance = norm(distances, Inf)
	#		println("HELD_clean_local_solutions distance: $(distance) $(tol_clean)")
    		if distance < tol_clean
    			iminfound[imins] = false
    		end
    	end
    end
	=#
#	if verbose
#		println("HELD_clean_local_solutions iminfound: $(iminfound)")
#	end
    fmins_unique = Vector{Float64}(undef,0)
    xmins_unique = Vector{Vector{Float64}}(undef,0)
    stable = true
    for ir = 1:length(fmins)
    	if iminfound[ir]
    		if fmins[ir] < G₀
    			stable = false
				push!(fmins_unique,fmins[ir])
    			push!(xmins_unique,xmins[ir])
    		end 
    	#	push!(fmins_unique,fmins[ir])
    	#	push!(xmins_unique,xmins[ir])
    	end
    end
#	if verbose
#		println("HELD_clean_local_solutions xmins_unique: $(xmins_unique)")
#	end
    return fmins_unique, xmins_unique, stable
end

function RadfordRice_BetaOnly_Solver!(z₀,β,xin,feasible,verbose)
	
	np = length(xin)
	nc = length(z₀)
	x = Vector{Vector{Float64}}(undef,0)
	for ip = 1:np
		xp = Vector{Float64}(undef,0)
		for ic = 1:nc-1
			push!(xp,xin[ip][ic])
		end
		push!(xp,1.0 - sum(xp))
		push!(x,xp)
	end

#	if verbose == true
#		println("HELD RadfordRice_Solver β - $(β)")
#		println("HELD RadfordRice_Solver np - $(np)")
#		println("HELD RadfordRice_Solver nc - $(nc)")
#		println("HELD RadfordRice_Solver x - $(x)")
#	end

	Kv = Vector{Vector{Float64}}(undef,0)
	for ip = 1:np-1
		Kp = Vector{Float64}(undef,0)
		for ic = 1:nc
			push!(Kp,x[ip][ic]/x[np][ic])
		end
		push!(Kv,Kp)
	end
#	if verbose == true
#		println("HELD RadfordRice_Solver Kv - $(Kv)")
#	end

	β0 = zeros(np-1)
	for ip = 1:np-1
		β0[ip] = β[ip]
	end

	error = 1.0
	iter = 0
	while(error > 0.0001*sqrt(eps(Float64)) && iter < 1000)

		iter += 1
		
		function funcβ(βs,Kv)
		#	if verbose == true
		#		println("HELD RadfordRice_Solver βs - $(βs)")
		#		println("HELD RadfordRice_Solver Kv - $(Kv)")
		#	end
			t = zeros(nc)
			for ic = 1:nc
				sumβsK = 0.0
				for ip = 1:np-1
					sumβsK += βs[ip]*(1.0 - Kv[ip][ic])
				end
				t[ic] = 1.0 - sumβsK
			end
			fβ = Vector{Float64}(undef,0)
			for ip = 1:np-1
				f = 0.0
				for ic = 1:nc
					f += z₀[ic]*(1.0 - Kv[ip][ic])/t[ic]
				end
				push!(fβ,f)
			end
			jac_fβ = Matrix{Float64}(undef,np-1,np-1)
			for mp = 1:np-1
				for ip = 1:np-1
					f = 0.0
					for ic = 1:nc
						f += z₀[ic]*(1.0 - Kv[ip][ic])*(1.0 - Kv[mp][ic])/t[ic]/t[ic]
					end
					jac_fβ[ip,mp] = f
				end
			end
			return fβ,jac_fβ
		end

		fβ,jac_fβ = funcβ(β0,Kv)

	#	if verbose == true
	#		println("HELD RadfordRice_Solver β0 - $(β0)")
	#		println("HELD RadfordRice_Solver fβ - $(fβ)")
	#		println("HELD RadfordRice_Solver jac_fβ - $(jac_fβ)")
	#	end

		if length(fβ) > 1
			Dβ = jac_fβ\fβ
			error = norm(Dβ,2)
			if verbose == true
		#		println("HELD RadfordRice_Solver Dβ - $(Dβ)")
			end
			for ip = 1:np-1
				β0[ip] -= Dβ[ip]
			end
		else
			Dβ = fβ[1]/jac_fβ[1]
			if verbose == true
		#		println("HELD RadfordRice_Solver Dβ - $(Dβ)")
			end
			error = abs(Dβ[1])
			β0[1] -= Dβ[1]
		end

	#	if verbose == true
	#		println("HELD RadfordRice_Solver β0 - $(β0)")
	#	end

	end

	# println("HELD RadfordRice_BetaOnly_Solver! iter, error - $(iter), $(error)")

	lb = 100.0*eps(Float64)
	feasible = true
	for ip = 1:np-1
		β[ip] = β0[ip]
		if β[ip] < 0.0
			feasible = false
			β[ip] = lb
		end
		if β[ip] > 1.0
			feasible = false
			β[ip] = 1.0 - lb
		end
	end
	sumβ = sum(β[1:np-1])
	if sumβ > 1.0
		feasible = false
		for ip = 1:np-1
			β[ip] *= (1.0 - lb)/sumβ
		end
		sumβ = sum(β[1:np-1])
	end
	β[np] = 1.0 - sumβ

	return iter, error
	
end

function RadfordRice_SS_Solver!(model,p,T,z₀,vref,β,xin,feasible,verbose)
	
	np = length(xin)
	nc = length(z₀)
	x = Vector{Vector{Float64}}(undef,0)
	ρ = Vector{Float64}(undef,0)
	for ip = 1:np
		xp = Vector{Float64}(undef,0)
		for ic = 1:nc-1
			push!(xp,xin[ip][ic])
		end
		push!(xp,1.0 - sum(xp))
		push!(x,xp)
		push!(ρ,0.0)
	end

	if verbose == true
		println("HELD RadfordRice_Solver β - $(β)")
		println("HELD RadfordRice_Solver np - $(np)")
		println("HELD RadfordRice_Solver nc - $(nc)")
		println("HELD RadfordRice_Solver x - $(x)")
	end

	β0 = zeros(np-1)
	for ip = 1:np-1
		β0[ip] = β[ip]
	end

	Kv = Vector{Vector{Float64}}(undef,0)
	for ip = 1:np-1
		Kp = Vector{Float64}(undef,0)
		for ic = 1:nc
			push!(Kp,0.0)
		end
		push!(Kv,Kp)
	end

	Kv_last = Vector{Vector{Float64}}(undef,0)
	for ip = 1:np-1
		Kp = Vector{Float64}(undef,0)
		for ic = 1:nc
			push!(Kp,0.0)
		end
		push!(Kv_last,Kp)
	end

	Kv_error = Vector{Float64}(undef,0)
	for ip = 1:np-1
		for ic = 1:nc
			push!(Kv_error,0.0)
		end
	end

	error = 1.0
	iter = 0
	while(error > 0.01*sqrt(eps(Float64)) && iter < 250)

		iter += 1

		phi = Vector{Vector{Float64}}(undef,0)
		for ip = 1:np
			ρ_ip = RVS_density(model,p,T,x[ip],vref)
			ρ[ip] = ρ_ip
			v = vref/ρ_ip
			phi_ip = Clapeyron.VT_fugacity_coefficient(model,v,T,x[ip])
			push!(phi,phi_ip)
		end
	#	Kv = Vector{Vector{Float64}}(undef,0)
		for ip = 1:np-1
			Kp = Vector{Float64}(undef,0)
			for ic = 1:nc
		#		push!(Kp,phi[np][ic]/phi[ip][ic])
				Kv[ip][ic] = phi[np][ic]/phi[ip][ic]
			end
		#	push!(Kv,Kp)
		end
	#	if verbose == true
	#		println("HELD RadfordRice_Solver Kv - $(Kv)")
	#	end

		ie = 0
		for ip = 1:np-1
			for ic = 1:nc
				ie += 1
				Kv_error[ie] = Kv[ip][ic] - Kv_last[ip][ic]
			end
		end

		error = norm(Kv_error,2)

		for ip = 1:np-1
			for ic = 1:nc
				Kv_last[ip][ic] = Kv[ip][ic]
			end
		end
		
		function funcβ(βs,Kv)
		#	if verbose == true
		#		println("HELD RadfordRice_Solver βs - $(βs)")
		#		println("HELD RadfordRice_Solver Kv - $(Kv)")
		#	end
			t = zeros(nc)
			for ic = 1:nc
				sumβsK = 0.0
				for ip = 1:np-1
					sumβsK += βs[ip]*(1.0 - Kv[ip][ic])
				end
				t[ic] = 1.0 - sumβsK
			end
			fβ = Vector{Float64}(undef,0)
			for ip = 1:np-1
				f = 0.0
				for ic = 1:nc
					f += z₀[ic]*(1.0 - Kv[ip][ic])/t[ic]
				end
				push!(fβ,f)
			end
			jac_fβ = Matrix{Float64}(undef,np-1,np-1)
			for mp = 1:np-1
				for ip = 1:np-1
					f = 0.0
					for ic = 1:nc
						f += z₀[ic]*(1.0 - Kv[ip][ic])*(1.0 - Kv[mp][ic])/t[ic]/t[ic]
					end
					jac_fβ[ip,mp] = f
				end
			end
			return fβ,jac_fβ
		end

		fβ,jac_fβ = funcβ(β0,Kv)

	#	if verbose == true
	#		println("HELD RadfordRice_Solver β0 - $(β0)")
	#		println("HELD RadfordRice_Solver fβ - $(fβ)")
	#		println("HELD RadfordRice_Solver jac_fβ - $(jac_fβ)")
	#	end

		if length(fβ) > 1
			Dβ = jac_fβ\fβ
		#	error = norm(Dβ,2)
			if verbose == true
		#		println("HELD RadfordRice_Solver Dβ - $(Dβ)")
			end
			for ip = 1:np-1
				β0[ip] -= Dβ[ip]
			end
		else
			Dβ = fβ[1]/jac_fβ[1]
			if verbose == true
		#		println("HELD RadfordRice_Solver Dβ - $(Dβ)")
			end
		#	error = abs(Dβ[1])
			β0[1] -= Dβ[1]
		end

		if verbose == true
			println("HELD RadfordRice_Solver β0 - $(β0)")
		end

		for ip = 1:np-1
			for ic = 1:nc
				sumβsK = 0.0
				for jp = 1:np-1
					sumβsK += β0[jp]*(Kv[jp][ic] - 1.0)
				end
				x[ip][ic] = z₀[ic]*Kv[ip][ic]/(1.0 + sumβsK)
			end
		end

		for ic = 1:nc
			sumβsK = 0.0
			for jp = 1:np-1
				sumβsK += β0[jp]*(Kv[jp][ic] - 1.0)
			end
			x[np][ic] = z₀[ic]/(1.0 + sumβsK)
		end

		if verbose == true
			println("HELD RadfordRice_Solver x - $(x)")
		end

	end

	if verbose == true
 		println("HELD RadfordRice_SS_Solver! iter, error - $(iter), $(error)")
	end

	lb = 100.0*eps(Float64)
	feasible = true
	for ip = 1:np-1
		β[ip] = β0[ip]
		if β[ip] < 0.0
			feasible = false
			β[ip] = lb
		end
		if β[ip] > 1.0
			feasible = false
			β[ip] = 1.0 - lb
		end
	end
	sumβ = sum(β[1:np-1])
	if sumβ > 1.0
		feasible = false
		for ip = 1:np-1
			β[ip] *= (1.0 - lb)/sumβ
		end
		sumβ = sum(β[1:np-1])
	end
	β[np] = 1.0 - sumβ

	if feasible
		for ip = 1:np
			for ic = 1:nc-1
				xin[ip][ic]= x[ip][ic]
			end
			xin[ip][nc] = ρ[ip]
		end
	end

	return iter, error
	
end

export HELDTPFlash
