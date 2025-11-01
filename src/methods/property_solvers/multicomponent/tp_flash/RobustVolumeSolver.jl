function RVS_func_rho(model,p,T,x₀,vref,rho)
    v = vref/rho
    A = eos(model,v,T,x₀)
    f = (A + p*v)/R̄/T
    return f
end

function RVS_func_p(model,p₀,T,x₀,vref,rho)
    v = vref/rho
	p,dpdV = p∂p∂V(model,v,T,x₀)
    return p₀ - p
end

function RVS_func_dpdV(model,T,x₀,vref,rho)
    v = vref/rho
	p,dpdV = p∂p∂V(model,v,T,x₀)
    return dpdV
end

function RVS_density(model,p,T,x₀,vref)

	G(x) = RVS_func_rho(model,p,T,x₀,vref,x)
	dG(x) = RVS_func_p(model,p,T,x₀,vref,x)
	ddG(x) = RVS_func_dpdV(model,T,x₀,vref,x)
	
	# calculate rho_ideal
	rho_ideal = vref/(R̄*T/p)
	rho_min = 1.0e-6
	rho_max = 15.0

	drho = rho_ideal/2.0
	if drho < rho_min
		rho_min = drho
		drho = 2.0*rho_min
	else
		drho = (rho_min + rho_ideal)/2.0
	end

	rho1 = rho_min
	rho2 = rho1+drho
	if abs(dG(rho2)) < sqrt(eps(Float64))
		rho2 += drho
	end

	rho_bracket = Vector{Vector{Float64}}(undef,0)
	while rho2 <= rho_max
		if ddG(rho1)*ddG(rho2) < 0.0
			push!(rho_bracket,[rho1,rho2])
		end
		if rho1 >= 0.01
			drho = 0.01
		end
		rho1=rho2
		rho2=rho1+drho
		if abs(ddG(rho2)) < sqrt(eps(Float64))
			rho2 += drho
		end
	end

#	println("rho_bracket  = $(rho_bracket)")

	rho_spinodial = Vector{Float64}(undef,0)
	for ib = eachindex(rho_bracket)
		ans = Roots.find_zero(ddG, rho_bracket[ib])
		for ia = eachindex(ans)
			push!(rho_spinodial,ans[ia])
		end
	end

	rho_sp_low = rho_min
	rho_sp_high = rho_max

	if length(rho_spinodial) > 1
		rho_sp_low = rho_spinodial[1]
		rho_sp_high = rho_spinodial[end]
	end

#	println("rho_spinodial  = $(rho_spinodial)")
#	println("rho_sp_low  = $(rho_sp_low)")
#	println("rho_sp_high  = $(rho_sp_high)")

	drho = rho_ideal/2.0
	if drho < rho_min
		rho_min = drho
		drho = 2.0*rho_min
	else
		drho = (rho_min + rho_ideal)/2.0
	end

	rho1 = rho_min
	rho2 = rho1+drho
	if abs(dG(rho2)) < sqrt(eps(Float64))
		rho2 += drho
	end

	rho_bracket = Vector{Vector{Float64}}(undef,0)
	while rho2 <= rho_max
		if dG(rho1)*dG(rho2) < 0.0
			push!(rho_bracket,[rho1,rho2])
		end
		if rho1 >= 0.01
			drho = 0.01
		end
		rho1=rho2
		rho2=rho1+drho
		if abs(dG(rho2)) < sqrt(eps(Float64))
			rho2 += drho
		end
	end

#	println("rho_bracket  = $(rho_bracket)")

	rho_found = Vector{Float64}(undef,0)
	for ib = eachindex(rho_bracket)
		ans = Roots.find_zero(dG, rho_bracket[ib])
		for ia = eachindex(ans)
			stab = ddG(ans[ia])
		#	if stab > 0.0
			if stab < 0.0
				push!(rho_found,ans[ia])
			end
		end
	end

#	println("rho_found  = $(rho_found)")

	rho_stable_set = Vector{Float64}(undef,0)
	if length(rho_found) > 1
		if length(rho_spinodial) > 1
			if rho_found[1] < rho_sp_low
				push!(rho_stable_set,rho_found[1])
			end
			if rho_found[end] > rho_sp_high
				push!(rho_stable_set,rho_found[end])
			end
		else
			push!(rho_stable_set,rho_found[1])
		#	push!(rho_stable_set,rho_found[end])
		end
	else
		push!(rho_stable_set,rho_found[1])
	end

#	println("rho_stable_set  = $(rho_stable_set)")

	rho_stable = rho_stable_set[1]
	if length(rho_stable_set) > 1
		if G(rho_stable_set[end]) < G(rho_stable)
			rho_stable = rho_stable_set[end]
		end
	end

	return rho_stable

end