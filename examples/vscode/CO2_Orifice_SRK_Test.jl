using Clapeyron, Plots

function bracketmin(func::Function,a::Float64, b::Float64, xmax::Float64, verbose::Bool)
    gold = 1.618034
    ax = a
    bx = b
    fa = func(ax)
    fb = func(bx)
    if fb > fa
        tmp = ax
        ax = bx
        bx = tmp
        tmp = fa
        fa = fb
        fb = tmp
    end
    if verbose
        println("start bracket...")
        println("ax = $(ax), fa = $(fa)")
        println("bx = $(bx), fb = $(fb)")
    end
    mx = (bx + ax)/2.0
    fm = func(mx)
    if (fm < fb)
        if verbose
            println("all ready a bracket...")
            println("ax = $(ax), fa = $(fa)")
            println("mx = $(mx), fm = $(fm)")
            println("bx = $(bx), fb = $(fb)")
        end
        return ax,bx,fa,fb
    end
    dx = (bx - ax)
    cx = bx + sign(dx)*min(gold*abs(bx - ax),xmax)
    fc = func(cx)
    while (fb > fc)
        ax = bx
        fa = fb
        bx = cx
        fb = fc
        dx = (bx - ax)
        cx =  bx + sign(dx)*min(gold*abs(bx - ax),xmax)
        fc = func(cx)
        if verbose
            println("searching...")
            println("ax = $(ax), fa = $(fa)")
            println("bx = $(bx), fa = $(fb)")
            println("cx = $(cx), fc = $(fc)")
        end
    end
    if verbose
        println("found...")
        println("ax = $(ax), fa = $(fa)")
        println("bx = $(bx), fb = $(fb)")
        println("cx = $(cx), fc = $(fc)")
    end
    return ax,cx,fa,fc
end

function brentmin(
    f,
    x_lower::T,
    x_upper::T,
    verbose::Bool,
    rel_tol::T = sqrt(eps(T)),
    abs_tol::T = eps(T),
    iterations::Integer = 100,
) where {T<:AbstractFloat}

    if x_lower > x_upper
    #    error("x_lower must be less than x_upper")
        tmp = x_lower
        x_lower = x_upper
        x_upper = tmp
    end

    # Save for later
    initial_lower = x_lower
    initial_upper = x_upper

    golden_ratio::T = T(1) / 2 * (3 - sqrt(T(5.0)))

    new_minimizer = x_lower + golden_ratio * (x_upper - x_lower)
    new_minimum = f(new_minimizer)
    best_bound = "initial"
    f_calls = 1 # Number of calls to f
    step = zero(T)
    old_step = zero(T)

    old_minimizer = new_minimizer
    old_old_minimizer = new_minimizer

    old_minimum = new_minimum
    old_old_minimum = new_minimum

    iteration = 0
    converged = false

    while iteration < iterations

        p = zero(T)
        q = zero(T)

        x_tol = rel_tol * abs(new_minimizer) + abs_tol

        x_midpoint = (x_upper + x_lower) / 2

        if abs(new_minimizer - x_midpoint) <= 2 * x_tol - (x_upper - x_lower) / 2
            converged = true
            break
        end

        iteration += 1

        if abs(old_step) > x_tol
            # Compute parabola interpolation
            # new_minimizer + p/q is the optimum of the parabola
            # Also, q is guaranteed to be positive

            r = (new_minimizer - old_minimizer) * (new_minimum - old_old_minimum)
            q = (new_minimizer - old_old_minimizer) * (new_minimum - old_minimum)
            p = (new_minimizer - old_old_minimizer) * q - (new_minimizer - old_minimizer) * r
            q = 2(q - r)
            if q > 0
                p = -p
            else
                q = -q
            end
        end

        if abs(p) < abs(q * old_step / 2) &&
           p < q * (x_upper - new_minimizer) &&
           p < q * (new_minimizer - x_lower)
            old_step = step
            step = p / q

            # The function must not be evaluated too close to x_upper or x_lower
            x_temp = new_minimizer + step
            if ((x_temp - x_lower) < 2 * x_tol || (x_upper - x_temp) < 2 * x_tol)
                step = (new_minimizer < x_midpoint) ? x_tol : -x_tol
            end
        else
            old_step = (new_minimizer < x_midpoint) ? x_upper - new_minimizer : x_lower - new_minimizer
            step = golden_ratio * old_step
        end

        # The function must not be evaluated too close to new_minimizer
        if abs(step) >= x_tol
            new_x = new_minimizer + step
        else
            new_x = new_minimizer + ((step > 0) ? x_tol : -x_tol)
        end

        new_f = f(new_x)
        f_calls += 1

        if new_f < new_minimum
            if new_x < new_minimizer
                x_upper = new_minimizer
                best_bound = "upper"
            else
                x_lower = new_minimizer
                best_bound = "lower"
            end
            old_old_minimizer = old_minimizer
            old_old_minimum = old_minimum
            old_minimizer = new_minimizer
            old_minimum = new_minimum
            new_minimizer = new_x
            new_minimum = new_f
        else
            if new_x < new_minimizer
                x_lower = new_x
            else
                x_upper = new_x
            end
            if new_f <= old_minimum || old_minimizer == new_minimizer
                old_old_minimizer = old_minimizer
                old_old_minimum = old_minimum
                old_minimizer = new_x
                old_minimum = new_f
            elseif new_f <= old_old_minimum || old_old_minimizer == new_minimizer || old_old_minimizer == old_minimizer
                old_old_minimizer = new_x
                old_old_minimum = new_f
            end
        end
        if verbose
            println("xstep = $(abs(step)) xtol = $(x_tol)")
            println("new f = $(new_f)")
            println("x_lower = $(x_lower), new_minimizer = $(new_minimizer), x_upper = $(x_upper)")
        end
    end

    return new_minimizer, new_minimum

end

fluid = ["methanol","nitrogen","carbon dioxide"]
nfluid = length(fluid)
model = SRK(fluid; idealmodel=AlyLeeIdeal)

function get_props(model,fr)
    T = fr.data.T
    β = fr.fractions
    v = fr.volumes
    x = fr.compositions
    mw = Vector{Float64}(undef,0)
    h = Vector{Float64}(undef,0)
    s = Vector{Float64}(undef,0)
    for ip = eachindex(β)
        push!(mw, Clapeyron.molecular_weight(model,x[ip]))
        push!(h, Clapeyron.VT_enthalpy(model,v[ip],T,x[ip]))
        push!(s, Clapeyron.VT_entropy(model,v[ip],T,x[ip]))
    end
    mwt = 0.0
    vt  = 0.0
    ht  = 0.0
    st  = 0.0
    for ip = eachindex(β)
        mwt += β[ip]*mw[ip]
        vt  += β[ip]*v[ip]
        ht  += β[ip]*h[ip]
        st  += β[ip]*s[ip]
    end
    return β,mw,mwt,v,vt,h,ht,s,st,x
end

pbase = 15.0e5

zCO2 = 0.997275
zMeOH = 0.0001
zin = [zMeOH,(1.0 - zCO2)*(1.0 - zMeOH),zCO2*(1.0 - zMeOH)]

println("Calculate T between Tbub and Tdew to give a inital volume fraction")
println("at 40 bara its 95% and at 27 bara its 5%")
pbubble = pbase
mix_Tbub = bubble_temperature(model,pbubble,zin)

println("mix_Tbub = $(mix_Tbub)")

println("T bubble = $(mix_Tbub[1]-273.15)")

pdew = pbase
mix_Tdew = dew_temperature(model,pdew,zin)

println("mix_Tdew = $(mix_Tdew)")

println("T dew = $(mix_Tdew[1]-273.15)")

if pbase/1e5 > (40+27)/2
    T_factor = 0.937823 # 40 bara
else
    T_factor = 0.5 # 27 bara
end
Tin = T_factor*mix_Tbub[1] + (1.0 - T_factor)*mix_Tdew[1]
pin = pbase

verbose = false
flash_result = Clapeyron.tp_flash_impl(model,pin,Tin,zin,HELDTPFlash(verbose = verbose))

Tf = flash_result.data.T
println("T flash = $(Tf-273.15)")
βf = flash_result.fractions
vf = flash_result.volumes
xf = flash_result.compositions
mwf = zeros(length(βf))
for ip = eachindex(mwf)
    mwf[ip] = Clapeyron.molecular_weight(model,xf[ip])
end
vt = 0.0
for ip = eachindex(βf)
    global vt += βf[ip]*vf[ip]
end

println("Volume Total = $(vt) m3/mol")
println("Volume Fraction Vapour = $(βf[1]*vf[1]/vt)")
println("Volume Fraction Liquid = $(βf[2]*vf[2]/vt)")

xvap = xf[1]
xliq = xf[2]

#function Qs(model,p,T,z,s)
#	fr = Clapeyron.tp_flash_impl(model,p,T,z,HELDTPFlash(verbose = false))
#	g = fr.data.g
#	return -(g*T + s*T/Clapeyron.R̄)
#end
function psflashmin(model,p,T,z,s)
	fr = Clapeyron.tp_flash_impl(model,p,T,z,HELDTPFlash(verbose = false))
	g = fr.data.g
	Q =  -(g*T + s*T/Clapeyron.R̄)
	return Q
end


flash = Clapeyron.tp_flash_impl(model,pin,Tin,xvap,HELDTPFlash(verbose = false))
βf,mwf,mwtf,vf,vtf,hf,htf,sf,stf,xf = get_props(model,flash)

function flowcalc(model,p,T,z,mw,h,s,rho,Cd,d0,d1)
	Ta = T - 0.1
	Tb = T + 0.1
	psfunc(x) = psflashmin(model,p,x,z,s)
	Ta,Tc,Qa,Qc = bracketmin(psfunc,Ta, Tb, 5.0, false)
	T, Qs = brentmin(psfunc,Ta,Tc,false)
	S0 = pi/4*d0^2
    S1 = pi/4*d1^2
	flash = Clapeyron.tp_flash_impl(model,p,T,z,HELDTPFlash(verbose = false))
	βss,mwss,mwtss,vss,vtss,hss,htss,sss,stss,xss  = get_props(model,flash)
	rhoss = 1.0/vtss
	Δh = h - htss
	flow = Cd*S0*rhoss*sign(Δh/mw/(1.0 - (rhoss/rho*S0/S1)^2))*sqrt(2.0*abs(Δh/mw/(1.0 - (rhoss/rho*S0/S1)^2)))
	Q = -flow
	return Q
end

sin = stf
hin = htf
mwin = mwtf
rhoin = 1.0/vtf
Cd = 0.875
d0 = 2*25.4/1000.0/sqrt(2)
d1 = 8.0*25.4/1000.0
ps = 0.585*pin
Ts = 0.934*Tin

Ta = Ts - 0.1
Tb = Ts + 0.1
psfunc(x) = psflashmin(model,ps,x,xvap,sin)
Ta,Tc,Qa,Qc = bracketmin(psfunc,Ta, Tb, 5.0, false)
Ts, Qs = brentmin(psfunc,Ta,Tc,false)

println("Ts = $(Ts)")

#flash = Clapeyron.tp_flash_impl(model,pin,Ts,xvap,HELDTPFlash(verbose = true))
#βf,mwf,mwtf,vf,vtf,hf,htf,sf,stf,xf = get_props(model,flash)

N=100
Tg = LinRange(Ts-0.25, Ts+0.25, N)
Sin = LinRange(sin, sin, N)
Sg =  zeros(N)
Betag =  zeros(N)
Qg =  zeros(N)
Qs_plot =  zeros(N)
Ts_plot = zeros(N)

for i in 1:N
	fg = Clapeyron.tp_flash_impl(model,ps,Tg[i],xvap,HELDTPFlash(verbose = false))
	βfg,mwfg,mwtfg,vfg,vtfg,hfg,htfg,sfg,stfg,xfg = get_props(model,fg)
	Sg[i] = stfg
	if length(βfg) > 1
	    Betag[i] = βfg[1]
    else
        Betag[i] = βfg[1]
    end
	Qg[i] = psflashmin(model,ps,Tg[i],xvap,sin)
	Qs_plot[i] = Qs
	Ts_plot[i] = Ts
end

p1 = plot([Tg,Ts_plot,Tg], [Qg,Qg,Qs_plot], xlabel = "T K", ylabel = "Q [J/mol/K]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,tickfontsize=12,size=(1600,1200))
display(p1)

p2 = plot([Tg,Tg], [Sg,Sin], xlabel = "T K", ylabel = "S [J/mol/K]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,tickfontsize=12,size=(1600,1200))
display(p2)

p3 = plot([Tg], [Betag], xlabel = "T K", ylabel = "beta [-]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,tickfontsize=12,size=(1600,1200))
display(p3)


function OrificeFlow(model,pin,Tin,zin,pout,pcrit_ratio,Tcrit_ratio,Cd,d0,d1,verbose)
	flashin = Clapeyron.tp_flash_impl(model,pin,Tin,zin,HELDTPFlash(verbose = false))
	βin,mwin,mwtin,vf,vtin,hin,htin,sin,stin,xin = get_props(model,flashin)
	pcrit = pcrit_ratio*pin
	Tcrit = Tcrit_ratio*Tin
	rhotin = 1.0/vtin
	flowFunc(px) = flowcalc(model,px,Tcrit,zin,mwtin,htin,stin,rhotin,Cd,d0,d1)
	pa = pcrit - (pin - pout)/100
	pb = pcrit + (pin - pout)/100
	pa,pc,Qa,Qc = bracketmin(flowFunc,pa, pb, 1e5, false)
	pcrit, Qcrit = brentmin(flowFunc,pa,pc,false)
	Ta = Tcrit - 0.1
	Tb = Tcrit + 0.1
	psfunc(Tx) = psflashmin(model,pcrit,Tx,zin,stin)
	Ta,Tc,Qa,Qc = bracketmin(psfunc,Ta, Tb, 5.0, false)
	Tcrit, Qs = brentmin(psfunc,Ta,Tc,false)
	return -Qcrit, pcrit, Tcrit
end

pout = 0.995e5

flowvap,pcritvap,Tcritvap  = OrificeFlow(model,pin,Tin,xvap,pout,0.585,0.934,Cd,d0,d1,verbose)
println("flow vap = $(flowvap) mol/s")
println("pcrit = $(pcritvap/1e5) bara")
println("pcrit ratio = $(pcritvap/pin)")
println("Tcrit = $(Tcritvap - 273.15) deg C")
println("Tcrit ratio = $(Tcritvap/Tin)")

flash = Clapeyron.tp_flash_impl(model,pin,Tin,xliq,HELDTPFlash(verbose = false))
βf,mwf,mwtf,vf,vtf,hf,htf,sf,stf,xf = get_props(model,flash)

sin = stf
hin = htf
mwin = mwtf
rhoin = 1.0/vtf
ps = 0.760*pin
Ts = 0.965*Tin

Ta = Ts - 0.1
Tb = Ts + 0.1
psfunc(x) = psflashmin(model,ps,x,xliq,sin)
Ta,Tc,Qa,Qc = bracketmin(psfunc,Ta, Tb, 5.0, false)
Ts, Qs = brentmin(psfunc,Ta,Tc,false)

println("Ts = $(Ts)")

N=100
Tg = LinRange(Ts-0.25, Ts+0.25, N)
Sin = LinRange(sin, sin, N)
Sg =  zeros(N)
Betag =  zeros(N)
Qg =  zeros(N)
Qs_plot =  zeros(N)
Ts_plot = zeros(N)

for i in 1:N
	fg = Clapeyron.tp_flash_impl(model,ps,Tg[i],xliq,HELDTPFlash(verbose = false))
	βfg,mwfg,mwtfg,vfg,vtfg,hfg,htfg,sfg,stfg,xfg = get_props(model,fg)
	Sg[i] = stfg
    if length(βfg) > 1
	    Betag[i] = βfg[2]
    else
        Betag[i] = βfg[1]
    end
	Qg[i] = psflashmin(model,ps,Tg[i],xliq,sin)
	Qs_plot[i] = Qs
	Ts_plot[i] = Ts
end

p1 = plot([Tg,Ts_plot,Tg], [Qg,Qg,Qs_plot], xlabel = "T K", ylabel = "Q [J/mol/K]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,tickfontsize=12,size=(1600,1200))
display(p1)

p2 = plot([Tg,Tg], [Sg,Sin], xlabel = "T K", ylabel = "S [J/mol/K]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,tickfontsize=12,size=(1600,1200))
display(p2)

p3 = plot([Tg], [Betag], xlabel = "T K", ylabel = "beta [-]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,tickfontsize=12,size=(1600,1200))
display(p3)

flowliq,pcritliq,Tcritliq  = OrificeFlow(model,pin,Tin,xliq,pout,0.763,0.966,Cd,d0,d1,verbose)
println("flow liq = $(flowliq) mol/s")
println("pcrit = $(pcritliq/1e5) bara")
println("pcrit ratio = $(pcritliq/pin)")
println("Tcrit = $(Tcritliq - 273.15) deg C")
println("Tcrit ratio = $(Tcritliq/Tin)")
