using Clapeyron, Plots

fluid = ["carbon dioxide","nitrogen"]
nfluid = length(fluid)
model = GERG2008(fluid)

modelCO2 = GERG2008(["carbon dioxide"])
mwCO2 = Clapeyron.molecular_weight(modelCO2,[1])

println("Calculate Triple Point for CO2 based on GERG2008:")

ptp = 5.1795e5

sat = saturation_temperature(modelCO2, ptp)
Ttp = sat[1]
vltp = sat[2]
vvtp = sat[3]

println("p tp = $(ptp/1e5) bara")
println("T tp = $(Ttp) K")
println("T tp = $(Ttp - 273.15) degC")

println("v liquid tp = $(vltp) m3/mol")
println("v vapour tp = $(vvtp) m3/mol")

hltp = Clapeyron.VT_enthalpy(modelCO2, vltp,Ttp,[1])
hvtp = Clapeyron.VT_enthalpy(modelCO2, vvtp,Ttp,[1])
println("h liquid tp = $(hltp) J/mol")
println("h liquid tp = $(hltp/mwCO2/1000) kJ/kg")
println("h vapour tp = $(hvtp) J/mol")
println("h vapour tp = $(hvtp/mwCO2/1000) kJ/kg")

println("h vaporisation tp = $(hvtp - hltp) J/mol")
println("h vaporisation tp = $((hvtp - hltp)/mwCO2/1000) kJ/kg")

sltp = Clapeyron.VT_entropy(modelCO2, vltp,Ttp,[1])
svtp = Clapeyron.VT_entropy(modelCO2, vvtp,Ttp,[1])
println("s liquid tp = $(sltp) J/mol/K")
println("s liquid tp = $(sltp/mwCO2/1000) kJ/kg/K")
println("s vapour tp = $(svtp) J/mol/K")
println("s vapour tp = $(svtp/mwCO2/1000) kJ/kg/K")

gltp = hltp - Ttp*sltp
gvtp = hvtp - Ttp*svtp
println("g liquid tp = $(gltp) J/mol")
println("g liquid tp = $(gltp/mwCO2/1000) kJ/kg")
println("g vapour tp = $(gvtp) J/mol")
println("g vapour tp = $(gvtp/mwCO2/1000) kJ/kg")

println("Triple Point for CO2 defined")
println("")

println("Calculate N2 content to give 1.5 bar bubble point above pure CO2 at 40 bara Tsat")
pbase = 27e5
psat = pbase
sat = saturation_temperature(modelCO2, psat)
Tsat = sat[1]
vl = sat[2]
vv = sat[3]

println("T sat = $(Tsat - 273.15)")

xCO2 = 0.997275
x = [xCO2,1.0-xCO2]

Tbubble = Tsat
mix_bub = bubble_pressure(model,Tbubble,x)

println("p bubble = $(mix_bub[1])")
println("p bubble - p sat = $(mix_bub[1]/1e5 - psat/1e5)")

println("Calculate N2 content complete")
println("")

println("Calculate T between Tbub and Tdew to give a inital volume fraction")
println("at 40 bara its 95% and at 27 bara its 5%")
pbubble = pbase
mix_Tbub = bubble_temperature(model,pbubble,x)

println("mix_Tbub = $(mix_Tbub)")

println("T bubble = $(mix_Tbub[1]-273.15)")

pdew = pbase
mix_Tdew = dew_temperature(model,pdew,x)

println("mix_Tdew = $(mix_Tdew)")

println("T dew = $(mix_Tdew[1]-273.15)")

if pbase/1e5 > (40+27)/2
    T_factor = 0.937823 # 40 bara
else
    T_factor = 0.0401 # 27 bara
end
Tin = T_factor*mix_Tbub[1] + (1.0 - T_factor)*mix_Tdew[1]
pin = pbase
zin = x

verbose = false
flash_result = Clapeyron.tp_flash_impl(model,pin,Tin,zin, HELDTPFlash(verbose = verbose))

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

println("Calculate complete")
println("")

function signab(a::Float64, b::Float64)
    return sign(b)*abs(a)
end

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
    mx = (bx + ax)/2.0
    fm = func(mx)
    if (fm < fb)
        if verbose
            println("all ready a braket...")
            println("ax = $(ax), fa = $(fa)")
            println("mx = $(ax), fm = $(fm)")
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
        println("bx = $(bx), fa = $(fb)")
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
    iterations::Integer = 1_000,
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
            p =
                (new_minimizer - old_old_minimizer) * q -
                (new_minimizer - old_minimizer) * r
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
            old_step =
                (new_minimizer < x_midpoint) ? x_upper - new_minimizer :
                x_lower - new_minimizer
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
            elseif new_f <= old_old_minimum ||
                   old_old_minimizer == new_minimizer ||
                   old_old_minimizer == old_minimizer
                old_old_minimizer = new_x
                old_old_minimum = new_f
            end
        end
        if verbose
            println("x_lower = $(x_lower), new_minimizer = $(new_minimizer), x_upper = $(x_upper)")
        end
    end

    return new_minimizer, new_minimum

end

# functions for fluid pressure entropy flash

function Qs(model,ps,Ts,zs,ss)
    fr = Clapeyron.tp_flash_impl(model,ps,Ts,zs, HELDTPFlash(verbose = false))
    gs = fr.data.g
    return (gs*Ts + ss*Ts/Clapeyron.R̄),fr.fractions
end

function psflashmin(T)
   Q,β = Qs(model,ps,T,zs,ss)
   return -Q,length(β)
end

function psflashmin2(T)
   Q,β = Qs(model,ps,T,zs,ss)
   return -Q
end

# CO2 Solid equation of state

println("CO2 Solid equation of state")
println("Check at calculayed triple point:")

function gibbs_dryice(p,T,Tref,href,sref,gint,sint,fach,facs)
    R̄  = 8.314462
    T0 = 150.0
    p0 = 101325.0
    hmelt = 8875.0
    g0  =  fach*((href - Tref*sref)/R̄/T0) - gint/R̄/T0
    g1  = -facs*((sref - hmelt/Tref)/R̄) + sint/R̄
    g2  = -2.0109135e+0
    g3  = -2.7976237e+0
    g4  =  2.6427834e-1
    g5  =  3.8259935e+0
    g6  =  3.1711996e-1
    g7  =  2.2087195e-3
    g8  = -1.1289668e+0
    g9  =  9.2923982e-3
    g10 =  3.3914617e+3
    gn  =  7.0
    ga0 =  3.9993365e-2
    ga1 =  2.3945101e-3
    ga2 =  3.2839467e-1
    ga3 =  5.7918471e-2
    ga4 =  2.3945101e-3
    ga5 = -2.6531689e-3
    ga6 =  1.6419734e-1
    ga7 =  1.7594802e-1
    ga8 =  2.6531689e-3
    gk0 =  2.2690751e-1
    gk1 = -7.5019750e-2
    gk2 =  2.6442913e-1
    tr = T/T0
    Dtr = tr - 1
    pr = p/p0
    Dpr = pr - 1
    g = g0 + g1*Dtr + g2*Dtr^2
    g += g3*(log((tr^2 + g4^2)/(1 + g4^2)) - 2*tr/g4*(atan(tr/g4) - atan(1/g4)))
    g += g5*(log((tr^2 + g6^2)/(1 + g6^2)) - 2*tr/g6*(atan(tr/g6) - atan(1/g6)))
    f =  ga0*(tr^2 - 1)
    f += ga1*log((tr^2 - ga2*tr + ga3)/(1 - ga2 + ga3)) + ga4*log((tr^2 + ga2*tr + ga3)/(1 + ga2 + ga3))
    f += ga5*(atan((tr - ga6)/ga7) - atan((1 - ga6)/ga7))
    f += ga8*(atan((tr + ga6)/ga7) - atan((1 + ga6)/ga7))
    K =  gk0*tr^2 + gk1*tr + gk2
    g += g7*Dpr*(exp(f) + K*g8)
    g += g9*K*((pr + g10)^((gn-1)/gn) - (1 + g10)^((gn-1)/gn))
    return R̄*T0*g
end

function entropy_dryice(p,T,Tref,href,sref,gint,sint,fach,facs)
    g(x) = gibbs_dryice(p,x,Tref,href,sref,gint,sint,fach,facs)
    s(x) = Clapeyron.Solvers.derivative(g,x)
    return -s(T)
end

function volume_dryice(p,T,Tref,href,sref,gint,sint,fach,facs)
    g(x) = gibbs_dryice(x,T,Tref,href,sref,gint,sint,fach,facs)
    v(x) = Clapeyron.Solvers.derivative(g,x)
    return v(p)
end

# we set the g0 and g1 values to match fluid triple point reference conditions

# stp0 calculated first
stp0 = entropy_dryice(ptp,Ttp,Ttp,hltp,sltp,0,0,0,0)
println("ss at tp = $(stp0) J/mol/K")

# gtp0 calculated next using stp0
gtp0 = gibbs_dryice(ptp,Ttp,Ttp,hltp,sltp,0,stp0,0,1)
println("gs at tp = $(gtp0) J/mol")

gtp = gibbs_dryice(ptp,Ttp,Ttp,hltp,sltp,gtp0,stp0,1,1)

println("gs at tp = $(gtp) J/mol")
println("gs at tp = $(gtp/mwCO2/1000) kJ/kg")

stp = entropy_dryice(ptp,Ttp,Ttp,hltp,sltp,gtp0,stp0,1,1)

println("ss at tp = $(stp) J/mol/K")
println("ss at tp = $(stp/mwCO2/1000) kJ/kg/K")

htp = gtp + Ttp*stp

println("hs at tp = $(htp) J/mol")
println("hs at tp = $(htp/mwCO2/1000) kJ/kg")

vtp = volume_dryice(ptp,Ttp,Ttp,hltp,sltp,gtp0,stp0,1,1)
println("vs at tp = $(vtp) m3/mol")
println("vs at tp = $(vtp/mwCO2) m3/kg")
println("rhos at tp = $(mwCO2/vtp) kg/m3")

println("Check complete")
println("")

# Start Blowdown

ps = 27.0e5
Ts = Tf

zCO2 = 0.997275
zs = [zCO2,1.0-zCO2]

tank_volume = 5500/9

tank_radius = 5.3/2
tank_xsarea = pi*tank_radius^2
tank_length = tank_volume/tank_xsarea

level = zeros(2)
level[1] = 0.95
level[2] = 1.0 - level[1]

mutable struct phase
    β::Float64
    volume::Float64
    moles::Float64
    mw::Float64
    v::Float64
    h::Float64
    s::Float64
    z::Vector{Float64}
end

mutable struct holdup
    p::Float64
    T::Float64
    volume::Float64         # total volume
    moles::Float64          # total moles
    mw::Float64
    v::Float64
    h::Float64
    s::Float64
    z::Vector{Float64}
    phases::Vector{phase}   # component moles
end

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

#=
println("Check ps flash at inital condition")

        ΔT = -0.33

        flash = Clapeyron.tp_flash_impl(model,ps,Ts,zs, HELDTPFlash(verbose = false))
        βf,mwf,mwtf,vf,vtf,hf,htf,sf,stf,xf = get_props(model,flash)

        ss = sf[1]
        xs = xf[1]
        ps -= 0.25e5

        mix_Tdew = dew_temperature(model,ps,xs)
        Tdew = mix_Tdew[1]

        flash = Clapeyron.tp_flash_impl(model,ps,Tdew,xs, HELDTPFlash(verbose = false))
        βf,mwf,mwtf,vf,vtf,hf,htf,sf,stf,xf = get_props(model,flash)

        println("Tdew = $(Tdew-273.15)")
        println("ss = $(ss) and sdew = $(stf)")

                    function Qs(model,p,T,x,sspec)
                        fr = Clapeyron.tp_flash_impl(model,p,T,x, HELDTPFlash(verbose = false))
                        g = fr.data.g
                        return (g*T + sspec*T/Clapeyron.R̄),fr.fractions
                    end

                    function psflashmin(T)
                        Q,β = Qs(model,ps,T,xs,ss)
                        return -Q,length(β)
                    end

                    function psflashmin2(T)
                        Q,β = Qs(model,ps,T,xs,ss)
                        return -Q
                    end

                    Ta = Ts + 2.0*ΔT
                    Qa,Npa = psflashmin(Ta)
                    Tb = Ts
                    Qb,Npb = psflashmin(Tb)
                    Tc = Ts - 1.0*ΔT
                    Qc,Npc = psflashmin(Tc)

                    println("Ta = $(Ta - 273.15) deg C, Qa = $(Qa) J/mol/K, Npa = $(Npa)")
                    println("Tb = $(Tb - 273.15) deg C, Qb = $(Qb) J/mol/K, Npb = $(Npb)")
                    println("Tc = $(Tc - 273.15) deg C, Qc = $(Qc) J/mol/K, Npc = $(Npc)")

                    bmres = brentmin(psflashmin2, Ta, Tc)

                    println("Tmin = $(bmres)")
                    println("ΔT = $(bmres[1] - Ts)")

                    N=50
                    Tm = LinRange(bmres[1]-273.15, bmres[1]-273.15,  N)
                    Tg = LinRange(Ta-273.15, Tc-273.15,  N)
                    Qg =  zeros(N)
                    βg =  zeros(N)
                    sg =  zeros(N)

                    for i in 1:N
                        Q, β = Qs(model,ps,Tg[i]+273.15,xs,ss)
                        Qg[i] = -Q
                        βg[i] =  β[1]
                        fg = Clapeyron.tp_flash_impl(model,ps,Tg[i]+273.15,xs, HELDTPFlash(verbose = false))
                        βfg,mwfg,mwtfg,vfg,vtfg,hfg,htfg,sfg,stfg,xfg = get_props(model,fg)
                        sg[i] = stfg - ss
                    end

                    p1 = plot([Tg,Tm], [Qg,Qg], xlabel = "T deg C", ylabel = "Qmod [J/mol/K]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,tickfontsize=12,size=(1600,1200))
                    display(p1)

                    p2 = plot([Tg], [βg], xlabel = "T deg C", ylabel = "βg [-]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,tickfontsize=12,size=(1600,1200))
                    display(p2)

                    p3 = plot([Tg,Tm], [sg,sg], xlabel = "T deg C", ylabel = "sg [J/mol/K]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,tickfontsize=12,size=(1600,1200))
                    display(p3)

                    flash = Clapeyron.tp_flash_impl(model,ps,bmres[1],xs, HELDTPFlash(verbose = false))
                    βf,mwf,mwtf,vf,vtf,hf,htf,sf,stf,xf = get_props(model,flash)

                    println("ss = $(ss) and sf = $(stf)")

println("Check complete")
=#

function orificeFlowOld(model,pin,Tin,zin,pout,d0,d1)

    S0 = 3.14159/4*d0^2
    S1 = 3.14159/4*d1^2

    flash = Clapeyron.tp_flash_impl(model,pin,Tin,zin, HELDTPFlash(verbose = false))
    βin,mwin,mwtin,vin,vtin,hin,htin,sin,stin,xin = get_props(model,flash)

    Cd = 0.98

    rhoin = 1.0/vtin

    Δp = 0.1e5
    pstart = 17.0e5
    pend = 16.0e5
    ns = Int(floor((pstart - pend)/Δp))
    println("ns = $(ns)")
    Δp = (pstart - pend)/ns
    ΔT = 0.6

    ps = pstart
    Ts = -25.0+273.15
    zs = zin
    ss = stin

    pstep = zeros(ns)
    flowstep =  zeros(ns)

    Ts_old = Ts
    for is = 1:ns
        ps -= Δp
        Ts -= ΔT

        function Qs(model,p,T,x,sspec)
            fr = Clapeyron.tp_flash_impl(model,p,T,x, HELDTPFlash(verbose = false))
            g = fr.data.g
            return (g*T + sspec*T/Clapeyron.R̄),fr.fractions
        end

    #    function psflashmin(T)
    #        Q,β = Qs(model,ps,T,zs,ss)
    #        return -Q,length(β)
    #    end

        function psflashmin(T)
            Q,β = Qs(model,ps,T,zs,ss)
            return -Q
        end

        ΔTs = -0.1
        Ta = Ts + ΔTs
     #   Qa,Npa = psflashmin(Ta)
        Tc = Ts - ΔTs
     #   Qc,Npc = psflashmin(Tc)

        Ta, Tc, Qa, Qc = bracketmin(psflashmin,Ta,Tc,5.0,false)

        println("Ta = $(Ta - 273.15) deg C, Qa = $(Qa) J/mol/K")
        println("Tc = $(Tc - 273.15) deg C, Qc = $(Qc) J/mol/K")

        bmres = brentmin(psflashmin, Ta, Tc, false)
        ΔT = bmres[1] - Ts_old
    #    println("ΔT = $(ΔT) deg C")
        Ts = bmres[1]
        Ts_old = Ts

        flash = Clapeyron.tp_flash_impl(model,ps,Ts,zs, HELDTPFlash(verbose = false))
        βss,mwss,mwtss,vss,vtss,hss,htss,sss,stss,xss  = get_props(model,flash)

        rhoss = 1.0/vtss

        Δh = htin - htss
        flow = Cd*S0*rhoss*sqrt(2.0*abs(Δh)/(1.0 - (rhoss/rhoin*S0/S1)^2))

        pstep[is] = ps/1e5
        flowstep[is] = flow

        println("βss = $(βss)")
        println("s in = $(ss)")
        println("s out = $(stss)")
        println("Δh = $(Δh) J/mol")
        println("Δv = $(vtss - vtin) m3/mol")
        println("(pin - pout)/rho = $(2.0*(pin - ps)/(rhoin+rhoss)) m3/mol")

        println("rhoin = $(rhoin) mol/m3")
        println("rhoss = $(rhoss) mol/m3")

        println("htin = $(htin) J/mol, htss = $(htss) J/mol, flow = $(flow*mwtin) kg/s")
        println("ΔT = $(ΔT) deg C, ps = $(ps/1e5) bara, flow = $(flow) mol/s")
        println("Ts = $(Ts - 273.15) deg C")

    end

    p1 = plot(
    label=["isentropic curve"], 
    [pstep],
    [flowstep],
    xlabel = "Pressure [bara]",
    ylabel = "Flow [mol/s]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

    display(p1)

end

function orificeFlow(model,pin,Tin,zin,pcrit,Tcrit,pout,Tout,d0,d1,Cd,verbose)

    S0 = 3.14159/4*d0^2
    S1 = 3.14159/4*d1^2

    flash = Clapeyron.tp_flash_impl(model,pin,Tin,zin, HELDTPFlash(verbose = false))
    βin,mwin,mwtin,vin,vtin,hin,htin,sin,stin,xin = get_props(model,flash)

    rhoin = 1.0/vtin
    Δp = (pin - pout)/50
 #   Ts = Tin

    function flowcalc(pf)
        function Qs(model,p,T,x,sspec)
            fr = Clapeyron.tp_flash_impl(model,p,T,x, HELDTPFlash(verbose = false))
            g = fr.data.g
            return (g*T + sspec*T/Clapeyron.R̄),fr.fractions
        end
        function psflashmin(T)
            Q,β = Qs(model,pf,T,zin,stin)
            return -Q
        end
        ΔTs = -0.1
        Ta = Ts + ΔTs
        Tc = Ts - ΔTs
        Ta,Tc,Qa,Qc = bracketmin(psflashmin,Ta,Tc,5.0,false)
        bmres = brentmin(psflashmin, Ta, Tc, false)
        Ts = bmres[1]
        flash = Clapeyron.tp_flash_impl(model,pf,Ts,zin, HELDTPFlash(verbose = false))
        βss,mwss,mwtss,vss,vtss,hss,htss,sss,stss,xss  = get_props(model,flash)
        rhoss = 1.0/vtss
        Δh = htin - htss
        flow = Cd*S0*rhoss*sqrt(2.0*abs(Δh)/(1.0 - (rhoss/rhoin*S0/S1)^2))
        Q = -flow
        return Q
    end

    flow1 = 0.0
    p2 = pcrit
    Ts = Tcrit
    Q2 = flowcalc(p2)
    T2 = Ts
#    println("T2 = $(T2 - 273.15) deg C")
    Ts = Tout
    Q3 = flowcalc(pout)
    T3 = Ts
#    println("T3 = $(T3 - 273.15) deg C")

#    println("flow1 = $(flow1) mol/s, pin = $(pin/1e5) bara")
#    println("flow2 = $(-Q2) mol/s, p2 = $(p2/1e5) bara")
#    println("flow3 = $(-Q3) mol/s, p3 = $(pout/1e5) bara")

    if Q2 < Q3
        if verbose
            println("critical flow detected")
        end
        Ts = T2
        pa = p2 + Δp
        pc = p2 - Δp
        pa,pc,Qa,Qc = bracketmin(flowcalc,pa,pc,1e5,verbose)
        bmres = brentmin(flowcalc,pa,pc,verbose)
        if verbose
            println("bracket...")
            println("flowa = $(-Qa) mol/s, pa = $(pa/1e5) bara")
            println("flowc = $(-Qc) mol/s, pc = $(pc/1e5) bara")
            println("critical flow...")
            println("flow_critical = $(-bmres[2]) mol/s, pcritical = $(bmres[1]/1e5) bara")
        end
        return -bmres[2],bmres[1],Ts,T3
    else
        if verbose
            println("flow = $(-Q3) mol/s")
        end
        return -Q3,p2,T2,T3
    end

end

pin = ps
Tin = Ts
zin = zs
pout = 1.01325e5
Tcrit = Tin
Tout = Tin
pcrit = (pin + pout)/2.0

#orificeFlowOld(model,pin,Tin,zin,pout,2*25.4/1000.0,8.0*25.4/1000.0)

flow,pcrit,Tcrit,Tout = orificeFlow(model,pin,Tin,zin,pcrit,Tcrit,pout,Tout,2*25.4/1000.0,8.0*25.4/1000.0,0.98,true)

println("orifice flow = $(flow) mol/s")
println("orifice critical pressure = $(pcrit/1e5) bara")
println("orifice critical Temperature = $(Tcrit - 273.15) deg C")
println("orifice critical ratio pf/pin = $(pcrit/pin)")
println("orifice critical pressure = $(pout/1e5) bara")
println("orifice critical Temperature = $(Tout - 273.15) deg C")


function BlowDown(model, ps, Ts, zs, tank_volume, tank_radius, tank_length, level, pe, nstep)

    mix_Tbub = bubble_temperature(model,ps,zs)
    Tbub = mix_Tbub[1]
    mix_Tdew = dew_temperature(model,ps,zs)
    Tdew = mix_Tdew[1]

    holdups = Vector{holdup}(undef,0)
    tank_xsa_holdups =Vector{Float64}(undef,0)

    if Ts <= Tbub

        y = zs
        for ic = eachindex(y)
            y[ic] = mix_Tbub[4][ic]
        end
        flashy = Clapeyron.tp_flash_impl(model,ps,Tbub,y, HELDTPFlash(verbose = false))

        βf,mwf,mwtf,vf,vtf,hf,htf,sf,stf,xf = get_props(model,flashy)

        tv = 0.0
        tm = tv/vtf
        cm = Vector{Float64}(undef,0)
        for ic = eachindex(y)
            push!(cm, y[ic]*tm)
        end
        push!(holdups, holdup(ps,Tbub,βf,mwf,vf,hf,sf,tv,tm,cm,zf,flashy))

        x = zs
        flashx = Clapeyron.tp_flash_impl(model,ps,Ts,x, HELDTPFlash(verbose = false))

        vf = flashx.volumes
        xf = flashx.compositions

        mw = Clapeyron.molecular_weight(model,x)
        v = vf[1]
        h = Clapeyron.VT_enthalpy(model,v,Ts,x)
        s = Clapeyron.VT_entropy(model,v,Ts,x)
        tv = t_v
        tm = tv/v
        cm = Vector{Float64}(undef,0)
        for ic = eachindex(x)
            push!(cm, x[ic]*tm)
        end
        push!(holdups, holdup(ps,Ts,mw,v,h,s,tv,tm,cm,x,flashx))

    elseif Ts >= Tdew

        y = zs
        flashy = Clapeyron.tp_flash_impl(model,ps,Ts,y, HELDTPFlash(verbose = false))

        βf = flashy.fractions
        vf = flashy.volumes
        xf = flashy.compositions

        println("Tdew βf = $(βf)")

        mw = Clapeyron.molecular_weight(model,y)
        v = vf[1]
        h = Clapeyron.VT_enthalpy(model,v,Ts,y)
        s = Clapeyron.VT_entropy(model,v,Ts,y)
        tv = tank_volume
        tm = tv/v
        cm = Vector{Float64}(undef,0)
        for ic = eachindex(y)
            push!(cm, y[ic]*tm)
        end
        push!(holdups, holdup(ps,Ts,mw,v,h,s,tv,tm,cm,y,flashy))

        x = zs
        for ic = eachindex(x)
            x[ic] = mix_Tdew[4][ic]
        end
        flashx = Clapeyron.tp_flash_impl(model,ps,Tdew,x, HELDTPFlash(verbose = false))

        vf = flashx.volumes
        xf = flashx.compositions

        mw = Clapeyron.molecular_weight(model,x)
        v = vf[1]
        h = Clapeyron.VT_enthalpy(model,v,Tdew,x)
        s = Clapeyron.VT_entropy(model,v,Tdew,x)
        tv = 0.0
        tm = tv/v
        cm = Vector{Float64}(undef,0)
        for ic = eachindex(x)
            push!(cm, x[ic]*tm)
        end
        push!(holdups, holdup(ps,Tdew,mw,v,h,s,tv,tm,cm,x,flashx))

    else
        flash = Clapeyron.tp_flash_impl(model,ps,Ts,zs, HELDTPFlash(verbose = false))
        βf,mwf,mwtf,vf,vtf,hf,htf,sf,stf,xf = get_props(model,flash) 
        for ip = eachindex(βf)
            phases = Vector{phase}(undef,0)
            push!(phases, phase(1.0,level[ip]*tank_volume,level[ip]*tank_volume/vf[ip],mwf[ip],vf[ip],hf[ip],sf[ip],xf[ip]))
            push!(holdups, holdup(ps,Ts,level[ip]*tank_volume,level[ip]*tank_volume/vf[ip],mwf[ip],vf[ip],hf[ip],sf[ip],xf[ip],phases))
            push!(tank_xsa_holdups, level[ip]*tank_volume/tank_length)
        end
    end

    # level
    theta = 80/180*pi
    error = 1.0
    while (error > 0.0001)
        f =  (tank_radius^2)*(theta - sin(theta)) - 2.0*tank_xsa_holdups[2]
        df_dtheta = (tank_radius^2)*(1.0 - cos(theta))
        theta = theta - f/df_dtheta
        error = abs(f)
        println("theta = $(theta), error = $(error)")
    end

    level[2] = 1.0 - cos(theta/2.0)
    level[1] = 1.0 - level[2]
    println("level[2] = $(level[2])")

    # set RO critical and outlet temperatures

    Tcrit = holdups[1].T
    Tout = holdups[1].T
    pout = 1.01325e5
    pcrit = (ps + pout)/2.0

    println("")
    println("Blowndown Function Start:")
    println("Initialization:")
    println("")
    println("holdups = $(holdups)")
    
    Δp = (ps - pe)/nstep
    ΔT = fill(-0.33,length(holdups))
    ΔQ = zeros(length(holdups))

    time_Plot = Vector{Float64}(undef,0)
    T1_Plot = Vector{Float64}(undef,0)
    p1_Plot = Vector{Float64}(undef,0)
    m1_Plot = Vector{Float64}(undef,0)

    T2_Plot = Vector{Float64}(undef,0)
    p2_Plot = Vector{Float64}(undef,0)
    m2_Plot = Vector{Float64}(undef,0)

    push!(T1_Plot, holdups[1].T-273.15)
    push!(T2_Plot, holdups[2].T-273.15)

    push!(p1_Plot, holdups[1].p/1e5)
    push!(p2_Plot, holdups[2].p/1e5)

    push!(m1_Plot, holdups[1].moles*holdups[1].mw/1000)
    push!(m2_Plot, holdups[2].moles*holdups[2].mw/1000)

    time = 0.0
    push!(time_Plot, time)

    for istep=1:nstep

        println("")
        println("Step: $(istep)")
        
        for ih = eachindex(holdups)

            println("")
            ps = holdups[ih].p - Δp
            println("holdup[$(ih)]")

            if holdups[ih].moles > 0.0

                ΔsLast = 0.0
                ΔsError = 1
                println("entropy iteration")
                is = 0
                Ts = holdups[ih].T + ΔT[ih]
                while(ΔsError > 0.001)
                    is += 1
                #    Ts = holdups[ih].T + ΔT[ih]
                    zs = holdups[ih].z
                    Δs = ΔQ[ih]*log(Ts/holdups[ih].T)/(Ts - holdups[ih].T)
                    ΔsError = abs(Δs - ΔsLast)
                    ΔsLast = Δs
                    ss = holdups[ih].s + Δs

                    println("entropy step $(is)")
                    println("ps = $(ps)")
                    println("Ts = $(Ts)")
                    println("ΔsError = $(ΔsError)")
                    println("ss = $(ss)")

                    function Qs(model,p,T,x,sspec)
                        fr = Clapeyron.tp_flash_impl(model,p,T,x, HELDTPFlash(verbose = false))
                        g = fr.data.g
                        return (g*T + sspec*T/Clapeyron.R̄),fr.fractions
                    end

                #    function psflashmin(T)
                #        Q,β = Qs(model,ps,T,zs,ss)
                #        return -Q,length(β)
                #    end

                    function psflashmin(T)
                        Q,β = Qs(model,ps,T,zs,ss)
                        return -Q
                    end

                    Ta = Ts + 2.0*ΔT[ih]
                    Tc = Ts - 2.0*ΔT[ih]

                    Ta, Tc, Qa, Qc = bracketmin(psflashmin,Ta,Tc,5.0,false)

                    bmres = brentmin(psflashmin, Ta, Tc, false)
                    ΔT[ih] = bmres[1] - holdups[ih].T
                    Ts = bmres[1]

                end

                    holdups[ih].p = ps
                    holdups[ih].T = Ts
                #    holdups[ih].z = zs
                    flash = Clapeyron.tp_flash_impl(model,holdups[ih].p,holdups[ih].T,holdups[ih].z, HELDTPFlash(verbose = false))
                #    println("flash[$(ih)] = $(flash)")
                    βs,mws,mwts,vs,vts,hs,hts,ss,sts,xs  = get_props(model,flash)
                    holdups[ih].mw = mwts
                    holdups[ih].v  = vts
                    holdups[ih].h  = hts
                    holdups[ih].s  = sts
                    holdups[ih].volume  = holdups[ih].moles*holdups[ih].v
                    phases = Vector{phase}(undef,0)
                    for ip = eachindex(βs)      
                        push!(phases, phase(βs[ip],βs[ip]*holdups[ih].moles*vs[ip],βs[ip]*holdups[ih].moles,mws[ip],vs[ip],hs[ip],ss[ip],xs[ip]))
                    end
                    holdups[ih].phases = phases

            end

            # mix produced vapour and liquid from each phase this is isenthalpic process
            # ih = 1 is always the lightest phase we assume lightest phases move up to previous phase and heaviest moves down
            # we need some bookkeeping here if we assume a two phase system 
            # then 
            # tm[1] = β[1][1]*holdups[1].tm + β[1][2]*holdups[2].tm
            # tm[2] = β[2][2]*holdups[2].tm + β[2][1]*holdups[1].tm
            #
            # three phase
            # tm[1] = β[1][1]*holdups[1].tm + β[1][2]*holdups[2].tm
            # tm[2] = β[2][2]*holdups[2].tm + β[2][1]*holdups[1].tm + β[3][1]*holdups[1].tm + β[1][3]*holdups[3].tm + β[2][3]*holdups[3].tm
            # tm[3] = β[3][3]*holdups[2].tm + β[3][2]*holdups[1].tm
            #
            # a phase can only exchange with its neighbours so in 4 phase system consider phase 3
            # tm[3] = β[3][3]*holdups[3].tm + β[3][2]*holdups[1].tm + β[4][2]*holdups[1].tm + β[1][4]*holdups[4].tm + β[2][4]*holdups[4].tm + β[3][4]*holdups[4].tm
            #
            # in general for a phase ih its everything from phase ih-1 β ih to nh and phase ih+1 evrything from 1 to ih
            #
            # if ih = 1 then we have no ih-1 contribution and if ih = nh we have no ih+1 contribution
            #
            #its possible that during a step a phase can disappear and conversly a phase can appear. So the β[ip = 1 to np][ih] we need to choose the phase with the 
            # cloest density to the density in the ih phase to select ip of the β[ip = 1 to np] that the separation is occur about. ip < ip selected go up  and ip > ip selected go down
            #
            # if the phase number of the nh hold up is great than nh we need a new holdup
            # ih = 1 is always the lightest phase
            #
            # its possible for a phase to disappear we set tm = 0 and it is kept in the holdup list in case it is reactivated by moles being added. 
            # if tm = 0 for a phase it contriutes nothing to the surrounding phases and does not need flashing. The T for the phase is assumed to be the average of the surrounding phases
            #
            # most likely way is we are liquid draining and at some point all the liquid is removed. However, even in top or vapour removal a light liquid may evaporate and level a 
            # heavier phase below.
            #
            # the phase that disappers or becomes single phase and more like the phase above or below it. We could check the density and mix in fully with the phase above or below and 
            # then remove it
            # 
            # for a phase to disappear then the flashes for all holupds must have produced less phases than current holdups. Likewise if a phase reappers then the number of phases or 
            # one of the holdups must be greater that the total holdups
            
            println("")
            println("Results:")
            println("ps[$(ih)] = $(holdups[ih].p/1e5) bara")
            println("ΔT[$(ih)] = $(ΔT[ih])")
            println("Ts[$(ih)] = $(holdups[ih].T-273.15) deg C")
            println("holdups[$(ih)] = $(holdups[ih])")

        end

        # recombine vapour with liquid phase vapour and liquid with vapour phase liquid
        cm = Vector{Float64}(undef,0)
        mt1 = holdups[1].moles
        β1 = holdups[1].phases[1].β
        z1 = holdups[1].phases[1].z
        h1 = holdups[1].phases[1].h
        h1m = 0.0
        for ic = eachindex(holdups[1].z)
            if length(holdups[2].phases) > 1
                mt2 = holdups[2].moles
                β2 = holdups[2].phases[1].β
                z2 = holdups[2].phases[1].z
                h2 = holdups[2].phases[1].h
                push!(cm,β1*mt1*z1[ic] + β2*mt2*z2[ic])
                h1m = β1*mt1*h1 + β2*mt2*h2
            else
                push!(cm,β1*mt1*z1[ic])
                h1m = β1*mt1*h1
            end
        end
        mt1 = sum(cm)
        z1 = Vector{Float64}(undef,0)
        for ic = eachindex(holdups[1].z)
            push!(z1,cm[ic]/mt1)
        end
        h1 = h1m/mt1
        holdups[1].moles = mt1
        holdups[1].h = h1
        holdups[1].z = z1
    #    println("mt1 = $(mt1) moles")
    #    println("z1 = $(z1)")
    #    println("h1 = $(h1) J/mol")
        cm = Vector{Float64}(undef,0)
        h2 = 0.0
        if length(holdups[2].phases) > 1
            mt2 = holdups[2].moles
            β2 = holdups[2].phases[2].β
            z2 = holdups[2].phases[2].z
            h2 = holdups[2].phases[2].h
        else
            mt2 = holdups[2].moles
            β2 = holdups[2].phases[1].β
            z2 = holdups[2].phases[1].z
            h2 = holdups[2].phases[1].h
        end
        h2m = 0.0
        for ic = eachindex(holdups[2].z)
            if length(holdups[1].phases) > 1
                mt1 = holdups[1].moles
                β1 = holdups[1].phases[2].β
                z1 = holdups[1].phases[2].z
                h1 = holdups[1].phases[2].h
                push!(cm,β2*mt2*z2[ic] + β1*mt1*z1[ic])
                h2m = β2*mt2*h2 + β1*mt1*h1
            else
                push!(cm,β2*mt2*z2[ic])
                h2m = β2*mt2*h2
            end
        end
        mt2 = sum(cm)
        z2 = Vector{Float64}(undef,0)
        for ic = eachindex(holdups[2].z)
            push!(z2,cm[ic]/mt2)
        end
        h2 = h2m/mt2
        holdups[2].moles = mt2
        holdups[2].h = h2
        holdups[2].z = z2
    #    println("mt2 = $(mt2) moles")
    #    println("z2 = $(z2)")
    #    println("h2 = $(h2) J/mol")

        # Ok we have all we need to do a phflash and reset the holpups with the new conditions after separation and mixing

                for ih = eachindex(holdups)

                    hs = holdups[ih].h
                    Ts = holdups[ih].T
                    ps = holdups[ih].p
                    zs = holdups[ih].z

                    function Qh(model,p,T,x,hspec)
                        fr = Clapeyron.tp_flash_impl(model,p,T,x, HELDTPFlash(verbose = false))
                        g = fr.data.g
                        return (g - hspec/Clapeyron.R̄/T),fr.fractions
                    end

                #    function phflashmin(T)
                #        Q,β = Qh(model,ps,T,zs,hs)
                #        return -Q,length(β)
                #    end

                    function phflashmin(T)
                        Q,β = Qh(model,ps,T,zs,hs)
                        return -Q
                    end

                    Ta = Ts + 2.0*ΔT[ih]
                    Tc = Ts - 2.0*ΔT[ih]

                    Ta, Tc, Qa, Qc = bracketmin(phflashmin,Ta,Tc,5.0,false)

                    bmres = brentmin(phflashmin, Ta, Tc, false)

                    holdups[ih].T = bmres[1]
                #    println("Ts[$(ih)] = $(bmres[1]) K")

                    flash = Clapeyron.tp_flash_impl(model,holdups[ih].p,holdups[ih].T,holdups[ih].z, HELDTPFlash(verbose = false))
                    βs,mws,mwts,vs,vts,hs,hts,ss,sts,xs  = get_props(model,flash)
                    holdups[ih].mw = mwts
                    holdups[ih].v  = vts
                    holdups[ih].h  = hts
                    holdups[ih].s  = sts
                    holdups[ih].volume  = holdups[ih].moles*holdups[ih].v
                    phases = Vector{phase}(undef,0)
                    for ip = eachindex(βs)      
                        push!(phases, phase(βs[ip],βs[ip]*holdups[ih].moles*vs[ip],βs[ip]*holdups[ih].moles,mws[ip],vs[ip],hs[ip],ss[ip],xs[ip]))
                    end
                    holdups[ih].phases = phases

                end
        #

        # assume vapour removal
        tank_volume_new = (holdups[1].volume + holdups[2].volume)
        Δv = tank_volume_new - tank_volume
        println("Δv = $(Δv)")
        holdups[1].volume -= Δv
        Δmoles = holdups[1].volume/holdups[1].v - holdups[1].moles
        Δmoles = holdups[1].moles - holdups[1].volume/holdups[1].v
        holdups[1].moles = holdups[1].volume/holdups[1].v

    #    pout = 1.01325e5
        Cd = 0.98
        d1 = 2.0*25.4/1000.0
        d2 = 8.0*25.4/1000.0
        flowout, pcrit, Tcrit, Tout = orificeFlow(model,holdups[1].p,holdups[1].T,holdups[1].z,pcrit,Tcrit,pout,Tout,d1,d2,Cd,false)

    #    println("orifice flow = $(flowout) mol/s")
    #    println("orifice critical pressure = $(pcrit/1e5) bara")
    #    println("orifice isentropic critical Temperature = $(Tcrit - 273.15) deg C")
    #    println("orifice critical ratio pf/pin = $(pcrit/pin)")
    #    println("orifice outlet pressure = $(pout/1e5) bara")
    #    println("orifice isentropic outlet Temperature = $(Tout - 273.15) deg C")

        # to find actual valve outlet condition we need an isenthalpy flash from pin to pout.
        
        Δt = Δmoles/flowout
        time += Δt
        println("time = $(time/3600)")
        push!(time_Plot, time/3600.0)

        # new level
        for ih = eachindex(holdups)
            tank_xsa_holdups[ih] = holdups[ih].volume/tank_length
        end

        error = 1.0
        while (error > 0.0001)
            f =  (tank_radius^2)*(theta - sin(theta)) - 2.0*tank_xsa_holdups[2]
            df_dtheta = (tank_radius^2)*(1.0 - cos(theta))
            theta = theta - f/df_dtheta
            error = abs(f)
            println("theta = $(theta), error = $(error)")
        end

        level[2] = 1.0 - cos(theta/2.0)
        level[1] = 1.0 - level[2]
        println("level[2] = $(level[2])")

        cord = 2.0*tank_radius*sin(theta/2.0)
        interface_area = cord*tank_length
        println("interface_area = $(interface_area) m2")

        # used for position along wall as an arc length used to define position of interface
        arc_length = pi*tank_radius
        arc_length_liquid = theta/2.0*tank_radius
        println("arc_length_liquid/arc_length = $(arc_length_liquid/arc_length)")

        wall_vapour_area = (arc_length - arc_length_liquid)*tank_length
        println("wall_vapour_area = $(wall_vapour_area) m2")

        wall_liquid_area = arc_length_liquid*tank_length
        println("wall_liquid_area = $(wall_liquid_area) m2")

        # need htc for vapour wall and liquid wall based on transport properties, Re, Pr, Gr etc
        htc_vap = 1000.0
        htc_liq_convective = 5000.0
        htc_boil = 20000.0
        htc_liq = sqrt(htc_boil^2 + htc_liq_convective^2)
        
        # then we use Ferrite to model finite element wall and do a transient step for the Δt

        # htc for interface ~ 1/alpha = 1/alpha_vap + 1/alpha_liq, alpha_liq is not boiling part is just convective
        htc_vap_liq = 1.0/(1.0/htc_vap + 1.0/htc_liq_convective)

        if holdups[1].moles > 0.0 && holdups[2].moles > 0.0
            ΔQ[1] = Δt*htc_vap_liq*interface_area*(holdups[2].T - holdups[1].T)/holdups[1].moles
            ΔQ[2] = Δt*htc_vap_liq*interface_area*(holdups[1].T - holdups[2].T)/holdups[2].moles
        else
            ΔQ[1] = 0.0
            ΔQ[2] = 0.0
        end

        push!(T1_Plot, holdups[1].T-273.15)
        push!(T2_Plot, holdups[2].T-273.15)

        push!(p1_Plot, holdups[1].p/1e5)
        push!(p2_Plot, holdups[2].p/1e5)

        push!(m1_Plot, holdups[1].moles*holdups[1].mw/1000)
        push!(m2_Plot, holdups[2].moles*holdups[2].mw/1000)

    end

    p1 = plot(
    label=["vapour" "liquid"], 
    [T1_Plot,T2_Plot],
    [p1_Plot,p2_Plot],
    xlabel = "Temperaure [deg C]",
    ylabel = "Pressure [bara]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

    display(p1)

    p2 = plot(
    label=["vapour" "liquid"], 
    [m1_Plot,m2_Plot],
    [p1_Plot,p2_Plot],
    xlabel = "Mass [tons]",
    ylabel = "Pressure [bara]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

    display(p2)

    p3 = plot(
    label=["vapour" "liquid"], 
    [time_Plot,time_Plot],
    [p1_Plot,p2_Plot],
    xlabel = "Time [hours]",
    ylabel = "Pressure [bara]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

    display(p3)

    p4 = plot(
    label=["vapour" "liquid"], 
    [time_Plot,time_Plot],
    [T1_Plot,T2_Plot],
    xlabel = "Time [hours]",
    ylabel = "Temperaure [deg C]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

    display(p4)

    p5 = plot(
    label=["vapour" "liquid"],
    [time_Plot,time_Plot], 
    [m1_Plot,m2_Plot],
    xlabel = "Time [hours]",
    ylabel = "Mass [tons]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

    display(p5)

end

#pe = 8.0e5
#nstep = 3*19
delta_P = 0.25e5
pe = 10.0e5
nstep = (ps - pe)/delta_P
BlowDown(model, ps, Ts, zs, tank_volume, tank_radius, tank_length, level, pe, nstep)
