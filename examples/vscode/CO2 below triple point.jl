using Clapeyron, Plots

fluid = ["carbon dioxide","nitrogen"]
nfluid = length(fluid)
#model = PR(fluid)
#model = PCSAFT(fluid)
#model = SAFTVRMie(fluid)
model = GERG2008(fluid)
#model = EOS_CG(fluid)

pbase = 27e5
modelCO2 = GERG2008(["carbon dioxide"])
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

pbubble = pbase
mix_Tbub = bubble_temperature(model,pbubble,x)

println("T bubble = $(mix_Tbub[1]-273.15)")

pdew = pbase
mix_Tdew = dew_temperature(model,pdew,x)

println("T dew = $(mix_Tdew[1]-273.15)")

Tin = -78.5+273.15
pin = 1.01325e5*0.9998965
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

R̄ = Clapeyron.R̄

println("Volume Total = $(vt) m3/mol")
println("Volume Total Ideal = $((R̄*Tin)/pin) m3/mol")
for ip = eachindex(βf)
    println("Volume Fraction Vapour = $(βf[ip]*vf[ip]/vt)")
end

Z = pin*vt/(R̄*Tin)
μ_res = Clapeyron.VT_chemical_potential_res(model, vt, Tin, zin)
ϕ = exp.(μ_res/(R̄*Tin) .- log(Z))

println("Z = $(Z)")
println("μ_res/RT = $(μ_res/(R̄*Tin))")
println("ϕ = $(ϕ)")

function psub(T)
    # CO2 T in K psub on Pa
    pt = 5.17958e5
    Tt = 216.592
    a1 = -14.740846
    a2 =   2.4327015
    a3 =  -5.3061778
    return pt*exp(Tt/T*(a1*(1 - T/Tt) + a2*(1 - T/Tt)^1.9 + a3*(1 - T/Tt)^2.9))
end

ppure = psub(Tin)

println("ppure = $(ppure/1e5) bara")

vpure = volume(modelCO2, ppure, Tin)

println("vpure = $(vpure) m3/mol")

μpure = Clapeyron.VT_chemical_potential_res(modelCO2, vpure, Tin, [1])

println("μpure = $(μpure)")

Zpure = ppure*vpure/(R̄*Tin)
ϕpure = exp.(μpure/(R̄*Tin) .- log(Zpure))

println("Zpure = $(Zpure)")
println("μpure/RT = $(μpure/(R̄*Tin))")
println("ϕpure = $(ϕpure)")

mws = Clapeyron.molecular_weight(modelCO2,[1])
println("mw solid = $(mws) kg/mol")
vs = 1.0/1560.0*mws
println("fvolume solid = $(vs) m3/mol")
fs = ϕpure[1]*ppure*exp(vs/(R̄*Tin)*(pin - ppure))
fv = zin[1]*ϕ[1]*pin

println("fugacity solid = $(fs/1e5) bara")
println("fugacity vapour = $(fv/1e5) bara")

println("error = $(fs-fv) Pa")

#=
g42.3945101 103
9.2923982 103g52.6531689 103
g103.3914617 103g 61.6419734 101
2.7976237 100n7.0000000 100g 71.7594802 101
g42.6427834 101g 03.9993365 102g82.6531689 103
g53.8259935 100g12.3945101 103g 02.2690751 101
g63.1711996 101g 23.2839467 101g17.5019750 102
g72.2087195 103g 35.7918471 102g 22.6442913 101
g0 2.6385478 100
g1
Value
g81.1289668 100
4.5088732 100g9g22.0109135 100g3
=#



Ttp = 216.592
ptp = 5.1795e5

sat = saturation_temperature(modelCO2, ptp)
Ttp = sat[1]
vltp = sat[2]
vvtp = sat[3]

println("T tp = $(Ttp) K")
println("T tp = $(Ttp - 273.15) degC")

println("v liquid tp = $(vltp) m3/mol")
println("v vapour tp = $(vvtp) m3/mol")

hltp = Clapeyron.VT_enthalpy(modelCO2, vltp,Ttp,[1])
hvtp = Clapeyron.VT_enthalpy(modelCO2, vvtp,Ttp,[1])
println("h liquid tp = $(hltp) J/mol")
println("h liquid tp = $(hltp/0.0440098/1000) kJ/kg")
println("h vapour tp = $(hvtp) J/mol")
println("h vapour tp = $(hvtp/0.0440098/1000) kJ/kg")

println("h vaporisation tp = $(hvtp - hltp) J/mol")
println("h vaporisation tp = $((hvtp - hltp)/0.0440098/1000) kJ/kg")

sltp = Clapeyron.VT_entropy(modelCO2, vltp,Ttp,[1])
svtp = Clapeyron.VT_entropy(modelCO2, vvtp,Ttp,[1])
println("s liquid tp = $(sltp) J/mol/K")
println("s liquid tp = $(sltp/0.0440098/1000) kJ/kg/K")
println("s vapour tp = $(svtp) J/mol/K")
println("s vapour tp = $(svtp/0.0440098/1000) kJ/kg/K")

gltp = hltp - Ttp*sltp
gvtp = hvtp - Ttp*svtp
println("g liquid tp = $(gltp) J/mol")
println("g vapour tp = $(gvtp) J/mol")


function gibbs_dryice(p,T,p0,T0,h0,s0)
    R̄  = 8.314462
    #T0  = 216.59242086
    #p0  = 5.1795e5
    #g0  = -2.6385478e+0
    hmelt = 8875.0
    g0  = (h0 - T0*s0)/R̄/T0
    #g1  =  4.5088732e+0
    g1 = -(s0 - hmelt/T0)/R̄
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

function entropy_dryice(p,T,p0,T0,h0,s0)
    g(x) = gibbs_dryice(p,x,p0,T0,h0,s0)
    s(x) = Clapeyron.Solvers.derivative(g,x)
    return -s(T)
end

gtp = gibbs_dryice(ptp,Ttp,ptp,Ttp,hltp,sltp)

println("gtp = $(gtp) J/mol")
println("gtp = $(gtp/0.0440098/1000) kJ/kg")

stp = entropy_dryice(ptp,Ttp,ptp,Ttp,hltp,sltp)

println("stp = $(stp) J/mol/K")
println("stp = $(stp/0.0440098/1000) kJ/kg/K")

htp = gtp + Ttp*stp

println("htp = $(htp) J/mol")
println("htp = $(htp/0.0440098/1000) kJ/kg")

gwLiq = (-426.74 + Ttp*2.2177)

println("gwLIQ = $(gwLiq) KJ/kg")

gwVap = (-76.364 + Ttp*0.59999)

println("gwVap = $(gwVap) KJ/kg")

sstp = sltp - 8875/Ttp

println("sstp = $(sstp) J/mol/K")

# check p = 101325 Pa

println("")
println("check")
Tin = -78.75+273.15
pin = 1.01325e5
xs = 0.64275
zin = [(x[1] - xs)/(1 - xs), x[2]/(1 - xs)]

verbose = false
flash_result = Clapeyron.tp_flash_impl(model,pin,Tin,zin, HELDTPFlash(verbose = verbose))

Tf = flash_result.data.T
println("T flash = $(Tf-273.15)")
βf = flash_result.fractions
println("βf = $(βf)")
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

hf = Clapeyron.VT_enthalpy(model, vt,Tf,zin)
println("hf = $(hf) J/mol/K")

sf = Clapeyron.VT_entropy(model, vt,Tf,zin)
println("sf = $(sf) J/mol/K")

gf = hf - Tf*sf
println("gf = $(gf) J/mol")

μf = Clapeyron.VT_chemical_potential(model, vf[1], Tf, xf[1])
println("μf = $(μf) J/mol")

μfmix = xf[1][1]*μf[1] + xf[1][2]*μf[2]
println("μfmix = $(μfmix) J/mol")

gs = gibbs_dryice(pin,Tin,ptp,Ttp,hltp,sltp)
println("gs = $(gs) J/mol")

gmix = xs*gs + (1 - xs)*gf
println("gmix = $(gmix) J/mol")

function gibbsmix(model,p,T,z,xs,iCO2)
    zf = zeros(length(z))
    zf[iCO2] = z[iCO2] - xs
    for i=eachindex(z)
        if i ≠ iCO2
            zf[i] = z[i]
        end
    end
    sumzf = sum(zf)
    for i=eachindex(z)
        zf[i] /= sumzf
    end
    f = Clapeyron.tp_flash_impl(model,p,T,zf, HELDTPFlash(verbose = false))
    βf = f.fractions
    vf = f.volumes
    xf = f.compositions
    vt = 0.0
    for ip = eachindex(βf)
        vt += βf[ip]*vf[ip]
    end
    hf = zeros(length(z))
    for ip = eachindex(βf)
        hf[ip] = Clapeyron.VT_enthalpy(model,vf[ip],T,xf[ip])
    end
    ht = 0.0
    for ip = eachindex(βf)
        ht += βf[ip]*hf[ip]
    end
    sf = zeros(length(z))
    for ip = eachindex(βf)
        sf[ip] = Clapeyron.VT_entropy(model,vf[ip],T,xf[ip])
    end
    st = 0.0
    for ip = eachindex(βf)
        st += βf[ip]*sf[ip]
    end
    gf = zeros(length(z))
    for ip = eachindex(βf)
        gf[ip] = hf[ip] - T*sf[ip]
    end
    gt = 0.0
    for ip = eachindex(βf)
        gt += βf[ip]*gf[ip]
    end
    gs = gibbs_dryice(p,T,ptp,Ttp,hltp,sltp)
    gm = xs*gs + (1-xs)*gt
    return gm
end

gm  = gibbsmix(model,pin,Tin,x,xs,1)
println("gm = $(gm) J/mol")

npoint = 21
xsmin = zeros(npoint)
gmmin = zeros(npoint)
xsolid = zeros(npoint)
gmixture = zeros(npoint)
for i = 1:npoint
    xsmin[i] = xs
    gmmin[i] = -796.33 + 0.04*(i-1)/(npoint-1)
    xsolid[i] = 0.6 + 0.075*(i-1)/(npoint-1)
    gmixture[i] = gibbsmix(model,pin,Tin,x,xsolid[i],1)
end

p1 = plot(
    [xsolid,xsmin], 
    [gmixture,gmmin], 
    xlabel = "xs", 
    ylabel = "gm",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))
display(p1)

function brentmin(
    f,
    x_lower::T,
    x_upper::T,
    rel_tol::T = sqrt(eps(T)),
    abs_tol::T = eps(T),
    iterations::Integer = 1_000,
) where {T<:AbstractFloat}

    if x_lower > x_upper
        error("x_lower must be less than x_upper")
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
    #    println("x_lower = $(x_lower), new_minimizer = $(new_minimizer), x_upper = $(x_upper)")
    end

    return new_minimizer, new_minimum

end

function pT_flash_CO2Solid(model,p,T,z)
    iCO2 = 0
    for i = eachindex(model.components)
        if model.components[i] == "carbon dioxide"
            iCO2 = i
        end
    end
    @assert iCO2 > 0 "carbon dioxide must me present in components"
    xs_eps = eps(Float64)
    xs_min = 0.0
    xs_max = z[iCO2]*(1 - xs_eps)
    f = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(verbose = false))
    βf = f.fractions
    vf = f.volumes
    xf = f.compositions
    imax = 1
    if length(βf) > 1
        βmax = 0.0
        for i = eachindex(βf)
            if βf[i] > βmax
                βmax = βf[i]
                imax = i
            end
        end
    end
    μf = Clapeyron.VT_chemical_potential(model, vf[imax], T, xf[imax])
    μs = gibbs_dryice(p,T,ptp,Ttp,hltp,sltp)
    if μs < μf[iCO2]
        fnc(xs) = gibbsmix(model,p,T,z,xs,iCO2)
        res = brentmin(fnc, xs_min, xs_max)
        xs = res[1]
        zf = zeros(length(z))
        zf[iCO2] = z[iCO2] - xs
        for i=eachindex(z)
            if i ≠ iCO2
                zf[i] = z[i]
            end
        end
        sumzf = sum(zf)
        for i=eachindex(z)
            zf[i] /= sumzf
        end
        f = Clapeyron.tp_flash_impl(model,p,T,zf, HELDTPFlash(verbose = false))
        return xs,f
    else
        return xsmin,f
    end
end

βCO2, flash_fluid = pT_flash_CO2Solid(model,pin,Tin,x)

println("moles of solid CO2 = $(βCO2)")
println("moles of fluid = $(1 - βCO2)")

println("fluid βf = $(flash_fluid.fractions)")
println("fluid vf = $(flash_fluid.volumes) m3/mol")
for i = eachindex(flash_fluid.fractions)
    println("fluid xf[$(i)] = $(flash_fluid.compositions[i])")
end


