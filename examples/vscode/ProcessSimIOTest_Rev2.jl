using ModelingToolkit, DifferentialEquations, Clapeyron, NLsolve
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: scalarize
using Suppressor

fluid = ["carbon dioxide","nitrogen","water"]
nfluid = length(fluid)
maxphases = 3
FlashpTSize = 1 + 4*maxphases + 3 + maxphases*nfluid
model = GERG2008(fluid)

Molecular_Weight(model::EoSModel,z) = Clapeyron.molecular_weight(model::EoSModel,z)
@register_symbolic Molecular_Weight(model::EoSModel,z::AbstractVector)

Enthalpy(model::EoSModel,p,T,z) = Clapeyron.enthalpy(model::EoSModel,p,T,z)
@register_symbolic Enthalpy(model::EoSModel,p,T,z::AbstractVector)

Entropy(model::EoSModel,p,T,z) = Clapeyron.entropy(model::EoSModel,p,T,z)
@register_symbolic Entropy(model::EoSModel,p,T,z::AbstractVector)

Density(model::EoSModel,p,T,z) = Clapeyron.mass_density(model::EoSModel,p,T,z)
@register_symbolic Density(model::EoSModel,p,T,z::AbstractVector)

function FlashpT(model::EoSModel,p,T,z)
    verbose = false
    xp,beta,vp,Gp = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(verbose = verbose))
    res = fill(0.0,(FlashpTSize))
    ires = 1
    res[ires] = length(beta)
    for ip = eachindex(beta)
        ires += 1
        res[ires] =  beta[ip]
    end
    vt = 0.0
    for ip = eachindex(beta)
        ires += 1
        res[ires] =  vp[ip]
        vt += beta[ip]*vp[ip]
    end
    ires += 1
    res[ires] =  vt
    hp = Vector{Float64}(undef,0)
    for ip = eachindex(beta)
        push!(hp, Clapeyron.VT_enthalpy(model,vp[ip],T,xp[ip]))
    end
    ht = 0.0
    for ip = eachindex(beta)
        ires += 1
        res[ires] = hp[ip]
        ht += beta[ip]*hp[ip]
    end
    ires += 1
    res[ires] =  ht
    sp = Vector{Float64}(undef,0)
    for ip = eachindex(beta)
        push!(sp, Clapeyron.VT_entropy(model,vp[ip],T,xp[ip]))
    end
    st = 0.0
    for ip = eachindex(beta)
        ires += 1
        res[ires] = sp[ip]
        st += beta[ip]*sp[ip]
    end
    ires += 1
    res[ires] =  st
    for ip = eachindex(beta)
        for ic = eachindex(xp[ip])
            ires += 1
            res[ires] = xp[ip][ic]
        end
    end
    return res
end
@register_array_symbolic FlashpT(model::EoSModel,p,T,z::AbstractVector) begin
    size = (FlashpTSize,)
    eltype= eltype(z)
end

function Get_phaseBeta_FlashpT(res,phase)
    np = round(Int64, res[1])
    ires = 1
    beta = Vector{Float64}(undef,0)
    for ip = 1:np
        ires += 1
        push!(beta, res[ires])    
    end
     if phase > np
        return 0.0
    end
    return beta[phase]
end
@register_symbolic Get_phaseBeta_FlashpT(res::AbstractVector, phase::Int64)

function Get_phaseVolume_FlashpT(res,phase)
    np = round(Int64, res[1])
    ires = 1+np
    hp = Vector{Float64}(undef,0)
    for ip = 1:np
        ires += 1
        push!(hp, res[ires])
    end
    if phase > np
        return 0.0
    end
    return hp[phase]
end
@register_symbolic Get_phaseVolume_FlashpT(res::AbstractVector, phase::Int64)

function Get_Volume_FlashpT(res)
    np = round(Int64, res[1])
    ires = 1+np+np+1
    return res[ires]
end
@register_symbolic Get_Volume_FlashpT(res::AbstractVector)

function Get_phaseEnthalpy_FlashpT(res,phase)
    np = round(Int64, res[1])
    ires = 1+np+np+1
    hp = Vector{Float64}(undef,0)
    for ip = 1:np
        ires += 1
        push!(hp, res[ires])
    end
    if phase > np
        return 0.0
    end
    return hp[phase]
end
@register_symbolic Get_phaseEnthalpy_FlashpT(res::AbstractVector, phase::Int64)

function Get_Enthalpy_FlashpT(res)
    np = round(Int64, res[1])
    ires = 1+np+np+1+np+1
    return res[ires]
end
@register_symbolic Get_Enthalpy_FlashpT(res::AbstractVector)

function Get_phaseEntropy_FlashpT(res,phase)
    np = round(Int64, res[1])
    ires = 1+np+np+1+np+1
    hp = Vector{Float64}(undef,0)
    for ip = 1:np
        ires += 1
        push!(hp, res[ires])
    end
    if phase > np
        return 0.0
    end
    return hp[phase]
end
@register_symbolic Get_phaseEntropy_FlashpT(res::AbstractVector, phase::Int64)

function Get_Entropy_FlashpT(res)
    np = round(Int64, res[1])
    ires = 1+np+np+1+np+1+np+1
    return res[ires]
end
@register_symbolic Get_Entropy_FlashpT(res::AbstractVector)

function Get_phaseComposition_FlashpT(res,phase)
    np = round(Int64, res[1])
    ires = 1+np+np+1+np+1+np+1
    xp = Vector{Vector{Float64}}(undef,0)
    for ip = 1:np
        x = fill(0.0,nfluid)
        for ic = 1:nfluid
            ires += 1
            x[ic] = res[ires]
        end
        push!(xp, x)
    end
    if phase > np
        return fill(0.0,nfluid)
    end
    return xp[phase]
end
@register_array_symbolic Get_phaseComposition_FlashpT(res::AbstractVector, phase::Int64) begin
    size = (nfluid,)
    eltype= eltype(res)
end

function Flashph(model::EoSModel,T_in,h_in,z_in,p_out)
    function Flash(F,x)
        res_out = FlashpT(model,p_out,x[1],z_in)
        h_out = Get_Enthalpy_FlashpT(res_out)
        F[1] = h_out - h_in
    end
    T_out = nlsolve(Flash , [T_in])
    return T_out.zero[1]
end
@register_symbolic Flashph(model::EoSModel,T_in::Float64,h_in::Float64,z_in::Array{Float64, nfluid},p_out::Float64)

p = 1.4e5
T = 25.0+273.15
zdry=[0.9981, 0.0019]
zwater=0.025
z=append!(zdry*(1.0-zwater),zwater)

FlashResult = FlashpT(model,p,T,z)
println("FlashTp $(FlashResult)")

for ip = 1:round(Int64, FlashResult[1])
    global bpp = Get_phaseBeta_FlashpT(FlashResult,ip)
    println("FlashpT beta[$(ip)] = $(bpp)")
end

for ip = 1:round(Int64, FlashResult[1])
    global vpp = Get_phaseVolume_FlashpT(FlashResult,ip)
    println("FlashpT Phase Volume[$(ip)] = $(vpp)")
end

vtp = Get_Volume_FlashpT(FlashResult)
println("FlashpT Total Volume = $(vtp)")

for ip = 1:round(Int64, FlashResult[1])
    global hpp = Get_phaseEnthalpy_FlashpT(FlashResult,ip)
    println("FlashpT Phase Enthalpy[$(ip)] = $(hpp)")
end

htp = Get_Enthalpy_FlashpT(FlashResult)
println("FlashpT Total Enthalpy = $(htp)")

for ip = 1:round(Int64, FlashResult[1])
    global spp = Get_phaseEntropy_FlashpT(FlashResult,ip)
    println("FlashpT Phase Entropy[$(ip)] = $(spp)")
end

stp = Get_Entropy_FlashpT(FlashResult)
println("FlashpT Total Entropy = $(stp)")

for ip = 1:round(Int64, FlashResult[1])
    global xpp = Get_phaseComposition_FlashpT(FlashResult,ip)
    println("FlashpT Phase Composition[$(ip)] = $(xpp)")
end

p_in = 90.0e5
T_in = 35.0+273.15
zdry=[0.9981, 0.0019]
zwater=0.00001
z_in=append!(zdry*(1.0-zwater),zwater)

FlashResult2 = FlashpT(model,p_in,T_in,z_in)
println("FlashTp2 $(FlashResult2)")
h_in = Get_Enthalpy_FlashpT(FlashResult2)
println("FlashpT2 Total Enthalpy = $(h_in)")

p_out = 32.0e5
T_out = Flashph(model,T_in,h_in,z_in,p_out)
println("Flashph T out = $(T_out-273.15) deg C")

FlashResult3 = FlashpT(model,p_out,T_out,z_in)
println("FlashTp3 $(FlashResult3)")

@connector function fluidPort(;name, pg, Tg, zg)
    vars = @variables begin
        p(t), [input = true, description = "Pressure (Pa)", guess = pg]
        T(t), [input = true, description = "Temperature (K)", guess = Tg]
        v(t), [input = true, description = "Vomume (m3/mol)"]
        h(t), [input = true, description = "Enthalpy (J/kg)"]
        s(t), [input = true, description = "Entropy (J/kg/)"]
        mdot(t), [input = true, description = "mass flow rate (kg/s)"]
        (z(t))[1:nfluid], [input = true,description  ="Mole fraction vector", guess = zg]
    end
    ODESystem(Equation[], t, vars, [];name=name)
end

@component function MassSource(;name, ps, Ts, zs, mdots)
    @named port = fluidPort(pg = ps, Tg = Ts, zg = zs)
    vars = @variables begin   
        T(t)
        p(t)
        v(t)
        h(t)
        s(t)
        (z(t))[1:nfluid], [description  = "Mole fraction vector",  guess = zs]
        (flash(t))[1:FlashpTSize], [guess = FlashpT(model,ps,Ts,zs)]
    end

    para = @parameters begin
    #   T_source, [description = "Temperature at source (K)"]
    #   p_source, [description = "pressure at source (Pa)"]
    #   mdot_source, [description = "mass flow rate at source (kg/s)"]
     end

    eqs = [
        scalarize(flash .~ FlashpT(model,ps,Ts,zs))
        scalarize(z .~ zs)
        T ~ Ts
        p ~ ps
        v ~ Get_Volume_FlashpT(flash)
        h ~ Get_Enthalpy_FlashpT(flash)
        s ~ Get_Entropy_FlashpT(flash)

        port.T ~ Ts
        port.p ~ ps
        port.v ~ v
        port.h ~ h
        port.s ~ s
        port.mdot ~ mdots
        scalarize(port.z .~ zs)
    ]

#    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))
    compose(ODESystem(eqs, t, collect(Iterators.flatten(vars)), para;name=name),port)

end

@component function MassSink(;name, ps, Ts, zs)
    @named port = fluidPort(pg = ps, Tg = Ts, zg = zs)
    vars = @variables begin
       T(t)
       p(t)
       v(t)
       h(t)
       s(t)
       (z(t))[1:nfluid], [description  = "Mole fraction vector",  guess = zs]
       (flash(t))[1:FlashpTSize], [guess = FlashpT(model,ps,Ts,zs)]
    end
    para = @parameters begin

     end
     eqs = [
        scalarize(flash .~ FlashpT(model,port.p,port.T,port.z))
        scalarize(z .~ port.z)
        T ~ port.T
        p ~ port.p
        v ~ Get_Volume_FlashpT(flash)
        h ~ Get_Enthalpy_FlashpT(flash)
        s ~ Get_Entropy_FlashpT(flash)

     ]

#     eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))
     compose(ODESystem(eqs, t, collect(Iterators.flatten(vars)), para;name=name),port)
end

function PolytropicCompression(model::EoSModel,T_in,p_in,z,πc,η)
    @assert πc >= 1
    @assert η  <= 1
    @assert η  >  0
    mw = Molecular_Weight(model,z)
    d_in = Density(model,p_in,T_in,z)
    p_out = p_in*πc
    h_in = Enthalpy(model,p_in,T_in,z)/mw
    function Compressor(F,x)
        d_out = Density(model,p_out,x[1],z)
        h_out = Enthalpy(model,p_out,x[1],z)/mw
        npoly = log(p_out/p_in)/log(d_out/d_in)
        F[1] = h_out - h_in - npoly/(npoly-1)*(p_out/d_out - p_in/d_in)/η
    end
    Tpoly = nlsolve(Compressor , [T_in])
    return Tpoly.zero[1]
end
@register_symbolic PolytropicCompression(model::EoSModel,T_in::Float64,p_in::Float64,z::Array{Float64, nfluid},πc::Float64,η::Float64)

@component function Compressor(;name, ps, Ts, zs, pd, Td)
    @named inport = fluidPort(pg = ps, Tg = Ts, zg = zs)
    @named outport = fluidPort(pg = pd, Tg = Td, zg = zs)
    vars = @variables begin
       s_in(t)
       T_in(t)
       p_in(t)
       h_in(t)
       (z_in(t))[1:nfluid]

       s_out(t)
       T_out(t)
       p_out(t)
       h_out(t)
       (z_out(t))[1:nfluid]

    end
    para = @parameters begin
        πc = 4.528, [description = "Pressure ratio (-)"]
        η = 0.84, [description = "Polytropic Efficiency (-)"]
    end
    eqs = [

        scalarize(z_in .~ inport.z)
        T_in ~ inport.T
        p_in ~ inport.p
        h_in ~ Enthalpy(model,p_in,T_in,z_in)
        s_in ~ Entropy(model,p_in,T_in,z_in)

        scalarize(z_out .~ z_in)
        p_out ~ πc*p_in
        T_out ~ PolytropicCompression(model,T_in,p_in,z_in,πc,η)
        s_out ~ Entropy(model,p_out,T_out,z_out)
        h_out ~ Enthalpy(model,p_out,T_out,z_out)

        outport.T ~ T_out
        outport.p ~ p_out
        outport.mdot ~ inport.mdot
        scalarize(outport.z .~ z_out)
    ]

#    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))    
    compose(ODESystem(eqs, t, collect(Iterators.flatten(vars)), para;name=name),inport,outport)

end

@component function Separator(;name, pg, Tg, zg)
    @named inport = fluidPort(pg = pg, Tg = Tg, zg = zg)
    @named outport1 = fluidPort(pg = pg, Tg = Tg, zg = zg) # on top
    @named outport2 = fluidPort(pg = pg, Tg = Tg, zg = zg) # on bottom
    vars = @variables begin
       mdot_in(t)
       s_in(t)
       T_in(t)
       p_in(t)
       h_in(t)
       (z_in(t))[1:nfluid]

       mdot_out1(t)
       s_out1(t)
       T_out1(t)
       p_out1(t)
       h_out1(t)
       (z_out1(t))[1:nfluid]

       mdot_out2(t)
       s_out2(t)
       T_out2(t)
       p_out2(t)
       h_out2(t)
       (z_out2(t))[1:nfluid]
       (flash(t))[1:FlashpTSize]

    end
    para = @parameters begin end
    eqs = [

        mdot_in ~ inport.mdot

        scalarize(flash .~ FlashpT(model,p_in,T_in,z_in))
        scalarize(z_in .~ inport.z)
        T_in ~ inport.T
        p_in ~ inport.p
        h_in ~ Get_Enthalpy_FlashpT(flash)
        s_in ~ Get_Entropy_FlashpT(flash)
    
    #   need to decide on number of phase we are hard coding 2 here we could have 1 or even 3
    #   we know enthalp and entropy of output but it not used yet would be needed for pure component as there are not TPFlash
        scalarize(z_out1 .~ Get_phaseComposition_FlashpT(flash,1))
        p_out1 ~ p_in
        T_out1 ~ T_in
        h_out1 ~ Get_phaseEnthalpy_FlashpT(flash,1)
        s_out1 ~ Get_phaseEntropy_FlashpT(flash,1)

    #   here if only 1 phase composition should be phase 1
    #   need a varable for nphases
        scalarize(z_out2 .~ Get_phaseComposition_FlashpT(flash,2))
        p_out2 ~ p_in
        T_out2 ~ T_in
        h_out2 ~ Get_phaseEnthalpy_FlashpT(flash,2)
        s_out2 ~ Get_phaseEntropy_FlashpT(flash,2)
        
        outport1.T ~ T_out1
        outport1.p ~ p_out1
        mdot_out1 ~ mdot_in*Get_phaseBeta_FlashpT(flash,1) 
        outport1.mdot ~ mdot_out1
        scalarize(outport1.z .~ z_out1)

        outport2.T ~ T_out2
        outport2.p ~ p_out2
        mdot_out2 ~ mdot_in*Get_phaseBeta_FlashpT(flash,2)
        outport2.mdot ~ mdot_out2
        scalarize(outport2.z .~ z_out2)
   
    ]

#    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))
    compose(ODESystem(eqs, t, collect(Iterators.flatten(vars)), para;name=name),inport,outport1,outport2)
end

p = 1.4e5
T = 25.0+273.15
zdry=[0.9981, 0.0019]
zwater=0.025
z=append!(zdry*(1.0-zwater),zwater)

@suppress_err begin

    local src, separator, compressor, sink, eqs, systems, flowsheet, sys, u0, para, prob

    @named src = MassSource(ps = p, Ts = T, mdots = 50, zs = z)
    @named separator = Separator(pg = p, Tg = T, zg = z)
    @named compressor = Compressor(ps = p, Ts = T, zs = z, pd = 4.528*p, Td = 161.0+273.15)
    @named sink =  MassSink(ps = p, Ts = T, zs = z)

    eqs =   [
            connect(src.port,separator.inport)
            connect(separator.outport1,compressor.inport)
            connect(compressor.outport,sink.port)
            ] # Define connections

    systems=[src,separator, compressor,sink] # Define system

    @named flowsheet = System(eqs, t, systems=systems)
    sys = structural_simplify(flowsheet)
    u0 = []
    para =  [   
                sys.compressor.πc =>4.528,
                sys.compressor.η =>0.84 
            ] # boundary conditions and parameters


    prob = SteadyStateProblem(sys,u0,para)
    sol = solve(prob)

    println("Compressor discharge temperature = $(round(sol[compressor.T_out][1]-273.15,digits=2)) deg C")
    comp_h_in = sol[compressor.h_in][1]
    println("Compressor h_in = $(round(comp_h_in,digits=2)) J/mol")
    comp_h_out = sol[compressor.h_out][1]
    println("Compressor h_out = $(round(comp_h_out,digits=2)) J/mol")
    comp_mdot = sol[compressor.inport.mdot][1]
    println("Compressor mdot = $(round(comp_mdot,digits=2)) kg/s")
    comp_mw =  Molecular_Weight(model,sol[compressor.z_out])
    println("Compressor mw = $(round(comp_mw*1000,digits=2)) kg/kmol")
    power = comp_mdot*(comp_h_out - comp_h_in)/comp_mw/1000.0/1000.0
    println("Compressor power = $(round(power,digits=2)) MW")

end