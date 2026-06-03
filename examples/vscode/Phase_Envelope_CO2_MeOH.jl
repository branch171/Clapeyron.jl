using Clapeyron, Plots

model = sCPA(["carbon dioxide","methanol"])
#model = SAFTVRMie(["carbon dioxide","nitrogen"])

xCO2 = 0.9975
x = [xCO2,1.0-xCO2]

idxstart = 10000

Tbub = zeros(idxstart)
pbub = zeros(idxstart)

idxend = idxstart

Tbub[1] = 240.
bub = bubble_pressure(model,Tbub[1],x)

pbub[1] = bub[1]
dT = 0.1
v0p = vcat(log10(bub[2]),log10(bub[3]),bub[4])
for j in 2:2:idxstart-1
    global dT,v0p,idxend,bub
    Tbub[j] = Tbub[j-1] + dT
    bub = bubble_pressure(model,Tbub[j],x;v0=v0p)
    if isnan(bub[1])
        idxend = j-1
        break
    end
    pbub[j] = bub[1]

    if abs(bub[3]-bub[2]) < 1e-8
        idxend = j-2
        break
    end

    dp = pbub[j]-pbub[j-1]
    pbub[j+1] = pbub[j] + dp
        
    v0T = vcat(Tbub[j],log10(bub[2]),log10(bub[3]),bub[4])
    bub = bubble_temperature(model,pbub[j+1],x;v0=v0T)
    Tbub[j+1] = bub[1]
    dT = Tbub[j+1]-Tbub[j]

    if abs(dT)<1e-3
        dT*=2.0
    end

    if abs(bub[3]-bub[2]) < 1e-8
        idxend = j-1
        break
    end

    v0p = vcat(log10(bub[2]),log10(bub[3]),bub[4])
end

Tbub = Tbub[1:idxend-1]
pbub = pbub[1:idxend-1]

idxstart = 10000

Tdew = zeros(idxstart)
pdew = zeros(idxstart)

idxend = idxstart

Tdew[1] = 240.
dew = dew_pressure(model,Tdew[1],x)
pdew[1] = dew[1]
dT = 0.1
v0p = vcat(log10(dew[2]),log10(dew[3]),dew[4])
for j in 2:2:idxstart-1
    global dT,v0p,idxend,dew
    Tdew[j] = Tdew[j-1] + dT
    dew = dew_pressure(model,Tdew[j],x;v0=v0p)
    if isnan(dew[1])
        idxend = j-1
        break
    end
    pdew[j] = dew[1]

    if abs(dew[3]-dew[2]) < 1e-8
        idxend = j-2
        break
    end

    dp = pdew[j]-pdew[j-1]
    pdew[j+1] = pdew[j] + dp
        
    v0T = vcat(Tdew[j],log10(dew[2]),log10(dew[3]),dew[4])
    dew = dew_temperature(model,pdew[j+1],x;v0=v0T)
    Tdew[j+1] = dew[1]
    dT = Tdew[j+1]-Tdew[j]

    if abs(dT)<1e-3
        dT*=2.0
    end

    if abs(dew[3]-dew[2]) < 1e-8
        idxend = j-1
        break
    end
    v0p = vcat(log10(dew[2]),log10(dew[3]),dew[4])
end

Tdew = Tdew[1:idxend-1]
pdew = pdew[1:idxend-1]

println("Tcrit = $((Tbub[end] + Tdew[end])/2.0)")
println("pcrit = $((pbub[end] + pdew[end])/1e5/2.0)")

T = vcat(Tbub,reverse(Tdew))
p = vcat(pbub,reverse(pdew))

TC = T .- 273.15
pbar = p ./ 1e5

p5 = plot(TC, pbar,
label=["GERG2008" ], 
xlabel = "Temperature [deg C]", 
ylabel = "Pressure [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p5)