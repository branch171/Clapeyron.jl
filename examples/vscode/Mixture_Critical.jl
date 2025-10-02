using Clapeyron, Plots

model1 = GERG2008(["carbon dioxide","nitrogen"])
model2 = PCSAFT(["carbon dioxide","nitrogen"])

#model3 = GERG2008(["carbon dioxide","ethane"])
#model4 = SAFTVRMie(["carbon dioxide","ethane"])
#model5 = sPCSAFT(["carbon dioxide","ethane"])
#models = [model1,model2,model3,model4,model5]

models = [model1,model2]

nmodels = length(models)

Trs = 0.95
p_1 = zeros(100,nmodels)
p_2 = zeros(100,nmodels)
T_1 = zeros(100,nmodels)
T_2 = zeros(100,nmodels)
vl_1 = zeros(100,nmodels)
vl_2 = zeros(100,nmodels)
vv_1 = zeros(100,nmodels)
vv_2 = zeros(100,nmodels)
for j=1:nmodels
    pure = split_model(models[j])
    crit = crit_pure.(pure)
    T = zeros(100,2)
    p = zeros(100,2)
    vl = zeros(100,2)
    vv = zeros(100,2)
    sat = []
    for i=1:2
        T[:,i] = range(Trs*crit[i][1],crit[i][1],length=100)
        sat = saturation_pressure.(pure[i],T[:,i])
        p[:,i] = [sat[i][1] for i in 1:100]
        vl[:,i] = [sat[i][2] for i in 1:100]
        vv[:,i] = [sat[i][3] for i in 1:100]
    end

    p_1[:,j] = p[:,1]./1e5
    p_2[:,j] = p[:,2]./1e5

    T_1[:,j] = T[:,1].-273.15
    T_2[:,j] = T[:,2].-273.15

    vl_1[:,j] = vl[:,1]
    vl_2[:,j] = vl[:,2]

    vv_1[:,j] = vl[:,1]
    vv_2[:,j] = vl[:,2]
end

Tc = zeros(nmodels)
Tc[1] = T_1[100,1]+273.15
Tc[2] = T_1[100,2]+273.15

pc = zeros(nmodels)
pc[1] = p_1[100,1]*1e5
pc[2] = p_1[100,2]*1e5

vc = zeros(nmodels)
vc[1] = vl_1[100,1]
vc[2] = vl_1[100,2]

npc = 250
x = range(1-1e-6,0.50,length=npc)
X = Clapeyron.FractionVector.(x)
T_crit = zeros(npc,nmodels)
p_crit = zeros(npc,nmodels)

#v0 = Clapeyron.x0_crit_mix(models[1],X[1])
#println("v0 = $(v0)")
#mix_crit0 = crit_mix(models[1],X[1];v0=v0)
#println("mix_crit0 = $(mix_crit0)")

for j=1:nmodels
#    v0 = Clapeyron.x0_crit_mix(models[j],X[1])
    v0 = [log10(vc[j]),Tc[j]]
    for i=1:npc
        mix_crit = crit_mix(models[j],X[i];v0=v0)
        v0 = [log10(mix_crit[3]),mix_crit[1]]
        
        T_crit[i,j] = mix_crit[1]-273.15
        p_crit[i,j] = mix_crit[2]/1e5
    end
end

println("T_crit = $(T_crit)")
println("p_crit = $(p_crit)")

p1 = plot([T_1,T_2,T_crit], [p_1,p_2,p_crit],
label=["GERG2008" "SAFTVRMie" "GERG2008" "SAFTVRMie" "GERG2008" "SAFTVRMie"], 
xlabel = "Temperature [deg C]", 
ylabel = "Pressure [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p1)
#=
p2 = plot([T_crit], [p_crit],
label=["GERG2008" "SAFTVRMie" "GERG2008" "SAFTVRMie" "GERG2008" "SAFTVRMie"], 
xlabel = "Temperature [deg C]", 
ylabel = "Pressure [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p2)
=#