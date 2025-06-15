using Clapeyron, Plots

model1 = GERG2008(["carbon dioxide","nitrogen"])
model2 = SAFTVRMieVTC(["carbon dioxide","nitrogen"])

models = [model1,model2];

p_1 = zeros(100,2)
p_2 = zeros(100,2)
T_1 = zeros(100,2)
T_2 = zeros(100,2)
for j=1:2
    pure = split_model(models[j])
    crit = crit_pure.(pure)
    T = zeros(100,2)
    p = zeros(100,2)
    sat = []
    for i=1:2
        T[:,i] = range(0.7*crit[i][1],crit[i][1],length=100)
        sat = saturation_pressure.(pure[i],T[:,i])
        p[:,i] = [sat[i][1] for i in 1:100]
    end

    p_1[:,j] = p[:,1]
    p_2[:,j] = p[:,2]

    T_1[:,j] = T[:,1]
    T_2[:,j] = T[:,2]
end

N = 20
x = range(1e-3,1.0-1e-3,length=N)
X = Clapeyron.FractionVector.(x)
T_crit = zeros(N,2)
p_crit = zeros(N,2)
for j=1:2
    T_crit[1,j] = T_2[end,j]
    p_crit[1,j] = p_2[end,j]
    T_crit[N,j] = T_1[end,j]
    p_crit[N,j] = p_1[end,j]
    v0 = Clapeyron.x0_crit_mix(models[j],X[1])
    for i=2:N-1
        mix_crit = crit_mix(models[j],X[i];v0=v0)
        v0 = [log10(mix_crit[3]),mix_crit[1]]
        
        T_crit[i,j] = mix_crit[1]
        p_crit[i,j] = mix_crit[2]
    end
end

#println("T_crit GERG2008 = $(T_crit[:,1])")
#println("T_crit SAFTVRMieVTC = $(T_crit[:,2])")

p1 = plot([T_1[:,1],T_1[:,2],T_2[:,1],T_2[:,2],T_crit[:,1],T_crit[:,2]], [p_1[:,1],p_1[:,2],p_2[:,1],p_2[:,2],p_crit[:,1],p_crit[:,2]],
label=["GERG2008 CO2" "SAFTVRMieVTC CO2" "GERG2008 N2" "SAFTVRMieVTC N2" "GERG2008 Critcal" "SAFTVRMieVTC Critical"], 
xlabel = "T [K]", 
ylabel = "p [Pa]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p1)