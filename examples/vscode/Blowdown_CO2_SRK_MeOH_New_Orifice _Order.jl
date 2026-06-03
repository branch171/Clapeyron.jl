using Clapeyron, Ferrite, FerriteGmsh, Gmsh, SparseArrays, WriteVTK, Plots

model_fluid = ["methanol","nitrogen","carbon dioxide"]
nfluid = length(model_fluid)
#model = sCPA(model_fluid; idealmodel=AlyLeeIdeal, assoc_options=AssocOptions(combining=:elliott))
model = SRK(model_fluid; idealmodel=AlyLeeIdeal)
pures_model = Clapeyron.split_pure_model(model)
crit_model = Clapeyron.crit_pure.(pures_model)

locs = ["properties/molarmass.csv","properties/critical.csv"]
model_params = Clapeyron.getparams(model_fluid, locs; userlocations = String[], verbose = false)
#println("model_params = $(model_params)")

Vc = model_params["Vc"].values
for i = 1:nfluid
	Tc,pc,vc = crit_model[i]
    println("vc = $(vc)")
    println("Vc = $(Vc[i])")
end

#modelCO2 = sCPA(["carbon dioxide"]; idealmodel=AlyLeeIdeal)
modelCO2 = SRK(["carbon dioxide"]; idealmodel=AlyLeeIdeal)
mwCO2 = Clapeyron.molecular_weight(modelCO2,[1])

#modelN2 = sCPA(model_fluid; idealmodel=AlyLeeIdeal, assoc_options=AssocOptions(combining=:elliott))
modelN2 = SRK(model_fluid; idealmodel=AlyLeeIdeal)
zN2 = [0.0005,0.0005,0.999]
mwN2 = Clapeyron.molecular_weight(modelN2,zN2)

modelAir_fluid = ["nitrogen","oxygen"]
modelAir = GERG2008(modelAir_fluid)
zAir = [0.79,0.21]
mwAir = Clapeyron.molecular_weight(modelAir,zAir)

modelAir_params = Clapeyron.getparams(modelAir_fluid, locs; userlocations = String[], verbose = false)
#println("modelAir_params = $(modelAir_params)")
pures_modelAir = Clapeyron.split_pure_model(modelAir)

modelC3H6_params = Clapeyron.getparams(["propane"], locs; userlocations = String[], verbose = false)
#println("modelC3H6_params = $(modelC3H6_params)")

#modelMeOH = sCPA(model_fluid; idealmodel=AlyLeeIdeal, assoc_options=AssocOptions(combining=:elliott))
modelMeOH = SRK(model_fluid; idealmodel=AlyLeeIdeal)
zMeOH = [0.999,0.0005,0.0005]
mwMeOH = Clapeyron.molecular_weight(modelMeOH,zMeOH)

z_solid = [0.0,0.0,1.0]

function SuperTRAPP_mu(model_params, modelC3H6_params, Tin, rhoin, zin)

    R̄ = Clapeyron.R̄

	E = [
	-14.113294896,
	968.22940153,
	13.686545032,
	-12511.628378,
	0.0168910864,
	43.527109444,
	7659.4543472]

	T = Tin
	rho = rhoin/1000.0 # kmol/m3
	z = zin

	nc = length(z)

    #=
	modelC3H6 = GERG2008(["propane"])

	MwRef = modelC3H6.params.Mw.values[1]
	TcRef = modelC3H6.params.Tc.values[1]
	PcRef = modelC3H6.params.Pc.values[1]
	VcRef = modelC3H6.params.Vc.values[1]*100.0^3 # cm3/mol
	rhocRef = 1000.0/VcRef
	OmegaRef = modelC3H6.params.acentricfactor.values[1]
	ZcRef = PcRef*VcRef/(100.0^3)/R̄/TcRef

	println("Mw, Tc, Pc, Vc, Omega = $(MwRef), $(TcRef), $(PcRef), $(VcRef), $(OmegaRef), $(ZcRef)")
    =#


	MwRef = modelC3H6_params["Mw"].values[1]
	TcRef = modelC3H6_params["Tc"].values[1]
	PcRef = modelC3H6_params["Pc"].values[1]
	VcRef = modelC3H6_params["Vc"].values[1]*100.0^3 # cm3/mol
	rhocRef = 1000.0/VcRef
	OmegaRef = modelC3H6_params["acentricfactor"].values[1]
	ZcRef = PcRef*VcRef/(100.0^3)/R̄/TcRef

#	println("Mw, Tc, Pc, Vc, Omega = $(MwRef), $(TcRef), $(PcRef), $(VcRef), $(OmegaRef), $(ZcRef)")

    #=
	Mw = model.params.Mw.values
	Tc = model.params.Tc.values
	Pc = model.params.Pc.values
	Vc = model.params.Vc.values*100.0^3 # cm3/mol
	Zc = zeros(nc)
	for i = eachindex(z)
		Zc[i] = Pc[i]*Vc[i]/(100.0^3)/R̄/Tc[i]
    end
	Omega = model.params.acentricfactor.values

	println("Mw, Tc, Pc, Vc, Omega = $(Mw), $(Tc), $(Pc), $(Vc), $(Omega), $(Zc)")
    =#

    Mw = model_params["Mw"].values
	Tc = model_params["Tc"].values
	Pc = model_params["Pc"].values
	Vc = model_params["Vc"].values*100.0^3 # cm3/mol
	Zc = zeros(nc)
	for i = eachindex(z)
		Zc[i] = Pc[i]*Vc[i]/(100.0^3)/R̄/Tc[i]
    end
	Omega = model_params["acentricfactor"].values

#	println("Mw, Tc, Pc, Vc, Omega = $(Mw), $(Tc), $(Pc), $(Vc), $(Omega), $(Zc)")

	Mw_mix = 0.0
	for i = eachindex(z)
		Mw_mix += z[i] * Mw[i]
    end

	f = zeros(nc)
	for i = eachindex(z)
		f[i] = Tc[i] / TcRef * (1.0 + (Omega[i] - OmegaRef) * (0.05202976 - 0.7498189 * log(T / Tc[i])))
	end

#	println("f = $(f)")

	h = zeros(nc)
	for i = eachindex(z)
		h[i] = Vc[i] / VcRef * ZcRef / Zc[i] * (1.0 + (Omega[i] - OmegaRef) * (0.1435971 - 0.2821562 * log(T / Tc[i])))
	end

#	println("h = $(h)")

	hm = 0.0;
	for i = eachindex(z)
		for j = eachindex(z)
			temp = h[i]^(1/3) + h[j]^(1/3)
			temp = temp * temp * temp / 8.0
			hm += z[i] * z[j] * temp
		end
	end

#	println("hm = $(hm)")

	fm = 0.0
	for i = eachindex(z)
		for j = eachindex(z)
			temp = h[i]^(1/3) + h[j]^(1/3)
			temp = temp * temp * temp / 8.0
			fm += z[i] * z[j] * temp * sqrt(f[i] * f[j])
		end
	end
	fm = fm / hm

#	println("fm = $(fm)")

	To = T / fm;
	rhoo = rho * hm;

	G1 = exp(E[1] + E[2] / To)
	temp = To^1.5
	G2 = E[3] + E[4] / temp;
	G3 = E[5] + E[6] / To + E[7] / To / To

#	println("G1 = $(G1)")
#	println("G2 = $(G2)")
#	println("G3 = $(G3)")

	DnetaR = G1 * (exp((rhoo^0.1) * G2 + sqrt(rhoo) * (rhoo / rhocRef - 1.0) * G3) - 1.0)

#	println("DnetaR = $(DnetaR)")

	s = zeros(nc)
	for i = eachindex(z)
		s[i] = 4.771 * h[i]^(1/3)
	end

	ss2 = 0.0
	for i = eachindex(z)
		ss2 += z[i] * s[i] * s[i]
	end

	ss3 = 0.0
	for i = eachindex(z)
		ss3 += z[i] * s[i] * s[i] * s[i]
	end

	zeta = 0.0
	for i = eachindex(z)
		zeta += z[i] * s[i] * s[i] * s[i];
	end
	zeta = 6.023e-4 * pi / 6.0 * rho * zeta

#	println("zeta = $(zeta)")

	Fnetam = 0.0;
	for i = eachindex(z)
		for j = eachindex(z)
			temp = h[i]^(1/3) + h[j]^(1/3)
			temp = temp * temp * temp / 8.0
			temp = temp^(4/3)
			Mwij = 2.0 * Mw[i] * Mw[j] / (Mw[i] + Mw[j])
			temp2 = sqrt(f[i] * f[j])
			Fnetam += z[i] * z[j] * temp * sqrt(temp2 * Mwij)
		end
	end
	Fnetam = Fnetam / sqrt(MwRef) / hm / hm

#	println("Fnetam = $(Fnetam)")

	g = zeros(nc,nc)
	temp = 1.0 - zeta
	for i = eachindex(z)
		for j = eachindex(z)
			theta = s[i] * s[j] / (s[i] + s[j]) * ss2 / ss3
			g[i,j] = 1.0 / temp + 3.0 * zeta / temp / temp * theta + 2.0 * zeta * zeta / temp / temp / temp * theta * theta
		end
	end

	Y = zeros(nc);
	for i = eachindex(z)
		temp = 0.0;
		for j = eachindex(z)
			Mwij = Mw[j] / (Mw[i] + Mw[j])
			sij = (s[i] + s[j]) / 2.0
			temp += z[j] * Mwij * sij * sij * sij * g[i,j]
		end
		Y[i] = z[i] * (1.0 + 8.0 * pi / 15.0 * 6.023e-4 * rho * temp)
	end

	neta0 = zeros(nc,nc)
	for i = eachindex(z)
		for j = eachindex(z)
			Mwij = 2.0 * Mw[i] * Mw[j] / (Mw[i] + Mw[j])
			sij = (s[i] + s[j]) / 2.0
			neta0[i,j] = 26.692 * sqrt(Mwij * T) / sij / sij
		end
	end

	B = zeros(nc,nc)
	for i = eachindex(z)
		for j = eachindex(z)
			temp = 0.0;
			deltaij = 0.0;
			if (i == j) 
				deltaij = 1.0
			end
			for k = eachindex(z)
				Mwik1 = Mw[k] / (Mw[i] + Mw[k])
				Mwik2 = Mw[i] / Mw[k]
				deltajk = 0.0
				if (j == k) 
					deltajk = 1.0
				end
				temp += z[i] * z[k] * g[i,k] / neta0[i,k] * Mwik1 * Mwik1 * ((1.0 + 5.0 / 3.0 * Mwik2) * deltaij - 2.0 / 3.0 * Mwik2 * deltajk)
			end
			B[i,j] = 2.0 * temp
		end
	end

	beta = zeros(nc)
	beta = B \ Y

#	println("beta = $(beta)")

	temp = 2.0 * pi / 3.0 * 6.023e-4
	temp *= temp
	temp *= 48.0 / 25.0 / pi * rho * rho

	neta_Enskog_m = 0.0;
	for i = eachindex(z)
		neta_Enskog_m += Y[i] * beta[i]
	end

#	println("neta_Enskog_m = $(neta_Enskog_m)")

	temp2 = 0.0
	for i = eachindex(z)
		for j = eachindex(z)
			sij = (s[i] + s[j]) / 2.0
			temp2 += z[i] * z[j] * sij * sij * sij * sij * sij * sij * neta0[i,j] * g[i,j]
		end
	end

	neta_Enskog_m += temp * temp2

	sx = 0.0
	for i = eachindex(z)
		for j = eachindex(z)
			sij = (s[i] + s[j]) / 2.0
			sx += z[i] * z[j] * sij * sij * sij
		end
	end
	sx = sx^(1/3)

	Mwx = 0.0
	for i = eachindex(z)
		for j = eachindex(z)
			Mwij = 2.0 * Mw[i] * Mw[j] / (Mw[i] + Mw[j])
			sij = (s[i] + s[j]) / 2.0
			Mwx += z[i] * z[j] * sqrt(Mwij) * sij * sij * sij * sij
		end
	end
	Mwx = Mwx / sx / sx / sx / sx
	Mwx = Mwx * Mwx

	zetax = 6.023e-4 * pi / 6.0 * rho * sx * sx * sx

	temp = (1.0 - zetax)
	gxx = 1.0 / temp + 3.0 * zetax / temp / temp * 0.5 + 2.0 * zetax * zetax / temp / temp / temp * 0.25

	Yx = 1.0 + 8.0 * pi / 15.0 * 6.023e-4 * rho * sx * sx * sx * gxx / 2.0

	neta0x = 26.692 * sqrt(Mwx * T) / sx / sx

	Bxx = gxx / neta0x

	betax = Yx / Bxx

	temp = 2.0 * pi / 3.0 * 6.023e-4
	temp *= temp
	temp *= 48.0 / 25.0 / pi * rho * rho

	neta_Enskog_x = betax * Yx + temp * sx * sx * sx * sx * sx * sx * neta0x * gxx

#	println("neta_Enskog_x = $(neta_Enskog_x)")

	Dneta_Enskog = 0.1 * (neta_Enskog_m - neta_Enskog_x)

#	println("Dneta_Enskog = $(Dneta_Enskog)")

	for i = eachindex(z)
		s[i] = 0.809 * Vc[i]^(1/3)
	end

	sm = 0.0
	for i = eachindex(z)
		for j = eachindex(z)
			sij = sqrt(s[i] * s[j])
			sm += z[i] * z[j] * sij * sij * sij
		end
	end
	sm = sm^(1/3)

#	println("sm = $(sm)")

	eta = zeros(nc)
	for i = eachindex(z)
		eta[i] = Tc[i] / 1.2593
	end

	etam = 0.0
	for i = eachindex(z)
		for j = eachindex(z)
			sij = sqrt(s[i] * s[j])
			etaij = sqrt(eta[i] * eta[j])
			etam += z[i] * z[j] * etaij * sij * sij * sij
		end
	end
	etam = etam / sm / sm / sm

#	println("etam = $(etam)")

	Mwm = 0.0
	for i = eachindex(z)
		for j = eachindex(z)
			sij = sqrt(s[i] * s[j])
			etaij = sqrt(eta[i] * eta[j])
			Mwij = 2.0 * Mw[i] * Mw[j] / (Mw[i] + Mw[j])
			Mwm += z[i] * z[j] * etaij * sij * sij * sqrt(Mwij)
		end
	end
	Mwm = Mwm / etam / sm / sm
	Mwm = Mwm * Mwm

	Omegam = 0.0
	for i = eachindex(z)
		for j = eachindex(z)
			sij = sqrt(s[i] * s[j])
			Omegaij = (Omega[i] + Omega[j]) / 2.0
			Omegam += z[i] * z[j] * Omegaij * sij * sij * sij
		end
	end
	Omegam = Omegam / sm / sm / sm

	Fc = 1.0 - 0.2756 * Omegam

#	println("Fc = $(Fc)")

	FF = 1.16145 * (T / etam)^(-0.14874) + 0.52487 * exp(-0.77320 * T / etam) + 2.16178 * exp(-2.43787 * T / etam)

#	println("FF = $(FF)")

#	println("Mwm = $(Mwm)")

	neta0m = 2.6692 * Fc * sqrt(Mwm * T) / sm / sm / FF

#	println("neta0m = $(neta0m)")

	neta = (neta0m + Fnetam * DnetaR + Dneta_Enskog) * 1.0e-6

	return neta # Pa.s

end


function SuperTRAPP_lambda(pures, model_params, modelC3H6_params, Tin, rhoin, zin)

	R̄ = Clapeyron.R̄

#	println("R̄ = $(R̄)")

	C = [
	15.2583985944,
	5.29917319127,
	-3.05330414748,
	0.450477583739,
	1.03144050679,
	-0.185480417707]

	T = Tin
	rho = rhoin/1000.0 # kmol/m3
	z = zin

	nc = length(z)

#    pures = Clapeyron.split_pure_model(model)
	
    #=
	modelC3H6 = GERG2008(["propane"])

	MwRef = modelC3H6.params.Mw.values[1]
	TcRef = modelC3H6.params.Tc.values[1]
	PcRef = modelC3H6.params.Pc.values[1]
	VcRef = modelC3H6.params.Vc.values[1]*100.0^3 # cm3/mol
	rhocRef = 1000.0/VcRef
	OmegaRef = modelC3H6.params.acentricfactor.values[1]
	ZcRef = PcRef*VcRef/(100.0^3)/R̄/TcRef

	println("Mw, Tc, Pc, Vc, Omega = $(MwRef), $(TcRef), $(PcRef), $(VcRef), $(OmegaRef), $(ZcRef)")
    =#

	MwRef = modelC3H6_params["Mw"].values[1]
	TcRef = modelC3H6_params["Tc"].values[1]
	PcRef = modelC3H6_params["Pc"].values[1]
	VcRef = modelC3H6_params["Vc"].values[1]*100.0^3 # cm3/mol
	rhocRef = 1000.0/VcRef
	OmegaRef = modelC3H6_params["acentricfactor"].values[1]
	ZcRef = PcRef*VcRef/(100.0^3)/R̄/TcRef

#	println("Mw, Tc, Pc, Vc, Omega = $(MwRef), $(TcRef), $(PcRef), $(VcRef), $(OmegaRef), $(ZcRef)")

    #=
	Mw = model.params.Mw.values
	Tc = model.params.Tc.values
	Pc = model.params.Pc.values
	Vc = model.params.Vc.values*100.0^3 # cm3/mol
	Zc = zeros(nc)
	for i = eachindex(z)
		Zc[i] = Pc[i]*Vc[i]/(100.0^3)/R̄/Tc[i]
    end
	Omega = model.params.acentricfactor.values

	println("Mw, Tc, Pc, Vc, Omega = $(Mw), $(Tc), $(Pc), $(Vc), $(Omega), $(Zc)")
    =#

    Mw = model_params["Mw"].values
	Tc = model_params["Tc"].values
	Pc = model_params["Pc"].values
	Vc = model_params["Vc"].values*100.0^3 # cm3/mol
	Zc = zeros(nc)
	for i = eachindex(z)
		Zc[i] = Pc[i]*Vc[i]/(100.0^3)/R̄/Tc[i]
    end
	Omega = model_params["acentricfactor"].values

#	println("Mw, Tc, Pc, Vc, Omega = $(Mw), $(Tc), $(Pc), $(Vc), $(Omega), $(Zc)")

	Mw_mix = 0.0
	for i = eachindex(z)
		Mw_mix += z[i] * Mw[i]
    end

	f = zeros(nc)
	for i = eachindex(z)
		f[i] = Tc[i] / TcRef * (1.0 + (Omega[i] - OmegaRef) * (0.05202976 - 0.7498189 * log(T / Tc[i])))
	end

#	println("f = $(f)")

	h = zeros(nc)
	for i = eachindex(z)
		h[i] = Vc[i] / VcRef * ZcRef / Zc[i] * (1.0 + (Omega[i] - OmegaRef) * (0.1435971 - 0.2821562 * log(T / Tc[i])))
	end

#	println("h = $(h)")

	hm = 0.0;
	for i = eachindex(z)
		for j = eachindex(z)
			temp = h[i]^(1/3) + h[j]^(1/3)
			temp = temp * temp * temp / 8.0
			hm += z[i] * z[j] * temp
		end
	end

#	println("hm = $(hm)")

	fm = 0.0
	for i = eachindex(z)
		for j = eachindex(z)
			temp = h[i]^(1/3) + h[j]^(1/3)
			temp = temp * temp * temp / 8.0
			fm += z[i] * z[j] * temp * sqrt(f[i] * f[j])
		end
	end
	fm = fm / hm

#	println("fm = $(fm)")

	To = T / fm;
	rhoo = rho * hm;

	Flambdam = 0.0;
	for i = eachindex(z)
		for j = eachindex(z)
			temp = h[i]^(1/3) + h[j]^(1/3)
			temp = temp * temp * temp / 8.0
			temp = temp^(-4/3)
			invMwij = (1.0 / Mw[i] + 1.0 / Mw[j]) / 2.0
			temp2 = sqrt(f[i] * f[j])
			Flambdam += z[i] * z[j] * temp * sqrt(temp2 * invMwij)
		end
	end
	Flambdam = Flambdam * sqrt(MwRef) * hm^(2.0 / 3.0)

	Omegam = 0.0;
	for i = eachindex(z)
		Omegam += z[i] * Omega[i]
	end

	Xm = sqrt((1.0 + 2.186634 * (Omegam - OmegaRef))) / (1.0 - 0.5050059 * (Omegam - OmegaRef))

#	println("Xm = $(Xm)")

	Tr = To / TcRef;
	rhor = rhoo / rhocRef;

	DlambdaR = C[1] * rhor + C[2] * rhor * rhor * rhor + (C[3] + C[4] / Tr) * rhor * rhor * rhor * rhor + (C[5] + C[6] / Tr) * rhor * rhor * rhor * rhor * rhor;

#	println("DlambdaR = $(DlambdaR)")

	s = zeros(nc)
	for i = eachindex(z)
		s[i] = 0.809 * Vc[i]^(1/3)
	end

	eta = zeros(nc)
	for i = eachindex(z)
		eta[i] = Tc[i] / 1.2593
	end

	Cp = zeros(nc)
	for i = eachindex(z)
	#	y = T / (T + params[i].CpIdeal[2]);
	#	F = params[i].CpIdeal[3] + params[i].CpIdeal[4] * y + params[i].CpIdeal[5] * y * y + params[i].CpIdeal[6] * y * y * y;
	#	Cp[i] = (params[i].CpIdeal[0] + y * y * (params[i].CpIdeal[1] - params[i].CpIdeal[0]) * (1.0 + (y - 1.0) * F)) * RGas / 100.0;
		Cp[i] = isobaric_heat_capacity.(pures[i],1e2,T)
	end

#	println("Cp = $(Cp)")

	neta0x = zeros(nc)
	for i = eachindex(z)
		Fc = 1.0 - 0.2756 * Omega[i];
		FF = 1.16145 * (T / eta[i])^(-0.14874) + 0.52487 * exp(-0.77320 * T / eta[i]) + 2.16178 * exp(-2.43787 * T / eta[i])
		neta0x[i] = 2.6692 * Fc * sqrt(Mw[i] * T) / s[i] / s[i] / FF
	end

#	println("neta0x = $(neta0x)")

	lambda0x = zeros(nc)
	for i = eachindex(z)
		lambda0x[i] = 15.0 / 4.0 * R̄ * neta0x[i] * 1.0e-6 / Mw[i] * 1000.0
	end

#	println("lambda0x = $(lambda0x)")

	lambdaintx = zeros(nc)
	for i = eachindex(z)
		lambdaintx[i] = 1.32 * neta0x[i] * 1.0e-6 / Mw[i] * (Cp[i] - 2.5 * R̄) * 1000.0
	end

#	println("lambdaintx = $(lambdaintx)")

	lambda0m = 0.0
	for i = eachindex(z)
		for j = eachindex(z)
			lambda0ij = 2.0 * lambda0x[i] * lambda0x[j] / (lambda0x[i] + lambda0x[j])
			lambda0m += z[i] * z[j] * lambda0ij
		end
	end

#	println("lambda0m = $(lambda0m)")

	lambdaintm = 0.0
	for i = eachindex(z)
		for j = eachindex(z)
			lambdaintij = 2.0 * lambdaintx[i] * lambdaintx[j] / (lambdaintx[i] + lambdaintx[j]);
			lambdaintm += z[i] * z[j] * lambdaintij;
		end
	end

#	println("lambdaintm = $(lambdaintm)")

	lambda = lambda0m + lambdaintm + Flambdam * Xm * DlambdaR / 1000.0

	return lambda # W/m/K

end

println("Calculate Triple Point for CO2 based on SAFTVRMie:")

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
pbase = 27.0e5
#pbase = 5.25e5
psat = pbase
sat = saturation_temperature(modelCO2, psat)
Tsat = sat[1]
vl = sat[2]
vv = sat[3]

println("T sat = $(Tsat - 273.15)")

zsCO2 = 0.997275
zsMeOH = 0.0001
x = [zsMeOH,(1.0-zsCO2)*(1.0 - zsMeOH),zsCO2*(1.0 - zsMeOH)]

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
#    T_factor = 0.037 # 27 bara
    T_factor = 0.24 # 27 bara
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
for ip = eachindex(βf)
    println("Volume Fraction $(ip) = $(βf[ip]*vf[ip]/vt)")
end
#println("Volume Fraction Vapour = $(βf[1]*vf[1]/vt)")
#println("Volume Fraction Liquid = $(βf[2]*vf[2]/vt)")

for ip = eachindex(βf)
    phase_type = Clapeyron.VT_identify_phase(model, vf[ip], Tin, zin)
    println("phase type = $(phase_type)")
    if phase_type == :vapour
        println("phase type = vapour") 
    elseif phase_type == :liquid
        println("phase type = liquid")
    else
        println("phase type = unknown")
    end
end

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
    abs_tol::T = rel_tol*rel_tol,
    iterations::Integer = 250,
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

function cp_dryice(p,T,Tref,href,sref,gint,sint,fach,facs)
    g(x) = gibbs_dryice(p,x,Tref,href,sref,gint,sint,fach,facs)
    dg(x) = Clapeyron.Solvers.derivative(g,x)
    ddg(x) = Clapeyron.Solvers.derivative(dg,x)
    return -T*ddg(T)
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

cptp = cp_dryice(ptp,Ttp,Ttp,hltp,sltp,gtp0,stp0,1,1)
println("cp at tp = $(cptp) J/mol/K")

println("Check complete")
println("")

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
    gs = gibbs_dryice(p,T,Ttp,hltp,sltp,gtp0,stp0,1,1)
    gm = xs*gs + (1-xs)*gt
    return gm/Clapeyron.R̄/T
end

function pT_flash_CO2Solid(model,p,T,z,verbose::Bool)
    iCO2 = 0
    for i = eachindex(model.components)
        if model.components[i] == "carbon dioxide"
            iCO2 = i
        end
    end
    @assert iCO2 > 0 "carbon dioxide must me present in components"
    xs_eps = eps(Float64)
    xs_min = 0.0
    xs_max = z[iCO2]*(1.0 - xs_eps)
    f = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(verbose = verbose))
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
    μs = gibbs_dryice(p,T,Ttp,hltp,sltp,gtp0,stp0,1,1)
    fnc(xs) = gibbsmix(model,p,T,z,xs,iCO2)
    if μs < μf[iCO2]
        res = brentmin(fnc, xs_min, xs_max, false)
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
        return xs,zf,f,fnc(xs)
    else
        return xs_min,z,f,fnc(xs_min)
    end
end

println("Check dry ice formation")
println("")

Tsolid = -78.85+273.15
psolid = 1.01325e5

zsCO2 = 0.997275
zsMeOH = 0.475
x = [zsMeOH,(1.0-zsCO2)*(1.0 - zsMeOH),zsCO2*(1.0 - zsMeOH)]

βCO2, z_fluid, flash_fluid, gibbs_min = pT_flash_CO2Solid(model,psolid,Tsolid,x,false)

println("moles of solid CO2 = $(βCO2)")
println("moles of fluid = $(1 - βCO2)")

println("fluid βf = $(flash_fluid.fractions)")
println("fluid vf = $(flash_fluid.volumes) m3/mol")
println("z_fluid = $(z_fluid)")
for i = eachindex(flash_fluid.fractions)
    println("fluid xf[$(i)] = $(flash_fluid.compositions[i])")
end
println("gibbs min = $(gibbs_min)")

println("")
println("Check dry ice formation complete")

#=
println("")
println("Check dry ice flash")

psf = 5.0e5

pbubble = psf
mix_Tbub = bubble_temperature(model,pbubble,x)

println("mix_Tbub = $(mix_Tbub)")
println("T bubble = $(mix_Tbub[1]-273.15)")

s_bubble = Clapeyron.VT_entropy(model,mix_Tbub[2],mix_Tbub[1],x)

println("s bubble = $(s_bubble)")

μf_bub = Clapeyron.VT_chemical_potential(model, mix_Tbub[2], mix_Tbub[1], x)

println("μf_bub = $(μf_bub)")

pdew = psf
mix_Tdew = dew_temperature(model,pdew,x)

println("mix_Tdew = $(mix_Tdew)")
println("T dew = $(mix_Tdew[1]-273.15)")

s_dew = Clapeyron.VT_entropy(model,mix_Tdew[3],mix_Tdew[1],x)

println("s dew = $(s_dew)")

μf_dew = Clapeyron.VT_chemical_potential(model, mix_Tdew[3], mix_Tdew[1], x)
μs = gibbs_dryice(psf,mix_Tdew[1],Ttp,hltp,sltp,gtp0,stp0,1,1)

println("μf_dew = $(μf_dew)")
println("μs = $(μs)")

Tsf = mix_Tdew[1]+0.312

βCO2, z_fluid, flash_fluid, gibbs_min = pT_flash_CO2Solid(model,psf,Tsf,x,false)

println("moles of solid CO2 = $(βCO2)")
println("moles of fluid = $(1 - βCO2)")

println("fluid βf = $(flash_fluid.fractions)")
println("fluid vf = $(flash_fluid.volumes) m3/mol")
println("z_fluid = $(z_fluid)")
for i = eachindex(flash_fluid.fractions)
    println("fluid xf[$(i)] = $(flash_fluid.compositions[i])")
end
println("gibbs min = $(gibbs_min)")

μf_dew = Clapeyron.VT_chemical_potential(model, flash_fluid.volumes[1], Tsf, z_fluid)
μs = gibbs_dryice(psf,Tsf,Ttp,hltp,sltp,gtp0,stp0,1,1)

println("μf_dew = $(μf_dew)")
println("μs = $(μs)")

sf = Clapeyron.VT_entropy(model,flash_fluid.volumes[1], Tsf,z_fluid)
ss = entropy_dryice(psf,Tsf,Ttp,hltp,sltp,gtp0,stp0,1,1)

st = βCO2*ss + (1.0 - βCO2)*sf
println("st = $(st)")

s_spec = -40.0
ps = psf
Tsf = mix_Tdew[1]
zs = x

# functions for fluid pressure entropy flash

function Qs(model,ps,Ts,zs,s_spec)
    β, zf, fr, gs = pT_flash_CO2Solid(model,ps,Ts,zs,false)
    return (gs*Ts + s_spec*Ts/Clapeyron.R̄)
end

function psflashmin(T)
   Q = Qs(model,ps,T,zs,s_spec)
   return -Q
end

function Qs2(model,ps,Ts,zs,s_spec)
    β, zf, fr, gs = pT_flash_CO2Solid(model,ps,Ts,zs,false)
    sf = Clapeyron.VT_entropy(model,fr.volumes[1],Ts,zf)
    ss = entropy_dryice(ps,Ts,Ttp,hltp,sltp,gtp0,stp0,1,1)
    st = β*ss + (1.0 - β)*sf
    return (gs*Ts + s_spec*Ts/Clapeyron.R̄),st
end

function psflashmin2(T)
   Q,st = Qs2(model,ps,T,zs,s_spec)
   return -Q,st
end

Ta = Tsf + 0.1
Tc = Tsf - 0.1

Ta, Tc, Qa, Qc = bracketmin(psflashmin,Ta,Tc,5.0,true)

bmres = brentmin(psflashmin, Ta, Tc, true)

Tsf = bmres[1]
println("Tsf = $(Tsf - 273.15)")

βCO2, z_fluid, flash_fluid, gibbs_min = pT_flash_CO2Solid(model,psf,Tsf,x,false)

println("moles of solid CO2 = $(βCO2)")
println("moles of fluid = $(1 - βCO2)")

println("fluid βf = $(flash_fluid.fractions)")
println("fluid vf = $(flash_fluid.volumes) m3/mol")
println("z_fluid = $(z_fluid)")
for i = eachindex(flash_fluid.fractions)
    println("fluid xf[$(i)] = $(flash_fluid.compositions[i])")
end
println("gibbs min = $(gibbs_min)")

μf_dew = Clapeyron.VT_chemical_potential(model, flash_fluid.volumes[1], Tsf, z_fluid)
μs = gibbs_dryice(psf,Tsf,Ttp,hltp,sltp,gtp0,stp0,1,1)

println("μf_dew = $(μf_dew)")
println("μs = $(μs)")

sf = Clapeyron.VT_entropy(model,flash_fluid.volumes[1], Tsf,z_fluid)
ss = entropy_dryice(psf,Tsf,Ttp,hltp,sltp,gtp0,stp0,1,1)

st = βCO2*ss + (1.0 - βCO2)*sf
println("st = $(st)")

N  = 10

Tplt = LinRange(mix_Tdew[1]-273.15+0.3, mix_Tdew[1]-273.15+0.32,N)
Qplt = zeros(N)
splt = zeros(N)

for i = eachindex(Qplt)
    Qplt[i],splt[i] = psflashmin2(Tplt[i]+273.15)
    println("Qplt[$(i)],splt[$(i)] = $(Qplt[i]),$(splt[i])")
end

p1 = plot(
    label=["Q"], 
    [Tplt,],
    [Qplt,],
    xlabel = "Temperaure [deg C]",
    ylabel = "Q",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

    display(p1)

p2 = plot(
    label=["s"], 
    [Tplt],
    [splt],
    xlabel = "Temperaure [deg C]",
    ylabel = "s",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

display(p2)


println("")
println("Check dry ice flash complete")

=#

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
    phases::Vector{phase}
    phase_solid::phase
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

# Start Blowdown

ps = pbase
Ts = Tf
Twall_initial = Tf
#Twall_initial = -56.19+273.15

u_ambient = 0.1
p_ambient = 0.995e5
T_ambient = Ts
#T_ambient = -9.54+273.15
flashAir = Clapeyron.tp_flash_impl(modelAir,p_ambient,T_ambient,zAir, HELDTPFlash(verbose = false))
βf,mwf,mwtf,vf,vtf,hf,htf,sf,stf,xf = get_props(modelAir,flashAir)

rho_ambient = 1.0/vtf
mu_ambient = SuperTRAPP_mu(modelAir_params, modelC3H6_params, T_ambient, rho_ambient, zAir)
lambda_ambient = SuperTRAPP_lambda(pures_modelAir, modelAir_params, modelC3H6_params, T_ambient, rho_ambient, zAir)
cp_ambient = Clapeyron.VT_isobaric_heat_capacity(modelAir,vtf,T_ambient,zAir)/mwAir # J/Kg/K

T_ambient -= 273.15
println("T_ambient = $(T_ambient)")
println("mu_ambient = $(mu_ambient)")
println("lambda_ambient = $(lambda_ambient)")
println("cp_ambient = $(cp_ambient)")

zsCO2 = 0.997275
zsMeOH = 0.0001
zs = [zsMeOH,(1.0-zsCO2)*(1.0 - zsMeOH),zsCO2*(1.0 - zsMeOH)]

tank_diamter = 5.2 #m
tank_radius = tank_diamter/2.0
tank_xsarea = pi*tank_radius^2
tank_length_TT = 26.882 #m
# spherical ends
tank_volume = tank_xsarea*tank_length_TT + 4/3*pi*tank_radius^3
tank_length = tank_volume/tank_xsarea

println("tank volume = $(tank_volume)")
println("tank length = $(tank_length)")

level = zeros(2)
level[1] = 0.95
#level[1] = 1.0 - 0.05/100
level[2] = 1.0 - level[1]

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

function psflashmin(model,p,T,z,s, verbose)
	fr = Clapeyron.tp_flash_impl(model,p,T,z,HELDTPFlash(verbose = verbose))
	g = fr.data.g
	Q =  -(g*T + s*T/Clapeyron.R̄)
	return Q
end

function flowcalc(model,p,T,z,pin,pout,mw,h,s,rho,Cd,d0,d1,verbose)
	Ta = T - 0.1
	Tb = T + 0.1
	psfunc(x) = psflashmin(model,p,x,z,s, verbose)
	Ta,Tc,Qa,Qc = bracketmin(psfunc,Ta, Tb, 5.0, verbose)
	T, Qs = brentmin(psfunc,Ta,Tc,verbose)
	S0 = pi/4*d0^2
    S1 = pi/4*d1^2
	flash = Clapeyron.tp_flash_impl(model,p,T,z,HELDTPFlash(verbose = verbose))
	βss,mwss,mwtss,vss,vtss,hss,htss,sss,stss,xss  = get_props(model,flash)
	rhoss = 1.0/vtss
	Δh = h - htss
	flow = Cd*S0*rhoss*sign(pin - pout)*sqrt(2.0*abs(Δh/mw/(1.0 - (rhoss/rho*S0/S1)^2)))
	Q = -flow
	return Q
end

function OrificeFlow(model,pin,Tin,zin,pout,pcrit_ratio,Tcrit_ratio,Cd,d0,d1,verbose)
	flashin = Clapeyron.tp_flash_impl(model,pin,Tin,zin,HELDTPFlash(verbose = false))
	βin,mwin,mwtin,vf,vtin,hin,htin,sin,stin,xin = get_props(model,flashin)
	pcrit = pcrit_ratio*pin
	Tcrit = Tcrit_ratio*Tin
	rhotin = 1.0/vtin
	flowFunc(px) = flowcalc(model,px,Tcrit,zin,pin,pout,mwtin,htin,stin,rhotin,Cd,d0,d1,false)
	pa = pcrit - (pin - pout)/100.0
	pb = pcrit + (pin - pout)/100.0
	pa,pc,Qa,Qc = bracketmin(flowFunc,pa, pb, (pin - pout)/10.0, verbose)
	pcrit, Qcrit = brentmin(flowFunc,pa,pc,verbose, 1.0e-4, 1.0e-8)
	Ta = Tcrit - 0.1
	Tb = Tcrit + 0.1
	psfunc(Tx) = psflashmin(model,pcrit,Tx,zin,stin, false)
	Ta,Tc,Qa,Qc = bracketmin(psfunc,Ta, Tb, 5.0, verbose)
	Tcrit, Qs = brentmin(psfunc,Ta,Tc,verbose)
	return -Qcrit, pcrit, Tcrit
end

function ValveFlow(pin,pout,pcrit,rhoin,CV)
	flow = CV*rhoin*sign(pin - max(pcrit*pin,pout))*sqrt(abs(pin - max(pcrit*pin,pout))/rhoin)
	return flow
end

#=
pin = 2.675e6
Tin = 262.9517719353827
zin = [0.9956678999065446, 0.004329942967552403, 2.157125903032932e-6]
pout = p_ambient

flow,pcrit,Tcrit  = OrificeFlow(model,pin,Tin,zin,pout,0.58,0.93,0.85,2.0*25.4/1000.0/sqrt(2),8.0*25.4/1000.0,true)
println("flow  = $(flow) mol/s")
println("pcrit = $(pcrit/1e5) bara")
println("pcrit ratio = $(pcrit/pin)")
println("Tcrit = $(Tcrit - 273.15) deg C")
println("Tcrit ratio = $(Tcrit/Tin)")

=#

#=

pin = 2.65e6
Tin = 262.6682856518778
zin = [0.998864211267771, 0.0002566652842950465, 0.0008791234479337743]
pout = p_ambient

flow,pcrit,Tcrit  = OrificeFlow(model,pin,Tin,zin,pout,0.7592642284676264,0.9650535130492818,0.85,2.0*25.4/1000.0/sqrt(2),8.0*25.4/1000.0,true)
println("flow  = $(flow) mol/s")
println("pcrit = $(pcrit/1e5) bara")
println("pcrit ratio = $(pcrit/pin)")
println("Tcrit = $(Tcrit - 273.15) deg C")
println("Tcrit ratio = $(Tcrit/Tin)")

=#

function BlowDown(model, ps, Ts, zs, tank_volume, tank_radius, tank_length, level, pe, nstep, nVTK)

#    println("zs = $(zs)")
    
    mix_Tbub = bubble_temperature(model,ps,zs)
    Tbub = mix_Tbub[1]
    mix_Tdew = dew_temperature(model,ps,zs)
    Tdew = mix_Tdew[1]

    #holdups_last = Vector{holdup}(undef,0)
    holdups = Vector{holdup}(undef,0)
    tank_xsa_holdups = Vector{Float64}(undef,0)

    vv = volume(modelCO2,ptp,Ttp,[1])
    hv = Clapeyron.VT_enthalpy(modelCO2,vv,Ttp,[1])
    sv = Clapeyron.VT_entropy(modelCO2,vv,Ttp,[1])

    v_solid = volume_dryice(ptp,Ttp,Ttp,hltp,sltp,gtp0,stp0,1,1)
    s_solid = entropy_dryice(ptp,Ttp,Ttp,hltp,sltp,gtp0,stp0,1,1)
    g_solid = gibbs_dryice(ptp,Ttp,Ttp,hltp,sltp,gtp0,stp0,1,1)
    h_solid = g_solid + Ttp*s_solid

    if Ts <= Tbub

        ybub = zeros(length(zs))
        for ic = eachindex(y)
            ybub[ic] = mix_Tbub[4][ic]
        end
        mw = Clapeyron.molecular_weight(model,ybub)
        v = mix_Tbub[3]
        h = Clapeyron.VT_enthalpy(model,v,Ts,ybub)
        s = Clapeyron.VT_entropy(model,v,Ts,ybub)
        tv = 0.0
        tm = 0.0
        phases = Vector{phase}(undef,0)
        push!(phases, phase(1.0,tv,tm,mw,v,h,s,ybub))
        phase_solid = phase(0.0,0.0,0.0,mwCO2,v_solid,h_solid,s_solid,z_solid)
        push!(holdups, holdup(ps,Ts,tv,tm,mw,v,h,s,ybub,phases,phase_solid))

        xbub = zeros(length(zs))
        for ic = eachindex(xbub)
            xbub[ic] = zs[ic]
        end
        flashx = Clapeyron.tp_flash_impl(model,ps,Ts,xbub, HELDTPFlash(verbose = false))
        βf,mwf,mwtf,vf,vtf,hf,htf,sf,stf,xf = get_props(model,flashx)
        tv = tank_volume
        tm = tv/v
        phases = Vector{phase}(undef,0)
        push!(phases, phase(1.0,tv,tm,mwtf,vtf,htf,stf,xbub))
        phase_solid = phase(0.0,0.0,0.0,mwCO2,v_solid,h_solid,s_solid,z_solid)
        push!(holdups, holdup(ps,Ts,tv,tm,mwtf,vtf,htf,stf,xbub,phases,phase_solid))

        level[1] = 0.0
        level[2] = 1.0

        push!(tank_xsa_holdups, level[1]*tank_volume/tank_length)
        push!(tank_xsa_holdups, level[2]*tank_volume/tank_length)

        println("Bubble point holdups = $(holdups)")

    elseif Ts >= Tdew

        ydew = zeros(length(zs))
        for ic = eachindex(ydew)
            ydew[ic] = zs[ic]
        end
        flashy = Clapeyron.tp_flash_impl(model,ps,Ts,ydew, HELDTPFlash(verbose = false))
        βf,mwf,mwtf,vf,vtf,hf,htf,sf,stf,xf = get_props(model,flashy)
        mw = Clapeyron.molecular_weight(model,ydew)
        tv = tank_volume
        tm = tv/v
        phases = Vector{phase}(undef,0)
        push!(phases, phase(1.0,tv,tm,mwtf,vtf,htf,stf,ydew))
        phase_solid = phase(0.0,0.0,0.0,mwCO2,v_solid,h_solid,s_solid,z_solid)
        push!(holdups, holdup(ps,Ts,tv,tm,mwtf,vtf,htf,stf,ydew,phases,phase_solid))

        xdew = zeros(length(zs))
        for ic = eachindex(x)
            xdew[ic] = mix_Tdew[4][ic]
        end
        mw = Clapeyron.molecular_weight(model,xdew)
        v = mix_Tdew[2]
        h = Clapeyron.VT_enthalpy(model,v,Ts,xdew)
        s = Clapeyron.VT_entropy(model,v,Ts,xdew)
        tv = 0.0
        tm = 0.0
        phases = Vector{phase}(undef,0)
        push!(phases, phase(1.0,tv,tm,mw,v,h,s,xdew))
        phase_solid = phase(0.0,0.0,0.0,mwCO2,v_solid,h_solid,s_solid,z_solid)
        push!(holdups, holdup(ps,Tdew,tv,tm,mw,v,h,s,xdew,phases,phase_solid))

        level[1] = 1.0
        level[2] = 0.0

        push!(tank_xsa_holdups, level[1]*tank_volume/tank_length)
        push!(tank_xsa_holdups, level[2]*tank_volume/tank_length)

        println("Dew point holdups = $(holdups)")

    else

        flash = Clapeyron.tp_flash_impl(model,ps,Ts,zs, HELDTPFlash(verbose = false))
        βf,mwf,mwtf,vf,vtf,hf,htf,sf,stf,xf = get_props(model,flash) 
        for ip = eachindex(βf)
            phases = Vector{phase}(undef,0)
            push!(phases, phase(1.0,level[ip]*tank_volume,level[ip]*tank_volume/vf[ip],mwf[ip],vf[ip],hf[ip],sf[ip],xf[ip]))
            phase_solid = phase(0.0,0.0,0.0,mwCO2,v_solid,h_solid,s_solid,z_solid)
            push!(holdups, holdup(ps,Ts,level[ip]*tank_volume,level[ip]*tank_volume/vf[ip],mwf[ip],vf[ip],hf[ip],sf[ip],xf[ip],phases,phase_solid))
            push!(tank_xsa_holdups, level[ip]*tank_volume/tank_length)
        end

        println("Two phase holdups = $(holdups)")

    end

    phases = Vector{phase}(undef,0)
    push!(phases, phase(0.0,0.0,0.0,mwCO2,vv,hv,sv,z_solid))
    phase_solid = phase(1.0,0.0,0.0,mwCO2,v_solid,h_solid,s_solid,z_solid)
    push!(holdups, holdup(ptp,Ttp,0.0,0.0,mwCO2,v_solid,h_solid,h_solid,z_solid,phases,phase_solid))
    push!(tank_xsa_holdups, 0.0)

    function phFlashFluids!(holdups,ΔT,tiny_holdup,verbose)
        
        for ih = 1:2

            if holdups[ih].moles > tiny_holdup

                hs = holdups[ih].h
                Ts = holdups[ih].T
                ps = holdups[ih].p
                zs = holdups[ih].z

                function Qh(model,p,T,x,hspec)
                    β, zf, fr, g = pT_flash_CO2Solid(model,p,T,x,false)
                    return (g - hspec/Clapeyron.R̄/T)
                end

                function phflashmin(T)
                    return -Qh(model,ps,T,zs,hs)
                end

                Ta = Ts + 2.0*ΔT[ih]
                Tc = Ts - 2.0*ΔT[ih]

                Ta, Tc, Qa, Qc = bracketmin(phflashmin,Ta,Tc,5.0,false)

                bmres = brentmin(phflashmin, Ta, Tc, false)

                Ts = bmres[1]
                holdups[ih].T = Ts

                β_solid, z_fluid, flash, g_total = pT_flash_CO2Solid(model,holdups[ih].p,holdups[ih].T,holdups[ih].z,verbose)
                v_solid = volume_dryice(ps,Ts,Ttp,hltp,sltp,gtp0,stp0,1,1)
                s_solid = entropy_dryice(ps,Ts,Ttp,hltp,sltp,gtp0,stp0,1,1)
                g_solid = gibbs_dryice(ps,Ts,Ts,hltp,sltp,gtp0,stp0,1,1)
                h_solid = g_solid + Ts*s_solid
                # if we form a solid then the composition of the fluid changes.
                holdups[ih].z = z_fluid
                βs,mws,mwts,vs,vts,hs,hts,ss,sts,xs  = get_props(model,flash)

                holdups[ih].volume  = (1.0 - β_solid)*holdups[ih].moles*vts + β_solid*holdups[ih].moles*v_solid
                holdups[ih].mw = (1.0 - β_solid)*mwts + β_solid*mwCO2
                holdups[ih].v  = (1.0 - β_solid)*vts + β_solid*v_solid
                holdups[ih].h  = (1.0 - β_solid)*hts + β_solid*h_solid
                holdups[ih].s  = (1.0 - β_solid)*sts + β_solid*s_solid

                moles_fluid = max((1.0 - β_solid)*holdups[ih].moles,0.0)
                phases = Vector{phase}(undef,0)
                for ip = eachindex(βs)      
                    push!(phases, phase(βs[ip]*(1.0 - β_solid),βs[ip]*moles_fluid*vs[ip],βs[ip]*moles_fluid,mws[ip],vs[ip],hs[ip],ss[ip],xs[ip]))
                end
                holdups[ih].phases = phases

                moles_solid = max(β_solid*holdups[ih].moles,0.0)
                holdups[ih].phase_solid =  phase(β_solid,moles_solid*v_solid,moles_solid,mwCO2,v_solid,h_solid,s_solid,z_solid)

                if verbose
                    if β_solid > 0.0
                        println("dry ice detected in phase[$(ih)] = $(β_solid)")
                    else
                        println("no dry ice detected in phase[$(ih)]")
                    end
                end

            else

                # no moles so correct volume of phases and overall
                holdups[ih].moles  = 0.0
                holdups[ih].volume  = 0.0
                for ip = eachindex(length(holdups[ih].phases))
                    holdups[ih].phases[ip].moles = 0.0
                    holdups[ih].phases[ip].volume = 0.0
                end
                holdups[ih].phase_solid.moles = 0.0
                holdups[ih].phase_solid.volume = 0.0

            end # moles > 0.0

        end # ih loop
    end

    function phFlashSolid!(holdups,tiny_holdup,verbose)

        βsolid = 0.0
        vs = 0.0
        hs = 0.0
        ss = 0.0
        vv = 0.0
        hv = 0.0
        sv = 0.0

        Tsv = (holdups[1].T + holdups[2].T)/2.0

        if holdups[3].moles > tiny_holdup

            # This assumes the solid form below the triple point so we have solid/vapour equilibria
            # need to check vv = volume(modelCO2,ps,T,[1]) returns the vapour root.

            function fnc_gibbs(T)
                μs = gibbs_dryice(ps,T,Ttp,hltp,sltp,gtp0,stp0,1,1)
                vv = volume(modelCO2,ps,T,[1])
                μv = Clapeyron.VT_chemical_potential(modelCO2, vv, T, [1])
                hv = Clapeyron.VT_enthalpy(modelCO2,vv,T,[1])
                sv = Clapeyron.VT_entropy(modelCO2,vv,T,[1])
                μv = hv - T*sv
                return μs - μv
            end

            f = fnc_gibbs(Tsv)
            error = abs(f)
            dTsv = 0.001
            while (error > sqrt(eps(Float64)))
                df_dT = (fnc_gibbs(Tsv+dTsv) - f)/dTsv
                Tsv = Tsv - f/df_dT
                f = fnc_gibbs(Tsv)
                error = abs(f)
            end

            vv = volume(modelCO2,ps,Tsv,[1])
            hv = Clapeyron.VT_enthalpy(modelCO2,vv,Tsv,[1])
            sv = Clapeyron.VT_entropy(modelCO2,vv,Tsv,[1])
            gv = hv - Tsv*sv
            μv = Clapeyron.VT_chemical_potential(modelCO2, vv, Tsv, [1])

            gs = gibbs_dryice(ps,Tsv,Ttp,hltp,sltp,gtp0,stp0,1,1)
            ss = entropy_dryice(ps,Tsv,Ttp,hltp,sltp,gtp0,stp0,1,1)
            hs = gs + Tsv*ss
            vs = volume_dryice(ps,Tsv,Ttp,hltp,sltp,gtp0,stp0,1,1)

            if holdups[3].h > hv

                βsolid = 0.0

                function fnc_hv(T)
                    vv = volume(modelCO2,ps,T,[1])
                    hv = Clapeyron.VT_enthalpy(modelCO2, vv, T, [1])
                    return hv - holdups[3].h
                end

                f = fnc_hv(Tsv)
                error = abs(f)
                dTsv = 0.001
                while (error > sqrt(eps(Float64)))
                    df_dT = (fnfnc_hvc(Tsv+dTsv) - f)/dTsv
                    Tsv = Tsv - f/df_dT
                    f = fnc_hv(Tsv)
                    error = abs(f)
                end
                if verbose
                    println("Superheated vapour")
                end

                vv = volume(modelCO2,ps,Tsv,[1])
                hv = Clapeyron.VT_enthalpy(modelCO2, vv, Tsv, [1])
                sv = Clapeyron.VT_entropy(modelCO2,vv,Tsv,[1])

            elseif holdups[3].h < hs

                βsolid = 1.0

                function fnc_hs(T)
                    gs = gibbs_dryice(ps,T,Ttp,hltp,sltp,gtp0,stp0,1,1)
                    ss = entropy_dryice(ps,T,Ttp,hltp,sltp,gtp0,stp0,1,1)
                    hs = gs + T*ss
                    return hs - holdups[3].h
                end

                f = fnc_hs(Tsv)
                error = abs(f)
                dTsv = 0.001
                while (error > sqrt(eps(Float64)))
                    df_dT = (fnc_hs(Tsv+dTsv) - f)/dTsv
                    Tsv = Tsv - f/df_dT
                    f = fnc_hs(Tsv)
                    error = abs(f)
                end
                if verbose
                    println("Subcooled solid")
                end

                gs = gibbs_dryice(ps,Tsv,Ttp,hltp,sltp,gtp0,stp0,1,1)
                ss = entropy_dryice(ps,Tsv,Ttp,hltp,sltp,gtp0,stp0,1,1)
                hs = gs + Tsv*ss
                vs = volume_dryice(ps,Tsv,Ttp,hltp,sltp,gtp0,stp0,1,1)

            else

                βsolid = (holdups[3].h - hv)/(hs - hv)
                if verbose
                    println("Equilibrium solid/vapour")
                end

            end

            if verbose
                println("βsolid = $(βsolid)")
            end

            holdups[3].T = Tsv

            moles_total = holdups[3].moles

        #    if βsolid < 1.0

                moles_vapour = max((1.0 - βsolid)*moles_total,0.0)
                phases = Vector{phase}(undef,0)
                push!(phases, phase((1.0 - βsolid),moles_vapour*vv,moles_vapour,mwCO2,vv,hv,sv,z_solid))
                holdups[3].phases = phases

                moles_solid = max(βsolid*moles_total,0.0)
                phase_solid = phase(βsolid,moles_solid*vs,moles_solid,mwCO2,vs,hs,ss,z_solid)
                holdups[3].phase_solid = phase_solid

                holdups[3].volume = moles_vapour*vv + moles_solid*vs
                holdups[3].h = (1.0 - βsolid)*hv + βsolid*hs
                holdups[3].s = (1.0 - βsolid)*sv + βsolid*ss

        #    else

        #        moles_solid = max(moles_total,0.0)
        #        phase_solid = phase(1.0,moles_solid*vs,moles_solid,mwCO2,vs,hs,ss,[1.0,0.0])
        #        holdups[3].phase_solid = phase_solid

        #        holdups[3].volume = moles_solid*vs
        #        holdups[3].h = hs
        #        holdups[3].s = ss

        #    end

        end

    end

    println("holdups = $(holdups)")

    holdups_last = deepcopy(holdups)

    Δvolume = Vector{Float64}(undef,0)
    for ih = eachindex(holdups)
        push!(Δvolume,0.0)
    end

    Δmoles = Vector{Float64}(undef,0)
    for ih = eachindex(holdups)
        moles = 0.0
        push!(Δmoles,moles)
    end

    theta_solid = 0.0
    height_solid = 0.0
    solidPercent = 0.0

    # tank liquid level based on volume percentage height based on tank diamter
    height = zeros(length(holdups))
    theta = sqrt(0.9*level[2])*360/180*pi
    if tank_xsa_holdups[2] > 10.0*eps(Float64)
        error = 1.0
        while (error > 0.0001)
            f =  (tank_radius^2)*(theta - sin(theta)) - 2.0*tank_xsa_holdups[2]
            df_dtheta = (tank_radius^2)*(1.0 - cos(theta))
            theta = theta - f/df_dtheta
            error = abs(f)
            println("theta = $(theta), error = $(error)")
        end
    else
        theta = 0.0
        error = 0.0
        println("theta = $(theta), error = $(error)")
    end

    # height based on tank diamter
    tank_diameter = 2.0*tank_radius
    levelPercent = (1.0 - cos(theta/2.0))/2.0*100.0
    height[2] = tank_diameter*levelPercent/100.0
    height[1] = tank_diameter*(1.0 - levelPercent/100.0)
    println("Level Percent = $(levelPercent)")

    # set RO critical and outlet temperatures

    Tcrit_vap = holdups[1].T-20.0
    Tout_vap = holdups[1].T-30.0
    pout_vap = 0.995e5
    pcrit_vap = 0.58*ps

    Tcrit_liq = holdups[2].T-20.0
    Tout_liq = holdups[2].T-30.0
    pout_liq = 0.995e5
    pcrit_liq = 0.71*ps

    println("")
    println("Blowndown Function Start:")
    println("Initialization:")
    println("")
    
    Δpbase = (ps - pe)/nstep
    ΔT = fill(-0.16,length(holdups))
    ΔQ = zeros(length(holdups))

    time_Plot = Vector{Float64}(undef,0)
    Twallmin_Plot = Vector{Float64}(undef,0)
    Twallmin_av_Plot = Vector{Float64}(undef,0)
    Twallvap_Plot = Vector{Float64}(undef,0)
    Twallliq_Plot = Vector{Float64}(undef,0)
    Twallnoz_Plot = Vector{Float64}(undef,0)
    Twall_MMDT_Plot = Vector{Float64}(undef,0)

    T1_Plot = Vector{Float64}(undef,0)
    p1_Plot = Vector{Float64}(undef,0)
    m1_Plot = Vector{Float64}(undef,0)

    T2_Plot = Vector{Float64}(undef,0)
    p2_Plot = Vector{Float64}(undef,0)
    m2_Plot = Vector{Float64}(undef,0)

    ps_Plot = Vector{Float64}(undef,0)
    ms_Plot = Vector{Float64}(undef,0)
    ms_kg_Plot = Vector{Float64}(undef,0)
    MeOH_kg_Plot = Vector{Float64}(undef,0)

    fv_Plot = Vector{Float64}(undef,0)
    fl_Plot = Vector{Float64}(undef,0)
    MeOH_Plot = Vector{Float64}(undef,0)

    Twallmin_out = Twall_initial - 273.15

    Twallmin = Twall_initial - 273.15
    push!(Twallmin_Plot, Twallmin)
    Twallmin_av = Twall_initial - 273.15
    push!(Twallmin_av_Plot, Twallmin_av)
    Twallvap = Twall_initial - 273.15
    push!(Twallvap_Plot, Twallvap)
    Twallliq = Twall_initial - 273.15
    push!(Twallliq_Plot, Twallliq)
    Twall_noz = Twall_initial - 273.15
    push!(Twallnoz_Plot, Twall_noz)
    push!(Twall_MMDT_Plot, -46.0)
    push!(T1_Plot, holdups[1].T-273.15)
    push!(T2_Plot, holdups[2].T-273.15)

    push!(p1_Plot, holdups[1].p/1e5)
    push!(p2_Plot, holdups[2].p/1e5)
    push!(ps_Plot, holdups[3].p/1e5)

    push!(m1_Plot, holdups[1].moles*holdups[1].mw/1000)
    push!(m2_Plot, holdups[2].moles*holdups[2].mw/1000)
    push!(ms_Plot, holdups[3].moles*holdups[3].mw/1000)
    push!(ms_kg_Plot, holdups[3].moles*holdups[3].mw)

    push!(fv_Plot, 0.0)
    push!(fl_Plot, 0.0)
    push!(MeOH_Plot, 0.0)
    push!(MeOH_kg_Plot, 0.0)

    time = 0.0
    push!(time_Plot, time)

    # Set up FEM wall

    # heat source
    heat_source = 0.0

    #properties
    rho_steel = 7850.0
    cp_steel = 490.0
    lambda_steel = 45.0
    rho_insulation = 30.0
    cp_insulation = 1500.0
    lambda_insulation = 0.025

    function property(x, x1, a1, a2)
        if x > x1
            a = a2
        else
            a = a1
        end
        return a
    end

    function fluid_profile(y, y1, y2, a1, a2, a3)
        if y > y2
            a = a3
        elseif y >= y1 && y <= y2
            a = a2
        else
            a = a1 
        end
        return a
    end

    Twall_vap = Twall_initial - 273.15
    Twall_liq = Twall_initial - 273.15
    Twall_solid = Twall_initial - 273.15
    Twall_noz_last = Twall_initial - 273.15
    Twall_outer = Twall_initial - 273.15
    #wall_thickness = 0.043
    #insulation_thickness = 0.015
    wall_thickness = 0.040
    insulation_thickness = 0.014
    overall_thickness = wall_thickness + insulation_thickness

    function setup_grid(h = 0.05)
        ## Initialize gmsh
        Gmsh.initialize()
        gmsh.option.set_number("General.Verbosity", 2)

        ## Add the points
        o  = gmsh.model.geo.add_point(0.0, 0.0, 0.0, h)
        p1 = gmsh.model.geo.add_point(0.0, -tank_radius, 0.0, h)
        p2 = gmsh.model.geo.add_point(0.0, -tank_radius - overall_thickness, 0.0, h)
        p3 = gmsh.model.geo.add_point(tank_radius + overall_thickness, 0.0, 0.0, h)
        p4 = gmsh.model.geo.add_point(0.0, tank_radius + overall_thickness, 0.0, h)
        p5 = gmsh.model.geo.add_point(0.0, tank_radius, 0.0, h)
        p6 = gmsh.model.geo.add_point(tank_radius, 0.0, 0.0, h)

        ## Add the lines
        l1 = gmsh.model.geo.add_line(p1, p2)
        l2 = gmsh.model.geo.add_circle_arc(p2, o, p3)
        l3 = gmsh.model.geo.add_circle_arc(p3, o, p4)
        l4 = gmsh.model.geo.add_line(p4, p5)
        l5 = gmsh.model.geo.add_circle_arc(p5, o, p6)
        l6 = gmsh.model.geo.add_circle_arc(p6, o, p1)

        ## Create the closed curve loop and the surface
        loop = gmsh.model.geo.add_curve_loop([l1, l2, l3, l4, l5, l6])
        surf = gmsh.model.geo.add_plane_surface([loop])

        ## Synchronize the model
        gmsh.model.geo.synchronize()

        ## Create the physical domains
        gmsh.model.add_physical_group(1, [l1], -1, "bottom")
        gmsh.model.add_physical_group(1, [l2,l3], -1, "right")
        gmsh.model.add_physical_group(1, [l4], -1, "top")
        gmsh.model.add_physical_group(1, [l5,l6], -1, "left")
        gmsh.model.add_physical_group(2, [surf])

        ## Generate a 2D mesh
        gmsh.model.mesh.generate(2)

        ## Save the mesh, and read back in as a Ferrite Grid
        grid = mktempdir() do dir
            path = joinpath(dir, "mesh.msh")
            gmsh.write(path)
            togrid(path)
        end

        ## Finalize the Gmsh library
        Gmsh.finalize()

        return grid
    end
    #nradius = 58
    nradius = 27
    h = overall_thickness/nradius # approximate element size
    println("approximate element size h = $(h)")
    grid = setup_grid(h)

    # ### Trial and test functions
    # Again, we define the structs that are responsible for the `shape_value` and `shape_gradient` evaluation.
    ip = Lagrange{RefTriangle, 1}()
    qr = QuadratureRule{RefTriangle}(2)
    cellvalues = CellValues(qr, ip);

    qr_face = FacetQuadratureRule{RefTriangle}(2);
    facetvalues = FacetValues(qr_face, ip);

    # ### Degrees of freedom
    # After this, we can define the `DofHandler` and distribute the DOFs of the problem.
    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh);

    # ### Assembling the linear system
    # As in the heat equation example we define a `doassemble!` function that assembles the diffusion parts of the equation:
    function doassemble_Kf!(K::SparseMatrixCSC, f::Vector, cellvalues::CellValues, dh::DofHandler)

        n_basefuncs = getnbasefunctions(cellvalues)
        Ke = zeros(n_basefuncs, n_basefuncs)
        fe = zeros(n_basefuncs)

        assembler = start_assemble(K, f)

        for cell in CellIterator(dh)

            fill!(Ke, 0)
            fill!(fe, 0)
            coords = getcoordinates(cell)

            reinit!(cellvalues, cell)

            for q_point in 1:getnquadpoints(cellvalues)
                dΩ = getdetJdV(cellvalues, q_point)
                coords_qp = spatial_coordinate(cellvalues, q_point, coords)
                x = coords_qp[1]
                y = coords_qp[2]
                r = sqrt(x^2 + y^2)
                for i in 1:n_basefuncs
                    v = shape_value(cellvalues, q_point, i)
                    ∇v = shape_gradient(cellvalues, q_point, i)
                    fe[i] += heat_source * v * dΩ
                    for j in 1:n_basefuncs
                        ∇u = shape_gradient(cellvalues, q_point, j)
                        Ke[i, j] += property(r,tank_radius+wall_thickness,lambda_steel,lambda_insulation) * (∇v ⋅ ∇u) * dΩ
                    end
                end
            end

            assemble!(assembler, celldofs(cell), Ke, fe)
        end

        return K, f

    end

    function doassemble_Kf_RobinBC!(K::SparseMatrixCSC, f::Vector, cellvalues::CellValues, dh::DofHandler,SL,LL,htcf1,Tf1,htcf2,Tf2,htcf3,Tf3,htcf4,Tf4)

        n_basefuncs = getnbasefunctions(cellvalues)
        Ke = zeros(n_basefuncs, n_basefuncs)
        fe = zeros(n_basefuncs)

        # Loop for the Robin boundary for inner wall
        for (cellid, faceid) in getfacetset(grid, "left")
            # Reset the buffers
            fill!(Ke, 0)
            fill!(fe, 0)

            coords = getcoordinates(grid, cellid)

            # Update FE basis
            reinit!(facetvalues, getcoordinates(grid, cellid), faceid)
            # Compute the local contribution
            for qp in 1:getnquadpoints(facetvalues)
                dΓ = getdetJdV(facetvalues, qp)
                coords_qp = spatial_coordinate(facetvalues, qp, coords)
                x = coords_qp[1]
                y = coords_qp[2]
                r = sqrt(x^2 + y^2)
                htcf = fluid_profile(y,SL,LL,htcf1,htcf2,htcf3)
                Tf = fluid_profile(y,SL,LL,Tf1,Tf2,Tf3)
                for i in 1:getnbasefunctions(facetvalues)
                    ϕᵢ = shape_value(facetvalues, qp, i)
                    fe[i] += ( ϕᵢ * htcf * Tf ) * dΓ
                    for j in 1:getnbasefunctions(facetvalues)
                        ϕⱼ = shape_value(facetvalues, qp, j)
                        Ke[i, j] += ( ϕᵢ * htcf * ϕⱼ ) * dΓ
                    end
                end
            end

            # Assemble local contribution
            dofs = celldofs(dh, cellid)
            for (i, I) in pairs(dofs)
                f[I] += fe[i]
                for (j, J) in pairs(dofs)
                    K[I, J] += Ke[i, j]
                end
            end

        end

        # Loop for the Robin boundary for outer wall
        for (cellid, faceid) in getfacetset(grid, "right")
            # Reset the buffers
            fill!(Ke, 0)
            fill!(fe, 0)

            coords = getcoordinates(grid, cellid)

            # Update FE basis
            reinit!(facetvalues, getcoordinates(grid, cellid), faceid)
            # Compute the local contribution
            for qp in 1:getnquadpoints(facetvalues)
                dΓ = getdetJdV(facetvalues, qp)
                coords_qp = spatial_coordinate(facetvalues, qp, coords)
                x = coords_qp[1]
                y = coords_qp[2]
                r = sqrt(x^2 + y^2)
                for i in 1:getnbasefunctions(facetvalues)
                    ϕᵢ = shape_value(facetvalues, qp, i)
                    fe[i] += ( ϕᵢ * htcf4 * Tf4 ) * dΓ
                    for j in 1:getnbasefunctions(facetvalues)
                        ϕⱼ = shape_value(facetvalues, qp, j)
                        Ke[i, j] += ( ϕᵢ * htcf4 * ϕⱼ ) * dΓ
                    end
                end
            end

            # Assemble local contribution
            dofs = celldofs(dh, cellid)
            for (i, I) in pairs(dofs)
                f[I] += fe[i]
                for (j, J) in pairs(dofs)
                    K[I, J] += Ke[i, j]
                end
            end

        end

        return K, f
    end

    # In addition to the diffusive part, we also need a function that assembles the mass matrix `M`.
    function doassemble_M!(M::SparseMatrixCSC, cellvalues::CellValues, dh::DofHandler)

        n_basefuncs = getnbasefunctions(cellvalues)
        Me = zeros(n_basefuncs, n_basefuncs)

        assembler = start_assemble(M)

        for cell in CellIterator(dh)

            fill!(Me, 0)
            coords = getcoordinates(cell)

            reinit!(cellvalues, cell)

            for q_point in 1:getnquadpoints(cellvalues)
                dΩ = getdetJdV(cellvalues, q_point)
                coords_qp = spatial_coordinate(cellvalues, q_point, coords)
                x = coords_qp[1]
                y = coords_qp[2]
                r = sqrt(x^2 + y^2)
                for i in 1:n_basefuncs
                    v = shape_value(cellvalues, q_point, i)
                    for j in 1:n_basefuncs
                        u = shape_value(cellvalues, q_point, j)
                        Me[i, j] += property(r,tank_radius+wall_thickness,rho_steel*cp_steel,rho_insulation*cp_insulation)*(v * u) * dΩ
                    end
                end
            end

            assemble!(assembler, celldofs(cell), Me)
        end
        return M
    end

    function compute_heat_fluxes(cellvalues::CellValues, dh::DofHandler, a::AbstractVector{T}) where {T}

        n = getnbasefunctions(cellvalues)
        cell_dofs = zeros(Int, n)
        nqp = getnquadpoints(cellvalues)

        # Allocate storage for the fluxes to store
        q = [Vec{2, T}[] for _ in 1:getncells(dh.grid)]

        for (cell_num, cell) in enumerate(CellIterator(dh))
            q_cell = q[cell_num]
            celldofs!(cell_dofs, dh, cell_num)
            aᵉ = a[cell_dofs]
            reinit!(cellvalues, cell)

            for q_point in 1:nqp
                q_qp = - function_gradient(cellvalues, q_point, aᵉ)
                push!(q_cell, q_qp)
            end
        end

        return q

    end

    uₙ = zeros(ndofs(dh))
    u = zeros(ndofs(dh))
    fill!(uₙ, Twall_initial - 273.15)

    wtite_VTK = false
    # To store the solution, we initialize a paraview collection (.pvd) file,
    if wtite_VTK
        pvd = paraview_collection("transient-heat")
        VTKGridFile("transient-heat-0", dh) do vtk
            write_solution(vtk, dh, uₙ)
            pvd[0.0] = vtk
        end
    end

    # FEM wall finished

    level_min = 0.01
    alpha_lag = 0.5
    volume_fraction_noz = 1.0

    Cd = 0.85
    d1 = 2.0*25.4/1000.0/sqrt(2)
#    d1 = 10/1000.0
    d2 = 8.0*25.4/1000.0

    number_tanks_parallel = 1
    println("number of tanks parallel = $(number_tanks_parallel)")

#    flowout_vap1, pcrit_vap, Tcrit_vap, Tout_vap = orificeFlow(model,holdups[1].p,holdups[1].T,holdups[1].z,pcrit_vap,Tcrit_vap,pout_vap,Tout_vap,d1,d2,Cdliq,Cdvap,false)
#    flowout_liq1, pcrit_liq, Tcrit_liq, Tout_liq = orificeFlow(model,holdups[2].p,holdups[2].T,holdups[2].z,pcrit_liq,Tcrit_liq,pout_liq,Tout_liq,d1,d2,Cdliq,Cdvap,false)

    flowout_vap = 0.0
    flowout_liq = 0.0

#    println("pcrit_vap = $(pcrit_vap)")
#    println("pcrit_liq = $(pcrit_liq)")

    pcrit_ratio_vap = 0.58
    Tcrit_ratio_vap = 0.934
    pcrit_ratio_liq = 0.70
    Tcrit_ratio_liq = 0.965

    istep = 0
    iVTK = 0
    Δp = Δpbase
    Δt = 0.0
    Δv = 0.0

    CONDENSING = 0.0
    BOILING = 0.0

    massMeOH_required = 0.0
    flowin_MeOH_mass = 0.0

    while (ps > pe)

        istep += 1
        iVTK += 1

        println("")
        println("Step: $(istep)")

        ps -= Δp

        Δv_old = Δv
        Δv_error = 1.0
        iter = 0

        for ih = eachindex(holdups)
            holdups[ih].p = ps
        end

    #    for ih = 1:2
    #        holdups[ih].T += ΔT[ih]
    #        holdups_last[ih].T += ΔT[ih]
    #    end

        tiny_holdup = 1.0e-6

    #    while(Δv_error > 1.0e-4 || iter <= 2)

    #        iter += 1
    #        println("Iteration: $(iter)")

            Twall_noz = Twall_noz_last

        #    println("holdups_last[$(1)] = $(holdups_last[1])")
        #    println("holdups_last[$(2)] = $(holdups_last[2])")
        #    println("holdups_last[$(3)] = $(holdups_last[3])")

            # zero the moles and enathalpy addition of N2
            Δcm = Vector{Vector{Float64}}(undef,0)
            for ih = eachindex(holdups)
                cm = zeros(length(holdups[ih].z))
                push!(Δcm,cm)
            end

            Δmh = Vector{Float64}(undef,0)
            for ih = eachindex(holdups)
                mh = 0.0
                push!(Δmh,mh)
            end
        
            holdups_mh = Vector{Float64}(undef,0)
            for ih = eachindex(holdups)
                push!(holdups_mh,holdups_last[ih].moles*holdups_last[ih].h)
            end

            # here we add in mh from flow into a phase from outside the vessel, like N2

            for ih = eachindex(holdups)
                holdups_mh[ih] += ΔQ[ih] - holdups_last[ih].volume*Δp
            end

            holdups_cm = Vector{Vector{Float64}}(undef,0)
            for ih = eachindex(holdups)
                cm = Vector{Float64}(undef,0)
                for ic = eachindex(holdups[ih].z)
                    push!(cm,holdups_last[ih].moles*holdups_last[ih].z[ic])
                end
                push!(holdups_cm,cm)
            end

            # zero the moles and enathalpy transfers between phases
            Δcm = Vector{Vector{Float64}}(undef,0)
            for ih = eachindex(holdups)
                cm = zeros(length(holdups[ih].z))
                push!(Δcm,cm)
            end

            Δmh = Vector{Float64}(undef,0)
            for ih = eachindex(holdups)
                mh = 0.0
                push!(Δmh,mh)
            end

            # assume vapour removal true, liquid removal false
            if holdups[1].p > 1.0e5
                topflow = false
            else
                topflow = true
            end

            # add in MeOH
            p_MeOH = 1.01325e5
            T_MeOH = 47.5+273.15
                
            # get enthalpy of inlet stream for mixing, don't need to know temperature at outlet. this would be the isenthalpic flash
            # but delta H is zero so just mix inlet enthalpy
            flashMeOH = Clapeyron.tp_flash_impl(modelMeOH,p_MeOH,T_MeOH,zMeOH, HELDTPFlash(verbose = false))
            βMeOH,mwMeOH,mwtMeOH,vMeOH,vtNMeOH,hMeOH,htNMeOH,sMeOH,stMeOH,xMeOH  = get_props(modelMeOH,flashMeOH)

            if holdups[1].p < 27.0e5 && holdups[1].p >= 5.25e5
                flowin_vap_MeOH = 0.125
                flowin_liq_MeOH = 0.125    
            elseif holdups[1].p < 5.25e5
                flowin_vap_MeOH = 0.125
                flowin_liq_MeOH = 0.875
            else
                flowin_vap_MeOH = 0.0
                flowin_liq_MeOH = 0.0
            end

            flowin_MeOH_mass = (flowin_vap_MeOH + flowin_liq_MeOH)*mwtMeOH

            println("MeOH Flow = $(flowin_vap_MeOH + flowin_liq_MeOH) mol/s")
            println("MeOH Flow = $(flowin_MeOH_mass) kg/s")

            # only add MeOH if there is holdup there
            if holdups_last[1].moles > tiny_holdup
                for ic = eachindex(zMeOH)
                    Δcm[1][ic] += Δt*flowin_vap_MeOH*zMeOH[ic]
                end
                Δmh[1] += Δt*flowin_vap_MeOH*htNMeOH
            end

             if holdups_last[2].moles > tiny_holdup
                for ic = eachindex(zMeOH)
                    Δcm[2][ic] += Δt*flowin_liq_MeOH*zMeOH[ic]
                end
                Δmh[2] += Δt*flowin_liq_MeOH*htNMeOH
            end

            # perform phase separation based on current holdups

            # fluid phases

            vapour_eff = 1.0
            liquid_eff = 1.0

            #=
            r_droptlet = 2.5/1000.0
            vol_droplet = 4.0/3.0*pi*(r_droptlet)^3
            vel_droplet = zeros(2)
            for ih = 1:2
                rho = zeros(length(holdups[ih].phases))
                mu = zeros(length(holdups[ih].phases))
                for ip = eachindex(holdups[ih].phases)
                    rho[ip] = 1.0/holdups[ih].phases[ip].v # mol/m3
                    mu[ip] = SuperTRAPP_mu(model_params, modelC3H6_params,holdups[ih].T,rho[ip],holdups[ih].phases[ip].z) # Ps.s
                    println("rho = $(rho[ip]*holdups[ih].phases[ip].mw) kg/m3")
                    println("mu = $(mu[ip]) Pa.s")
                end
                
                if length(holdups[ih].phases) > 1
                    vel_droplet[ih] = vol_droplet*(rho[2]*holdups[ih].phases[2].mw - rho[1]*holdups[ih].phases[1].mw)*9.81/(6.0*pi*mu[ih]*r_droptlet)
                end
            end

            println("vel_droplet = $(vel_droplet) m/s")

            distance_droplet = zeros(2)
            for ih = 1:2
                distance_droplet[ih] = vel_droplet[ih]*Δt
            end

            if height[1] > 0.0
                vapour_eff = max(min(distance_droplet[1]/height[1],1.0),0.001)
            else
                vapour_eff = 1.0
            end
            println("vapour_eff = $(vapour_eff)")
            if height[2] > 0.0
                liquid_eff = max(min(distance_droplet[2]/height[2],1.0),0.001)
            else
                liquid_eff = 1.0
            end
            println("liquid_eff = $(liquid_eff)")
            =#

            mt1 = holdups[1].moles
            if length(holdups[1].phases) > 1
                β2 = holdups[1].phases[2].β*vapour_eff
                z2 = holdups[1].phases[2].z
                h2 = holdups[1].phases[2].h
                for ic = eachindex(holdups[1].z)
                    Δcm[1][ic] -= β2*mt1*z2[ic]
                    Δcm[2][ic] += β2*mt1*z2[ic]
                end
                Δmh[1] -= β2*mt1*h2
                Δmh[2] += β2*mt1*h2
            else
                phase_type = Clapeyron.VT_identify_phase(model, holdups[1].v, holdups[1].T, holdups[1].z)
                if phase_type == :liquid
                    # all liquid β2 will be 1 and phase disappears
                    β2 = holdups[1].phases[2].β*vapour_eff
                    z2 = holdups[1].phases[2].z
                    h2 = holdups[1].phases[2].h
                    for ic = eachindex(holdups[1].z)
                        Δcm[1][ic] -= β2*mt1*z2[ic]
                        Δcm[2][ic] += β2*mt1*z2[ic]
                    end
                    Δmh[1] -= β2*mt1*h2
                    Δmh[2] += β2*mt1*h2
                end
            end

        #    println("Δcm = $(Δcm)")
        #    println("Δmh = $(Δmh)")

            if holdups[1].phase_solid.moles > 0.0
                βs = holdups[1].phase_solid.β
                zs = holdups[1].phase_solid.z
                hs = holdups[1].phase_solid.h
                for ic = eachindex(holdups[1].phase_solid.z)
                    Δcm[1][ic] -= βs*mt1*zs[ic]
                    if holdups[2].moles > 0.0
                        Δcm[2][ic] += βs*mt1*zs[ic]
                    else
                        Δcm[3][ic] += βs*mt1*zs[ic]
                    end
                end
                Δmh[1] -= βs*mt1*hs
                if holdups[2].moles > 0.0
                    Δmh[2] += βs*mt1*hs
                else
                    Δmh[3] += βs*mt1*hs
                end
            end

        #    println("Δcm = $(Δcm)")
        #    println("Δmh = $(Δmh)")

            mt2 = holdups[2].moles
            if length(holdups[2].phases) > 1
                β1 = holdups[2].phases[1].β*liquid_eff
                z1 = holdups[2].phases[1].z
                h1 = holdups[2].phases[1].h
                for ic = eachindex(holdups[2].z)
                    Δcm[1][ic] += β1*mt2*z1[ic]
                    Δcm[2][ic] -= β1*mt2*z1[ic]
                end
                Δmh[1] += β1*mt2*h1
                Δmh[2] -= β1*mt2*h1
            else
                phase_type = Clapeyron.VT_identify_phase(model, holdups[2].v, holdups[2].T, holdups[2].z)
                if phase_type == :vapour
                    # all vapour β1 will be 1 and phase disappears
                    β1 = holdups[2].phases[1].β*liquid_eff
                    z1 = holdups[2].phases[1].z
                    h1 = holdups[2].phases[1].h
                    for ic = eachindex(holdups[2].z)
                        Δcm[1][ic] += β1*mt2*z1[ic]
                        Δcm[2][ic] -= β1*mt2*z1[ic]
                    end
                    Δmh[1] += β1*mt2*h1
                    Δmh[2] -= β1*mt2*h1
                end
            end

        #    println("Δcm = $(Δcm)")
        #    println("Δmh = $(Δmh)")

            if holdups[2].phase_solid.moles > 0.0
                βs = holdups[2].phase_solid.β
                zs = holdups[2].phase_solid.z
                hs = holdups[2].phase_solid.h
                for ic = eachindex(holdups[2].phase_solid.z)
                    Δcm[2][ic] -= βs*mt2*zs[ic]
                    Δcm[3][ic] += βs*mt2*zs[ic]
                end
                Δmh[2] -= βs*mt2*hs
                Δmh[3] += βs*mt2*hs
            end

        #    println("Δcm = $(Δcm)")
        #    println("Δmh = $(Δmh)")

            mt3 = holdups[3].moles
            if holdups[3].phases[1].moles > 0.0
                β1 = holdups[3].phases[1].β
                z1 = holdups[3].phases[1].z
                h1 = holdups[3].phases[1].h
                for ic = eachindex(holdups[3].phases[1].z)
                    Δcm[1][ic] += β1*mt3*z1[ic]
                    Δcm[3][ic] -= β1*mt3*z1[ic]
                end
                Δmh[1] += β1*mt3*h1
                Δmh[3] -= β1*mt3*h1
            end

        #    println("Δcm = $(Δcm)")
        #    println("Δmh = $(Δmh)")

            # mix phases
            for ih = eachindex(holdups)
                holdups_mh[ih] += Δmh[ih]
            end

            for ih = eachindex(holdups)
                for ic = eachindex(holdups_cm[ih])
                    holdups_cm[ih][ic] += Δcm[ih][ic]
                end
            end

            for ih = eachindex(holdups)
                sum_cm = 0.0
                for ic = eachindex(holdups[ih].z)
                    sum_cm += holdups_cm[ih][ic]
                end
                holdups[ih].moles = sum_cm
            end

            for ih = eachindex(holdups)
                if abs(holdups[ih].moles) > tiny_holdup
                    for ic = eachindex(holdups[ih].z)
                        holdups[ih].z[ic] = holdups_cm[ih][ic]/holdups[ih].moles
                    end
                end
            end

            for ih = eachindex(holdups)
                if abs(holdups[ih].moles) > tiny_holdup
                    holdups[ih].h = holdups_mh[ih]/holdups[ih].moles
                end
            end

        #    println("flash in2 holdups[$(1)] = $(holdups[1])")
        #    println("flash in2 holdups[$(2)] = $(holdups[2])")
        #    println("flash in2 holdups[$(3)] = $(holdups[3])")

            # flash after separation and mixing to get consistent holdups. Note this may results in some phases appearing
            # we don't separate further, pragmatically we leave for next step

            # phflash fluid phases
            phFlashFluids!(holdups,ΔT,tiny_holdup,false)
            for ih = eachindex(holdups)
                ΔT[ih] = holdups[ih].T - holdups_last[ih].T
            end

            if length(holdups[1].phases) > 1
                CONDENSING = 1.0
                if Twall_vap+273.15 < holdups[1].T
                    println("CONDENSING in bulk and on wall")
                else
                    println("CONDENSING in bulk")
                end
            else
                CONDENSING = 0.0
                println("No CONDENSING")
            end

            if length(holdups[2].phases) > 1
                BOILING = 1.0
                if Twall_liq+273.15 > holdups[2].T
                    println("BOILING in bulk and on wall")
                else
                    println("BOILING in bulk")
                end
            else
                BOILING = 0.0
                println("No BOILING")
            end

            # phflash solids
            phFlashSolid!(holdups,tiny_holdup,false)

        #    println("flash out2 holdups[$(1)] = $(holdups[1])")
        #    println("flash out2 holdups[$(2)] = $(holdups[2])")
        #    println("flash out2 holdups[$(3)] = $(holdups[3])")

            tank_volume_new = holdups[1].volume + holdups[2].volume +  holdups[3].volume

            Δv = tank_volume_new - tank_volume

            println("Δv = $(Δv)")

       #     if Δv_old < 1.0
       #         Δv_error = abs(Δv - Δv_old)
       #     else
       #         Δv_error = abs(Δv/Δv_old - 1.0)
       #    end
       #     println("Δv_error = $(Δv_error)")
       #     Δv_old = Δv

    #    end # end of iteration to set volume of phases = tank volume

    #    println("flow in holdups[1] = $(holdups[1])")

    #    println("pcrit_ratio_vap = $(pcrit_ratio_vap)")
    #    flowout_vap, pcrit_vap, Tcrit_vap = OrificeFlow(model,holdups[1].p,holdups[1].T,holdups[1].z,pout_vap,pcrit_ratio_vap,Tcrit_ratio_vap,Cd,d1,d2,true)

        CV_vap = 0.019
        flowout_vap = ValveFlow(holdups[1].p,pout_vap,pcrit_ratio_vap,holdups[1].mw/holdups[1].v,CV_vap)
        println("rhoin_vap = $(holdups[1].mw/holdups[1].v)")
        println("flowout_vap = $(flowout_vap)")
    #    println("CV_vap = $(flowout_vap/flowoutvalve_vap)")

        flowout_vap /= number_tanks_parallel
        println("flowout_vap = $(flowout_vap)")

    #    pcrit_ratio_vap = pcrit_vap/holdups[1].p
    #    Tcrit_ratio_vap = Tcrit_vap/holdups[1].T
    #    println("p,T crit ratio vap = $(pcrit_ratio_vap),$(Tcrit_ratio_vap) ")

    #    println("flow in holdups[2] = $(holdups[2])")
        
    #    println("pcrit_ratio_liq = $(pcrit_ratio_liq)")
    #    flowout_liq, pcrit_liq, Tcrit_liq = OrificeFlow(model,holdups[2].p,holdups[2].T,holdups[2].z,pout_liq,pcrit_ratio_liq,Tcrit_ratio_liq,Cd,d1,d2,true)

        CV_liq = 0.014
        flowout_liq = ValveFlow(holdups[2].p,pout_liq,pcrit_ratio_liq,holdups[2].mw/holdups[2].v,CV_liq)
        println("rhoin_liq = $(holdups[2].mw/holdups[2].v)")
        println("flowout_liq = $(flowout_liq)")
    #    println("CV_liq = $(flowout_liq/flowoutvalve_liq)")
        
        flowout_liq /= number_tanks_parallel
        println("flowout_liq = $(flowout_liq)")
        
    #    pcrit_ratio_liq = pcrit_liq/holdups[2].p
    #    Tcrit_ratio_liq = Tcrit_liq/holdups[2].T
    #    println("p,T crit ratio liq = $(pcrit_ratio_liq),$(Tcrit_ratio_liq) ")

        # new level
        for ih = eachindex(holdups)
            tank_xsa_holdups[ih] = holdups[ih].volume/tank_length
        end

        volumeFraction = tank_xsa_holdups[3]/(pi*tank_radius^2)

        if tank_xsa_holdups[3] > 10.0*eps(Float64)
            if theta_solid == 0.0
                theta_solid = sqrt(0.9*volumeFraction)*360/180*pi
            end
            error = 1.0
            while (error > 0.0001)
                f =  (tank_radius^2)*(theta_solid - sin(theta_solid)) - 2.0*tank_xsa_holdups[3]
                df_dtheta = (tank_radius^2)*(1.0 - cos(theta_solid))
                theta_solid = theta_solid - f/df_dtheta
                error = abs(f)
            end
        else
            theta_solid = 0.0
        end

        # height based on tank diamter
        tank_diameter = 2.0*tank_radius
        solidPercent = (1.0 - cos(theta_solid/2.0))/2.0*100.0
        height_solid = tank_diameter*solidPercent/100.0

        cord_solid = 2.0*tank_radius*sin(theta_solid/2.0)
        interface_area_solid = cord_solid*tank_length

        volumeFraction = (tank_xsa_holdups[2] + tank_xsa_holdups[3])/(pi*tank_radius^2)

        if tank_xsa_holdups[2] + tank_xsa_holdups[3] > 10.0*eps(Float64)
            if theta == 0.0
                theta = sqrt(0.9*volumeFraction)*360/180*pi
            end
            error = 1.0
            while (error > 0.0001)
                f =  (tank_radius^2)*(theta - sin(theta)) - 2.0*(tank_xsa_holdups[2] + tank_xsa_holdups[3])
                df_dtheta = (tank_radius^2)*(1.0 - cos(theta))
                theta = theta - f/df_dtheta
                error = abs(f)
            end
        else
            theta = 0.0
        end

        # height based on tank diamter
        tank_diameter = 2.0*tank_radius
        levelPercent = (1.0 - cos(theta/2.0))/2.0*100.0
        height[2] = tank_diameter*levelPercent/100.0
        height[1] = tank_diameter*(1.0 - levelPercent/100.0)

        cord = 2.0*tank_radius*sin(theta/2.0)
        interface_area = cord*tank_length

        if topflow
            level_noz = levelPercent/100.0
            #volume_fraction_noz = alpha_lag*max(min((1.0 - level_noz)/level_min,1.0),0.0) + (1.0 - alpha_lag)*volume_fraction_noz
            volume_fraction_noz = max(min((1.0 - level_noz)/level_min,1.0),0.0)
            println("volume_fraction_noz = $(volume_fraction_noz)")
            Δvolume[1] = min(Δv,holdups[1].volume*volume_fraction_noz)
            Δvolume[2] = Δv - Δvolume[1]
            Δt_liq = Δvolume[2]/(flowout_liq*holdups[2].v)
            Δt_vap = Δvolume[1]/(flowout_vap*holdups[1].v)
            holdups[2].volume -= Δvolume[2]
            Δmoles[2] = max(holdups[2].moles - holdups[2].volume/holdups[2].v,0.0)
            holdups[2].moles -= Δmoles[2]
            holdups[1].volume -= Δvolume[1]
            Δmoles[1] = max(holdups[1].moles - holdups[1].volume/holdups[1].v,0.0)
            holdups[1].moles -= Δmoles[1]
        else
            level_noz = (levelPercent - solidPercent)/100.0
            #volume_fraction_noz = alpha_lag*max(min(level_noz/level_min,1.0),0.0) + (1.0 - alpha_lag)*volume_fraction_noz
            volume_fraction_noz = max(min(level_noz/level_min,1.0),0.0)
            println("volume_fraction_noz = $(volume_fraction_noz)")
            Δvolume[2] = min(Δv,holdups[2].volume*volume_fraction_noz)
            Δvolume[1] = Δv - Δvolume[2]
            Δt_liq = Δvolume[2]/(flowout_liq*holdups[2].v)
            Δt_vap = Δvolume[1]/(flowout_vap*holdups[1].v)
            holdups[2].volume -= Δvolume[2]
            Δmoles[2] = max(holdups[2].moles - holdups[2].volume/holdups[2].v,0.0)
            holdups[2].moles -= Δmoles[2]
            holdups[1].volume -= Δvolume[1]
            Δmoles[1] = max(holdups[1].moles - holdups[1].volume/holdups[1].v,0.0)
            holdups[1].moles -= Δmoles[1]
        end

        # update moles and volume
        for ih = 1:2
            for ip = eachindex(holdups[ih].phases)
                holdups[ih].phases[ip].moles  = holdups[ih].phases[ip].β*holdups[ih].moles
                holdups[ih].phases[ip].volume = holdups[ih].phases[ip].moles*holdups[ih].phases[ip].v
            end
            holdups[ih].phase_solid.moles  = holdups[ih].phase_solid.β*holdups[ih].moles
            holdups[ih].phase_solid.volume = holdups[ih].phase_solid.moles*holdups[ih].phase_solid.v
        end

    #    println("holdups[$(1)] = $(holdups[1])")
    #    println("holdups[$(2)] = $(holdups[2])")
    #    println("holdups[$(3)] = $(holdups[3])")
    #    println("holdup volume = $(holdups[1].volume + holdups[2].volume + holdups[3].volume)")

        Δt = Δt_vap + Δt_liq

        println("Δt_vap = $(Δt_vap)")
        println("Δt_liq = $(Δt_liq)")

        dia_min = 0.010
        dia_vap = max(sqrt(4.0/pi*tank_xsa_holdups[1]), dia_min)
        if tank_xsa_holdups[1] > pi/4.0*dia_min*dia_min
            vel_vap = (Δmoles[1]/Δt)*holdups[1].v/tank_xsa_holdups[1]
        else
            vel_vap = (Δmoles[1]/Δt)*holdups[1].v/(pi/4.0*dia_min*dia_min)
        end

        mu_vap = 0.0
        lambda_vap = 0.0
        cp_vap = 0.0
        for ip = eachindex(holdups[1].phases)
            rho_vap = 1.0/holdups[1].phases[ip].v # mol/m3
            beta_vap = holdups[1].phases[ip].β
            mu_vap += beta_vap*(SuperTRAPP_mu(model_params, modelC3H6_params,holdups[1].T,rho_vap,holdups[1].phases[ip].z)) # Ps.s
            lambda_vap += beta_vap*(SuperTRAPP_lambda(pures_model, model_params, modelC3H6_params,holdups[1].T,rho_vap,holdups[1].phases[ip].z)) # W/m/K
            cp_vap += beta_vap*(Clapeyron.VT_isobaric_heat_capacity(model,holdups[1].phases[ip].v,holdups[1].T,holdups[1].phases[ip].z)) # J/Kg/K
        end

        cp_vap /= holdups[1].mw

    #    println("mu_vap = $(mu_vap)")
    #    println("lambda_vap = $(lambda_vap)")
    #    println("cp_vap = $(cp_vap)")

        function alpha_vap(Tw, T, u, d, L, rho, mu, cp , lambda)
            g = 9.8065
            Pr = max(cp*mu/lambda,0.01)
        #    println("Pr = $(cp*mu/lambda)")
            Re = rho*u*d/mu
        #    println("Re = $(rho*u*d/mu)")
            NuLam = 3.657 + 0.065*Re*Pr*d/L/(1.0 + 0.04*(Re*Pr*d/L)^(2/3))
            NuTurb = 0.0
            if Re > 1000.0
                f = 1.82 * log10(Re) - 1.64
                f = 1.0/f/f
                NuTurb = f/8.0*(Re - 1000.0)*Pr/(1.0 + 12.7*sqrt(f/8.0)*(Pr^(2/3) - 1.0))*(1.0 + (d/L)^(2/3))
            end
            beta = 2.0/(Tw + T)
            Gr = d*d*d*g*rho*rho*beta*abs(Tw - T)/mu/mu
            NuNC = 0.825 + 0.387*(Gr*Pr / (1 + (0.492/Pr)^(9/16))^(16/9))^(1/6)
            NuNC = NuNC * NuNC
            Nu = sqrt(NuNC * NuNC + NuLam * NuLam + NuTurb * NuTurb)
            return lambda*Nu/d
        end

        density_vap = 1.0/holdups[1].v*holdups[1].mw
        htc_vap = alpha_vap(Twall_vap+273.15, holdups[1].T, vel_vap, dia_vap, tank_length/2, density_vap, mu_vap, cp_vap ,lambda_vap)

        if holdups[1].volume < pi/4.0*dia_min*dia_min*tank_length
            htc_vap = 0.0
        end
        
        println("htc_vap = $(htc_vap)")

        dia_liq = max(sqrt(4.0/pi*tank_xsa_holdups[2]), dia_min)
        if tank_xsa_holdups[2] > pi/4.0*dia_min*dia_min
            vel_liq = (Δmoles[2]/Δt)*holdups[2].v/tank_xsa_holdups[2]
        else
            vel_liq = (Δmoles[2]/Δt)*holdups[2].v/(pi/4.0*dia_min*dia_min)
        end

        mu_liq = 0.0
        lambda_liq = 0.0
        cp_liq = 0.0
        for ip = eachindex(holdups[2].phases)
            rho_liq = 1.0/holdups[2].phases[ip].v # mol/m3
            beta_liq = holdups[2].phases[ip].β
            mu_liq += beta_liq*(SuperTRAPP_mu(model_params, modelC3H6_params,holdups[2].T,rho_liq,holdups[2].phases[ip].z)) # Ps.s
            lambda_liq += beta_liq*(SuperTRAPP_lambda(pures_model, model_params, modelC3H6_params,holdups[2].T,rho_liq,holdups[2].phases[ip].z)) # W/m/K
            cp_liq += beta_liq*(Clapeyron.VT_isobaric_heat_capacity(model,holdups[2].phases[ip].v,holdups[2].T,holdups[2].phases[ip].z)) # J/Kg/K
        end

        cp_liq /= holdups[2].mw

    #    println("mu_liq = $(mu_liq)")
    #    println("lambda_liq = $(lambda_liq)")
    #    println("cp_liq = $(cp_liq)")

        function alpha_liq(Tw, T, u, d, L, rho, mu, cp , lambda, BOILING)
            g = 9.8065
            Pr = max(cp*mu/lambda,0.01)
        #    println("Pr = $(cp*mu/lambda)")
            Re = rho*u*d/mu
        #    println("Re = $(rho*u*d/mu)")
            NuLam = 3.657 + 0.065*Re*Pr*d/L/(1.0 + 0.04*(Re*Pr*d/L)^(2/3))
            NuTurb = 0.0
            if Re > 1000.0
                f = 1.82*log10(Re) - 1.64
                f = 1.0/f/f
                NuTurb = f/8.0*(Re - 1000.0)*Pr/(1.0 + 12.7*sqrt(f/8.0)*(Pr^(2/3) - 1.0))*(1.0 + (d/L)^(2/3))
            end
            deltaT0 = 20000.0 / 3500.0
            deltaT = Tw - T
            if deltaT < 0.0 
                deltaT = 0.0
            end
            NuBOIL = d / lambda * 3500.0 * BOILING * (deltaT / deltaT0)^0.25
            beta = 2.0/(Tw + T)
            Gr = d*d*d*g*rho*rho*beta*abs(Tw - T)/mu/mu
            NuNC = 0.825 + 0.387*(Gr*Pr / (1 + (0.492/Pr)^(9/16))^(16/9))^(1 /6)
            NuNC = NuNC * NuNC
            Nu = sqrt(NuNC * NuNC + NuLam * NuLam + NuTurb * NuTurb + NuBOIL * NuBOIL)
            return lambda*Nu/d
        end

        density_liq = 1.0/holdups[2].v*holdups[2].mw
        htc_liq = alpha_liq(Twall_liq+273.15, holdups[2].T, vel_liq, dia_liq, tank_length/2, density_liq, mu_liq, cp_liq ,lambda_liq, BOILING)
    
        if holdups[2].volume < pi/4.0*dia_min*dia_min*tank_length
            htc_liq = 0.0
        end

        println("htc_liq = $(htc_liq)")

        dia_solid = max(sqrt(4.0/pi*tank_xsa_holdups[3]), dia_min)
        Nu_solid = 1.0
        pratio = ps/200/1e6
        A0 =  12.1860*pratio*pratio - 75.579*pratio + 166.6
        A1 = -(9.2052*pratio*pratio - 46.508*pratio + 98.634)
        A2 =  2.2499*pratio*pratio - 9.5547*pratio + 19.565
        A3 = -(0.1788*pratio*pratio - 0.6566*pratio + 1.3047)
        lambda_solid = exp(A0 + A1*log(Twall_solid+273.15) + A2*log(Twall_solid+273.15)^2 + A3*log(Twall_solid+273.15)^3) # W/m/K
    #    println("lambda_solid = $(lambda_solid)")
        htc_solid = Nu_solid*lambda_solid/dia_solid

        if holdups[3].volume < pi/4.0*dia_min*dia_min*tank_length
            htc_solid = 0.0
        end

        println("htc_solid = $(htc_solid)")

        function alpha_ambient(Tw, T, u, d, rho, mu, cp , lambda)
            g = 9.8065
            Re = rho*u*d/mu
            beta = 2.0/(Tw + T)
            Gr = d*d*d*g*rho*rho*beta*abs(Tw - T)/mu/mu
            Pr = cp*mu/lambda
            Re = sqrt(Re*Re + Gr/2.5)
            NuLam = 0.664*sqrt(Re)*Pr^(1/3)
            NuTurb = 0.037*(Re^0.9)*Pr/(Re^0.1 + 2.443 * (Pr^(2/3) - 1.0))
            Nu = 0.3 + sqrt(NuLam * NuLam + NuTurb * NuTurb)
            return lambda*Nu/d
        end

        
        htc_ambient = alpha_ambient(Twall_outer+273.15, T_ambient+273.15, u_ambient, tank_diameter, rho_ambient*mwAir, mu_ambient, cp_ambient , lambda_ambient)

        println("htc_ambient = $(htc_ambient)")

        function alpha_noz(Tw, T, u, d, rho, mu, cp , lambda)
            g = 9.8065
            Pr = cp*mu/lambda
            Re = rho*u*d/mu
            NuLam = 0.0
            NuTurb = 0.0
            if Re > 1.0
                NuLam = 3.657
                f = 1.82 * log10(Re) - 1.64
                f = 1.0/f/f
                NuTurb = f/8.0*Re*Pr/(1.0 + 12.7*sqrt(f/8.0)*(Pr^(2/3) - 1.0))
            end
            beta = 2.0/(Tw + T)
            Gr = d*d*d*g*rho*rho*beta*abs(Tw - T)/mu/mu
            NuNC = 0.825 + 0.387*(Gr*Pr / (1 + (0.492/Pr)^(9/16))^(16/9))^(1/6)
            NuNC = NuNC * NuNC
            Nu = sqrt(NuNC * NuNC + NuLam * NuLam + NuTurb * NuTurb)
        return lambda*Nu/d

        end

        if Δmoles[1] > tiny_holdup
            vel_noz_vap = flowout_vap*holdups[1].v/(pi/4*d2^2)
        else
            vel_noz_vap = 0.0
        end
        htc_noz_vap = alpha_noz(Twall_noz+273.15, holdups[1].T, vel_noz_vap, d2, density_vap, mu_vap, cp_vap ,lambda_vap)
    #    println("htc_noz_vap = $(htc_noz_vap)")
        if Δmoles[2] > tiny_holdup
            vel_noz_liq = flowout_liq*holdups[2].v/(pi/4*d2^2)
        else
            vel_noz_liq = 0.0
        end
        htc_noz_liq = alpha_noz(Twall_noz+273.15, holdups[2].T, vel_noz_liq, d2, density_liq, mu_liq, cp_liq ,lambda_liq)
    #    println("htc_noz_liq = $(htc_noz_liq)")
        htc_ambient_noz = alpha_ambient(Twall_noz+273.15, T_ambient+273.15, u_ambient, d2+2.0*overall_thickness, rho_ambient*mwAir, mu_ambient, cp_ambient , lambda_ambient)
        htc_insulation_noz = lambda_insulation/((d2/2.0)*log((d2/2.0+insulation_thickness)/(d2/2.0)))
        htc_ambient_overall = htc_ambient_noz*htc_insulation_noz/(htc_ambient_noz + htc_insulation_noz)
    #    println("htc_ambient_overall = $(htc_ambient_overall)")
        Twall_noz = (d2*Δt_vap*htc_noz_vap*(holdups[1].T - 273.15) + d2*Δt_liq*htc_noz_liq*(holdups[2].T - 273.15) + d2*Δt*htc_ambient_overall*T_ambient + cp_steel*rho_steel*((d2 + 2.0*wall_thickness)^2 - d2^2)/4*Twall_noz)/(cp_steel*rho_steel*((d2 + 2.0*wall_thickness)^2 - d2^2)/4 + d2*Δt_vap*htc_noz_vap + d2*Δt_liq*htc_noz_liq + d2*Δt*htc_ambient_overall)

        # then we use Ferrite to model finite element wall and do a transient step for the Δt

        # Start wall calc using Ferrite
        htcf1 = htc_solid
        Tf1 = holdups[3].T-273.15
        htcf2 = htc_liq
        Tf2 = holdups[2].T-273.15
        htcf3 = htc_vap
        Tf3 = holdups[1].T-273.15
        htcf4 = htc_ambient
        Tf4 = T_ambient
        LL = height[2] - tank_radius
        SL = height_solid - tank_radius

        K = allocate_matrix(dh)
        f = zeros(ndofs(dh))
        M = allocate_matrix(dh)

        K, f = doassemble_Kf!(K, f, cellvalues, dh)
        # need to modify RobinBC for solid and liquid levels SL and LL
        K, f = doassemble_Kf_RobinBC!(K, f, cellvalues,  dh, SL, LL, htcf1, Tf1, htcf2, Tf2, htcf3, Tf3, htcf4, Tf4)
        M = doassemble_M!(M, cellvalues, dh)

        A = (Δt .* K) + M
        b = Δt .* f .+ M * uₙ

        u = A \ b

        length_inner = pi*tank_diamter
        length_outer = pi*(tank_diamter + overall_thickness)

        npoints = Int(floor((length_inner + length_outer)/2/h))

        points = [Vec(((tank_radius+overall_thickness-0.0001)*cos(theta), (tank_radius+overall_thickness-0.0001)*sin(theta))) for theta in range(-pi/2, pi/2, length = npoints)];
        ph = PointEvalHandler(grid, points);
        u_points = evaluate_at_points(ph, dh, u, :u);

        Twall_outer = 0.0
        icount = 0
        for ipoints = 1:npoints
            icount += 1
            Twall_outer += u_points[ipoints]
        end

        if icount > 0
            Twall_outer /= icount
        else
            Twall_outer = T_ambient
        end

    #    println("Twall_outer = $(Twall_outer)")

        q_gp = compute_heat_fluxes(cellvalues, dh, u);
        projector = L2Projector(ip, grid);
        q_projected = project(projector, q_gp, qr);
        points = [Vec((tank_radius*cos(theta), tank_radius*sin(theta))) for theta in range(-pi/2, pi/2, length = npoints)];
        ph = PointEvalHandler(grid, points);
        q_points = evaluate_at_points(ph, projector, q_projected);
        u_points = evaluate_at_points(ph, dh, u, :u);

        # this is inner wall temperature
        Twallmin = u_points[1]
        ipoints_min = 1
        for ipoints = 2:npoints
            if u_points[ipoints] < Twallmin
                Twallmin = u_points[ipoints]
                ipoints_min = ipoints
            end
        end

        Twall_solid = 0.0
        icount_solid = 0
        Twall_liq = 0.0
        icount_liq = 0
        Twall_vap = 0.0
        icount_vap = 0
        for ipoints = 1:npoints
            if points[ipoints][2] < SL
                icount_solid += 1
                Twall_solid += u_points[ipoints]
            elseif points[ipoints][2] >= SL && points[ipoints][2] <= LL
                icount_liq += 1
                Twall_liq += u_points[ipoints]
            else
                icount_vap += 1
                Twall_vap += u_points[ipoints]
            end
        end

        if icount_solid > 0
            Twall_solid /= icount_solid
        else
            Twall_solid = u_points[1]
        end

        if icount_liq > 0
            Twall_liq /= icount_liq
        else
            Twall_liq = u_points[1]
        end

        if icount_vap > 0
            Twall_vap /= icount_vap
        else
            Twall_vap = u_points[npoints]
        end

    #    println("Twall_vap = $(Twall_vap)")
    #    println("Twall_liq = $(Twall_liq)")
    #    println("Twall_solid = $(Twall_solid)")

        perimiter = [tank_radius*theta for theta in range(0.0, pi, length = npoints)];

        heat_flux_mag = zeros(npoints)
        for ipoints = 1:npoints
            heat_flux_mag[ipoints] = sqrt((q_points[ipoints][1])^2 + (q_points[ipoints][2])^2)
        end

        Q_wallsolid = 0.0
        Q_wallliq = 0.0
        Q_wallvap = 0.0

        for ipoints = 2:npoints
            dp = perimiter[ipoints] - perimiter[ipoints-1]
            if points[ipoints][2] < SL
                Q_wallsolid += 0.5*(heat_flux_mag[ipoints-1] + heat_flux_mag[ipoints])*dp*tank_length
            elseif points[ipoints][2] >= SL && points[ipoints][2] <= LL
                Q_wallliq += 0.5*(heat_flux_mag[ipoints-1] + heat_flux_mag[ipoints])*dp*tank_length
            else
                Q_wallvap += 0.5*(heat_flux_mag[ipoints-1] + heat_flux_mag[ipoints])*dp*tank_length
            end
        end

        # factor 2 for both sides of wall
        Q_wallsolid *= 2.0*Δt
        Q_wallliq   *= 2.0*Δt
        Q_wallvap   *= 2.0*Δt

    #    println("Q_wallvap   = $(Q_wallvap)")
    #    println("Q_wallliq   = $(Q_wallliq)")
    #    println("Q_wallsolid = $(Q_wallsolid)")

        # End wall calc

        # find T at outside of steel

        points = [Vec(((tank_radius+wall_thickness)*cos(theta), (tank_radius+wall_thickness)*sin(theta))) for theta in range(-pi/2, pi/2, length = npoints)];
        ph = PointEvalHandler(grid, points);
        u_points = evaluate_at_points(ph, dh, u, :u);

        Twallmin_out = u_points[ipoints_min]
        Twallmin_av = 0.5*(Twallmin_out + Twallmin)

        # htc for interface ~ 1/alpha = 1/alpha_vap + 1/alpha_liq, alpha_liq is not boiling part is just convective
        if htc_vap + htc_liq > 0.0
            htc_vap_liq = htc_vap*htc_liq/(htc_vap + htc_liq)
        else
            htc_vap_liq = 0.0
        end
        if htc_liq + htc_solid > 0.0
            htc_liq_solid = htc_liq*htc_solid/(htc_liq + htc_solid)
        else
            htc_liq_solid = 0.0
        end

        if holdups[1].moles > 0.0 && holdups[2].moles > 0.0
            Q_liqvap = Δt*htc_vap_liq*interface_area*(holdups[2].T - holdups[1].T)
            ΔQ[1] = Q_liqvap
            Q_vapliq = Δt*htc_vap_liq*interface_area*(holdups[1].T - holdups[2].T)
            ΔQ[2] = Q_vapliq
        else
            ΔQ[1] = 0.0
            ΔQ[2] = 0.0
        end

        if holdups[1].moles > 0.0
            ΔQ[1] += Q_wallvap
        end

        if holdups[2].moles > 0.0
            ΔQ[2] += Q_wallliq
        end

        if holdups[3].moles > 0.0 && holdups[2].moles > 0.0
            Q_liqsolid = Δt*htc_liq_solid*interface_area_solid*(holdups[2].T - holdups[3].T)
            ΔQ[3] = Q_liqsolid
            Q_solidliq = Δt*htc_liq_solid*interface_area_solid*(holdups[3].T - holdups[2].T)
            ΔQ[2] += Q_solidliq
        else
            ΔQ[3] = 0.0
        end

        if holdups[3].moles > 0.0
            ΔQ[3] += Q_wallsolid
        end

        println("holdups[$(1)] = $(holdups[1])")
        println("holdups[$(2)] = $(holdups[2])")
        println("holdups[$(3)] = $(holdups[3])")
        println("holdup volume = $(holdups[1].volume + holdups[2].volume + holdups[3].volume)")

        time += Δt
        println("time = $(time/3600)")
        push!(time_Plot, time/3600.0)

        println("p = $(holdups[1].p/1e5)")
        println("Tvap = $(holdups[1].T -273.15)")
        println("Tliq = $(holdups[2].T -273.15)")
        println("Twall_vap = $(Twall_vap)")
        println("Twall_liq = $(Twall_liq)")
        println("Twall_solid = $(Twall_solid)")
        println("Twall_outer = $(Twall_outer)")
        println("Twallmin = $(Twallmin)")
        println("Twallmin outer = $(Twallmin_out)")
        println("Twallmin average = $(Twallmin_av)")
        println("Twall_noz = $(Twall_noz)")

        println("Level Percent = $(levelPercent)")
        println("Solid Percent = $(solidPercent)")

        holdups_last = deepcopy(holdups)

        Twall_noz_last = Twall_noz

        massMeOH_required += flowin_MeOH_mass*Δt

        println("Mass MeOH required = $(massMeOH_required) kg")

        #At the end of the time loop, we set the previous solution to the current one and go to the next time step.
        uₙ .= u

        if wtite_VTK
            if iVTK == nVTK
                VTKGridFile("transient-heat-$istep", dh) do vtk
                        write_solution(vtk, dh, uₙ)
                        pvd[time] = vtk
                end
                iVTK = 0
            end
        end

        println("Δmoles[1] = $(Δmoles[1])")
        println("Δmoles[2] = $(Δmoles[2])")

        println("ΔQ[1] = $(ΔQ[1]) J")
        println("ΔQ[2] = $(ΔQ[2]) J")
        println("ΔQ[3] = $(ΔQ[3]) J")

        if Δmoles[1] > tiny_holdup
            fv_plot = flowout_vap*holdups[1].mw
        else
            fv_plot = 0.0
        end

        if Δmoles[2] > tiny_holdup
            fl_plot = flowout_liq*holdups[2].mw
        else
            fl_plot = 0.0
        end

        push!(fv_Plot, fv_plot)
        push!(fl_Plot, fl_plot)
        push!(MeOH_Plot, flowin_MeOH_mass)

        push!(Twallmin_Plot, Twallmin)
        push!(Twallmin_av_Plot, Twallmin_av)
        push!(Twall_MMDT_Plot, -46.0)
        push!(Twallvap_Plot, Twall_vap)
        push!(Twallliq_Plot, Twall_liq)
        push!(Twallnoz_Plot, Twall_noz)

        push!(T1_Plot, holdups[1].T-273.15)
        push!(T2_Plot, holdups[2].T-273.15)

        push!(p1_Plot, holdups[1].p/1e5)
        push!(p2_Plot, holdups[2].p/1e5)
        push!(ps_Plot, holdups[3].p/1e5)

        push!(m1_Plot, holdups[1].moles*holdups[1].mw/1000)
        push!(m2_Plot, holdups[2].moles*holdups[2].mw/1000)
        push!(ms_Plot, holdups[3].moles*holdups[3].mw/1000)
        push!(ms_kg_Plot, holdups[3].moles*holdups[3].mw)
        push!(MeOH_kg_Plot,massMeOH_required)

        println("solid mass = $(holdups[3].moles*holdups[3].mw) kg")

    end

    if wtite_VTK
        vtk_save(pvd)
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
    savefig(p1,"~/julia/dev/Clapeyron.jl/fig1.png")

    p2 = plot(
    label=["vapour" "liquid" "solid"], 
    [m1_Plot,m2_Plot,ms_Plot],
    [p1_Plot,p2_Plot,ps_Plot],
    xlabel = "Mass [tons]",
    ylabel = "Pressure [bara]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

    display(p2)
    savefig(p2,"~/julia/dev/Clapeyron.jl/fig2.png")

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
    savefig(p3,"~/julia/dev/Clapeyron.jl/fig3.png")

    p4 = plot(
    label=["vapour" "liquid" "vapour wall" "liquid wall" "minimum inner wall" "minimum average wall" "nozzle" "LTCS MMDT"], 
    [time_Plot,time_Plot,time_Plot,time_Plot,time_Plot,time_Plot,time_Plot,time_Plot],
    [T1_Plot,T2_Plot,Twallvap_Plot,Twallliq_Plot,Twallmin_Plot,Twallmin_av_Plot,Twallnoz_Plot,Twall_MMDT_Plot],
    xlabel = "Time [hours]",
    ylabel = "Temperature [deg C]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

    display(p4)
    savefig(p4,"~/julia/dev/Clapeyron.jl/fig4.png")

    p5 = plot(
    label=["vapour" "liquid" "solid"],
    [time_Plot,time_Plot,time_Plot], 
    [m1_Plot,m2_Plot,ms_Plot],
    xlabel = "Time [hours]",
    ylabel = "Mass [tons]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

    display(p5)
    savefig(p5,"~/julia/dev/Clapeyron.jl/fig5.png")

    p5b = plot(
    label=["solid"],
    [time_Plot], 
    [ms_kg_Plot],
    xlabel = "Time [hours]",
    ylabel = "Mass [kg]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

    display(p5b)
    savefig(p5b,"~/julia/dev/Clapeyron.jl/fig5b.png")

    p6 = plot(
    label=["vapour" "liquid"],
    [time_Plot,time_Plot], 
    [fv_Plot,fl_Plot],
    xlabel = "Time [hours]",
    ylabel = "Flow [kg/s]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

    display(p6)
    savefig(p6,"~/julia/dev/Clapeyron.jl/fig6.png")

    p7 = plot(
    label=["MeOH"],
    [time_Plot], 
    [MeOH_Plot],
    xlabel = "Time [hours]",
    ylabel = "Flow [kg/s]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

    display(p7)
    savefig(p7,"~/julia/dev/Clapeyron.jl/fig7.png")

    # write data pout to file

    io=open("tank01_upper_27bara_5percent.dat","w") do io
        for i = eachindex(time_Plot)
            t = time_Plot[i]
            p = p1_Plot[i]
            m1 = m1_Plot[i]
            m2 = m2_Plot[i]
            ms = ms_Plot[i]
            MeOH_kg = MeOH_kg_Plot[i]
            MeOH = MeOH_Plot[i]
            fv = fv_Plot[i]
            fl = fl_Plot[i]
            T1 = T1_Plot[i]
            T2 = T2_Plot[i]
            Tw1 = Twallvap_Plot[i]
            Tw2 = Twallliq_Plot[i]
            Twm = Twallmin_Plot[i]
            Twav = Twallmin_av_Plot[i]
            Twnoz = Twallnoz_Plot[i]
            TMMDT = Twall_MMDT_Plot[i]
            println(io, "$(t),$(p),$(m1),$(m2),$(ms),$(MeOH_kg),$(fv),$(fl),$(MeOH),$(T1),$(T2),$(Tw1),$(Tw2),$(Twm),$(Twav),$(Twnoz),$(TMMDT)")
        end
    end

end

delta_P = 0.125e5

pe = 1.125e5
#pe = 26.75e5

nstep = Int((ps - pe)/delta_P)
println("nstep = $(nstep)")
nVTK = 2
BlowDown(model, ps, Ts, zs, tank_volume, tank_radius, tank_length, level, pe, nstep, nVTK)