using Clapeyron, Ferrite, FerriteGmsh, Gmsh, SparseArrays, WriteVTK, Plots

fluid = ["carbon dioxide","nitrogen","methanol"]
nfluid = length(fluid)
model = sCPA(fluid; idealmodel=AlyLeeIdeal, assoc_options=AssocOptions(combining=:elliott))

modelCO2 = sCPA(["carbon dioxide"];idealmodel=AlyLeeIdeal)
mwCO2 = Clapeyron.molecular_weight(modelCO2,[1])

modelAir = GERG2008(["nitrogen","oxygen"])
zAir = [0.79,0.21]
mwAir = Clapeyron.molecular_weight(modelAir,zAir)

function SuperTRAPP_mu(model, Tin, rhoin, zin)

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

	modelC3H6 = GERG2008(["propane"])

	MwRef = modelC3H6.params.Mw.values[1]
	TcRef = modelC3H6.params.Tc.values[1]
	PcRef = modelC3H6.params.Pc.values[1]
	VcRef = modelC3H6.params.Vc.values[1]*100.0^3 # cm3/mol
	rhocRef = 1000.0/VcRef
	OmegaRef = modelC3H6.params.acentricfactor.values[1]
	ZcRef = PcRef*VcRef/(100.0^3)/R̄/TcRef
#	println("Mw, Tc, Pc, Vc, Omega = $(MwRef), $(TcRef), $(PcRef), $(VcRef), $(OmegaRef), $(ZcRef)")

	Mw = model.params.Mw.values
	Tc = model.params.Tc.values
	Pc = model.params.Pc.values
	Vc = model.params.Vc.values*100.0^3 # cm3/mol
	Zc = zeros(nc)
	for i = eachindex(z)
		Zc[i] = Pc[i]*Vc[i]/(100.0^3)/R̄/Tc[i]
    end
	Omega = model.params.acentricfactor.values

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


function SuperTRAPP_lambda(model, Tin, rhoin, zin)

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

	pures = Clapeyron.split_pure_model(model)

	modelC3H6 = GERG2008(["propane"])

	MwRef = modelC3H6.params.Mw.values[1]
	TcRef = modelC3H6.params.Tc.values[1]
	PcRef = modelC3H6.params.Pc.values[1]
	VcRef = modelC3H6.params.Vc.values[1]*100.0^3 # cm3/mol
	rhocRef = 1000.0/VcRef
	OmegaRef = modelC3H6.params.acentricfactor.values[1]
	ZcRef = PcRef*VcRef/(100.0^3)/R̄/TcRef
#	println("Mw, Tc, Pc, Vc, Omega = $(MwRef), $(TcRef), $(PcRef), $(VcRef), $(OmegaRef), $(ZcRef)")

	Mw = model.params.Mw.values
	Tc = model.params.Tc.values
	Pc = model.params.Pc.values
	Vc = model.params.Vc.values*100.0^3 # cm3/mol
	Zc = zeros(nc)
	for i = eachindex(z)
		Zc[i] = Pc[i]*Vc[i]/(100.0^3)/R̄/Tc[i]
    end
	Omega = model.params.acentricfactor.values

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
xMeOH = 0.0001
x = [xCO2*(1.0 - xMeOH),(1.0-xCO2)*(1.0 - xMeOH),xMeOH]
#x = [xCO2,(1.0-xCO2)]

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
    return gm
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
        return xs_min,xf,f,fnc(xs_min)
    end
end

println("Check dry ice formation")
println("")

Tsolid = -78.85+273.15
psolid = 1.01325e5

xCO2 = 0.997275
xMeOH = 0.475
x = [xCO2*(1.0 - xMeOH),(1.0-xCO2)*(1.0 - xMeOH),xMeOH]

βCO2, z_fluid, flash_fluid, gibbs_min = pT_flash_CO2Solid(model,psolid,Tsolid,x)

println("moles of solid CO2 = $(βCO2)")
println("moles of fluid = $(1 - βCO2)")

println("fluid βf = $(flash_fluid.fractions)")
println("fluid vf = $(flash_fluid.volumes) m3/mol")
for i = eachindex(z_fluid)
    println("z_fluid[$(i)] = $(z_fluid[i])")
end
for i = eachindex(flash_fluid.fractions)
    println("fluid xf[$(i)] = $(flash_fluid.compositions[i])")
end
println("gibbs min = $(gibbs_min)")

println("")
println("Check dry ice formation complete")