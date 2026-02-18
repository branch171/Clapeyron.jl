using Clapeyron, Ferrite, FerriteGmsh, Gmsh, SparseArrays, WriteVTK, Plots

fluid = ["carbon dioxide","nitrogen"]
nfluid = length(fluid)
model = GERG2008(fluid)

modelCO2 = GERG2008(["carbon dioxide"])
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
pbase = 27.0e5
#pbase = 5.25e5
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
        println("start braket...")
        println("ax = $(ax), fa = $(fa)")
        println("bx = $(bx), fb = $(fb)")
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
    return -T*T*ddg(T)
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
        return xs_min,z,f,fnc(xs_min)
    end
end

println("Check dry ice formation")
println("")

Tsolid = -78.75+273.15
psolid = 1.01325e5

βCO2, z_fluid, flash_fluid, gibbs_min = pT_flash_CO2Solid(model,psolid,Tsolid,x)

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

βCO2, z_fluid, flash_fluid, gibbs_min = pT_flash_CO2Solid(model,psf,Tsf,x)

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
    β, zf, fr, gs = pT_flash_CO2Solid(model,ps,Ts,zs)
    return (gs*Ts + s_spec*Ts/Clapeyron.R̄)
end

function psflashmin(T)
   Q = Qs(model,ps,T,zs,s_spec)
   return -Q
end

function Qs2(model,ps,Ts,zs,s_spec)
    β, zf, fr, gs = pT_flash_CO2Solid(model,ps,Ts,zs)
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

βCO2, z_fluid, flash_fluid, gibbs_min = pT_flash_CO2Solid(model,psf,Tsf,x)

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
mu_ambient = SuperTRAPP_mu(modelAir, T_ambient, rho_ambient, zAir)
lambda_ambient = SuperTRAPP_lambda(modelAir, T_ambient, rho_ambient, zAir)
cp_ambient = Clapeyron.VT_isobaric_heat_capacity(modelAir,vtf,T_ambient,zAir)/mwAir # J/Kg/K

T_ambient -= 273.15
println("T_ambient = $(T_ambient)")
println("mu_ambient = $(mu_ambient)")
println("lambda_ambient = $(lambda_ambient)")
println("cp_ambient = $(cp_ambient)")

zCO2 = 0.997275
zs = [zCO2,1.0-zCO2]

tank_diamter = 0.1937 #m
tank_radius = tank_diamter/2.0
tank_xsarea = pi*tank_radius^2
tank_length_TT = 76.80/4 #m
# spherical ends
tank_volume = tank_xsarea*tank_length_TT + 4/3*pi*tank_radius^3
tank_length = tank_volume/tank_xsarea

println("tank volume = $(tank_volume)")
println("tank length = $(tank_length)")

level = zeros(2)
level[1] = 0.05
#level[1] = 1.0 - 0.041/611.11
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

function orificeFlow(model,pin,Tin,zin,pcrit,Tcrit,pout,Tout,d0,d1,Cdliq,Cdvap,verbose)

    S0 = pi/4*d0^2
    S1 = pi/4*d1^2

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
        if length(βss) > 1
            vf = βss[1]*vss[1]/vtss
            Cd = Cdliq + vf*(Cdvap - Cdliq)
        else
            phase_type = Clapeyron.VT_identify_phase(model, vtss, Ts, zin)
            if phase_type == :liquid
                Cd = Cdliq
            elseif phase_type == :vapour
                Cd = Cdvap
            else
                Cd = 0.5*(Cdliq + Cdvap)
            end
        end
        flow = Cd*S0*rhoss*sign(Δh/mwtin)*sqrt(2.0*abs(Δh/mwtin)/(1.0 - (rhoss/rhoin*S0/S1)^2))
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
        if verbose
            println("bracket complete...")
            println("flowa = $(-Qa) mol/s, pa = $(pa/1e5) bara")
            println("flowc = $(-Qc) mol/s, pc = $(pc/1e5) bara")
        end
    #    pa = pin
    #    pc = pout
        bmres = brentmin(flowcalc,pa,pc,verbose)
        if verbose
            println("critical flow found ...")
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

#=

pin = ps
Tin = Ts
zin = zs
pout = 1.01325e5
Tcrit = Tin
Tout = Tin
pcrit = (pin + pout)/2.0

flash = Clapeyron.tp_flash_impl(model,pin,Tin,zin, HELDTPFlash(verbose = false))

xyin = flash.compositions

flow1,pcrit1,Tcrit1,Tout1 = orificeFlow(model,pin,Tin,xyin[1],pcrit,Tcrit,pout,Tout,2*25.4/1000.0,8.0*25.4/1000.0,0.67,0.98,true)

println("orifice flow1 = $(flow1) mol/s")
println("orifice critical pressure = $(pcrit1/1e5) bara")
println("orifice critical Temperature = $(Tcrit1 - 273.15) deg C")
println("orifice critical ratio pf/pin = $(pcrit1/pin)")
println("orifice critical pressure = $(pout/1e5) bara")
println("orifice critical Temperature = $(Tout1 - 273.15) deg C")

flow2,pcrit2,Tcrit2,Tout2 = orificeFlow(model,pin,Tin,xyin[2],pcrit,Tcrit,pout,Tout,2*25.4/1000.0,8.0*25.4/1000.0,0.67,0.98,true)

println("orifice flow2 = $(flow2) mol/s")
println("orifice critical pressure = $(pcrit2/1e5) bara")
println("orifice critical Temperature = $(Tcrit2 - 273.15) deg C")
println("orifice critical ratio pf/pin = $(pcrit2/pin)")
println("orifice critical pressure = $(pout/1e5) bara")
println("orifice critical Temperature = $(Tout2 - 273.15) deg C")

=#

function BlowDown(model, ps, Ts, zs, tank_volume, tank_radius, tank_length, level, pe, Δtmax, nstep, nVTK)

#    println("zs = $(zs)")
    
    mix_Tbub = bubble_temperature(model,ps,zs)
    Tbub = mix_Tbub[1]
    mix_Tdew = dew_temperature(model,ps,zs)
    Tdew = mix_Tdew[1]

    holdups = Vector{holdup}(undef,0)
    tank_xsa_holdups = Vector{Float64}(undef,0)

    v_solid = volume_dryice(ptp,Ttp,Ttp,hltp,sltp,gtp0,stp0,1,1)
    s_solid = entropy_dryice(ptp,Ttp,Ttp,hltp,sltp,gtp0,stp0,1,1)
    g_solid = gibbs_dryice(ptp,Ttp,Ttp,hltp,sltp,gtp0,stp0,1,1)
    h_solid = g_solid + Ttp*s_solid
    phases = Vector{phase}(undef,0)
    push!(phases, phase(1.0,0.0,0.0,mwCO2,v_solid,h_solid,s_solid,[1.0,0.0]))
    holdup_solid = holdup(ptp,Ttp,0.0,0.0,mwCO2,v_solid,h_solid,s_solid,[1.0,0.0],phases)

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
        push!(holdups, holdup(ps,Ts,tv,tm,mw,v,h,s,ybub,phases))

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
        push!(holdups, holdup(ps,Ts,tv,tm,mwtf,vtf,htf,stf,xbub,phases))

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
        push!(holdups, holdup(ps,Ts,tv,tm,mwtf,vtf,htf,stf,ydew,phases))

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
        push!(holdups, holdup(ps,Tdew,tv,tm,mw,v,h,s,xdew,phases))

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
            push!(holdups, holdup(ps,Ts,level[ip]*tank_volume,level[ip]*tank_volume/vf[ip],mwf[ip],vf[ip],hf[ip],sf[ip],xf[ip],phases))
            push!(tank_xsa_holdups, level[ip]*tank_volume/tank_length)
        end

        println("Two phase holdups = $(holdups)")

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

    Tcrit_vap = holdups[1].T
    Tout_vap = holdups[1].T
    pout_vap = 0.995e5
    pcrit_vap = (ps + pout_vap)/2.0

    Tcrit_liq = holdups[1].T
    Tout_liq = holdups[1].T
    pout_liq = 0.995e5
    pcrit_liq = (ps + pout_liq)/2.0

    println("")
    println("Blowndown Function Start:")
    println("Initialization:")
    println("")
#    println("holdups = $(holdups)")
    
    Δpbase = (ps - pe)/nstep
    ΔT = fill(-0.33,length(holdups))
    ΔQ = zeros(length(holdups))
    ΔQ_solid = 0.0

    time_Plot = Vector{Float64}(undef,0)
    Twallmin_Plot = Vector{Float64}(undef,0)
    Twallmin_av_Plot = Vector{Float64}(undef,0)
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

    fv_Plot = Vector{Float64}(undef,0)
    fl_Plot = Vector{Float64}(undef,0)

    push!(Twallmin_Plot, Twall_initial - 273.15)
    push!(Twallmin_av_Plot, Twall_initial - 273.15)
    push!(Twallnoz_Plot, Twall_initial - 273.15)
    push!(Twall_MMDT_Plot, -46.0)
    push!(T1_Plot, holdups[1].T-273.15)
    push!(T2_Plot, holdups[2].T-273.15)

    push!(p1_Plot, holdups[1].p/1e5)
    push!(p2_Plot, holdups[2].p/1e5)
    push!(ps_Plot, holdup_solid.p/1e5)

    push!(m1_Plot, holdups[1].moles*holdups[1].mw/1000)
    push!(m2_Plot, holdups[2].moles*holdups[2].mw/1000)
    push!(ms_Plot, holdup_solid.moles*holdup_solid.mw/1000)

    push!(fv_Plot, 0.0)
    push!(fl_Plot, 0.0)

    time = 0.0
    push!(time_Plot, time)

    # Set up FEM wall

    # heat source
    heat_source = 0.0

    #properties
    rho_steel = 7850.0
    cp_steel = 490.0
    lambda_steel = 45.0
#    alpha_steel = lambda_steel/rho_steel/cp_steel
    rho_insulation = 30.0
    cp_insulation = 1500.0
    lambda_insulation = 0.025
#    alpha_insulation = lambda_insulation/rho_insulation/cp_insulation

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
    Twall_noz = Twall_initial - 273.15
    Twall_outer = Twall_initial - 273.15
    #T_ambient = holdups[1].T - 273.15
    wall_thickness = 0.013
    insulation_thickness = 0.015
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
    nradius = 28
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

    uₙ = zeros(ndofs(dh));
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

    # assume vapour removal
    topflow = true
#    topflow = false

    level_min = 0.01
    alpha_lag = 0.5
    volume_fraction_noz = 1.0

    
    Cdliq = 0.67
    Cdvap = 0.98
#    d1 = 2.0*25.4/1000.0/sqrt(2)
    d1 = 5/1000.0
    d2 = 8.0*25.4/1000.0
    number_tanks_parallel = 1
    println("number of tanks parallel = $(number_tanks_parallel)")

    flowout_vap1, pcrit_vap, Tcrit_vap, Tout_vap = orificeFlow(model,holdups[1].p,holdups[1].T,holdups[1].z,pcrit_vap,Tcrit_vap,pout_vap,Tout_vap,d1,d2,Cdliq,Cdvap,false)
    flowout_liq1, pcrit_liq, Tcrit_liq, Tout_liq = orificeFlow(model,holdups[2].p,holdups[2].T,holdups[2].z,pcrit_liq,Tcrit_liq,pout_liq,Tout_liq,d1,d2,Cdliq,Cdvap,false)

    flowout_vap1 /= number_tanks_parallel
    flowout_liq1 /= number_tanks_parallel

    pcrit_ratio_vap = pcrit_vap/holdups[1].p
    pcrit_ratio_liq = pcrit_liq/holdups[2].p

    println("flowout vapour = $(flowout_vap1*holdups[1].mw) kg/s")
    println("flowout liquid = $(flowout_liq1*holdups[2].mw) kg/s")
        
    println("flowout vapour = $(flowout_vap1*holdups[1].v) m3/s")
    println("flowout liquid = $(flowout_liq1*holdups[2].v) m3/s")

    istep = 0
    iVTK = 0
    Δt = Δtmax
    Δp = Δpbase

    while (ps >= pe)

        istep += 1
        iVTK += 1

        println("")
        println("Step: $(istep)")
        
        moles_solid = 0.0
        volume_solid = 0.0
        hm_solid = 0.0

        #=    
            if ps >= 4.75e5 && ps <= 5.5e5
               if Δt > Δtmax/5
                    Δp = Δp*Δtmax/5/Δt
                else
                    Δp = 0.5*Δpbase + 0.5*Δp
                end
            else
                if Δt > Δtmax
                    Δp = Δp*Δtmax/Δt
                else
                    Δp = 0.5*Δpbase + 0.5*Δp
                end
            end
        =#

        for ih = eachindex(holdups)

            println("")
            println("Δp = $(Δp)")

            ps = holdups[ih].p - Δp
        #    println("holdup[$(ih)]")

        #    if holdups[ih].moles > 0.0

                ΔsLast = 0.0
                ΔsError = 1
                println("entropy iteration")
                is = 0
                Ts = holdups[ih].T + ΔT[ih]
                while(ΔsError > 1.0e-6)
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
                    #    fr = Clapeyron.tp_flash_impl(model,p,T,x, HELDTPFlash(verbose = false))
                    #    g = fr.data.g
                        β, zf, fr, g = pT_flash_CO2Solid(model,p,T,x)
                        return (g*T + sspec*T/Clapeyron.R̄)
                    end

                    function psflashmin(T)
                        Q = Qs(model,ps,T,zs,ss)
                        return -Q
                    end

                    Ta = Ts + 2.0*ΔT[ih]
                    Tc = Ts - 2.0*ΔT[ih]

                    Ta, Tc, Qa, Qc = bracketmin(psflashmin,Ta,Tc,5.0,false)

                    bmres = brentmin(psflashmin, Ta, Tc, false)
                    ΔT[ih] = bmres[1] - holdups[ih].T
                    Ts = bmres[1]

                end
                
                moles_total = holdups[ih].moles

                if moles_total > 0.0

                    holdups[ih].p = ps
                    holdups[ih].T = Ts
                #    flash = Clapeyron.tp_flash_impl(model,holdups[ih].p,holdups[ih].T,holdups[ih].z, HELDTPFlash(verbose = false))
                    β_solid, z_fluid, flash, g_total = pT_flash_CO2Solid(model,holdups[ih].p,holdups[ih].T,holdups[ih].z)
                    # if we form a solid then the composition of the fluid changes.
                    holdups[ih].z = z_fluid
                    βs,mws,mwts,vs,vts,hs,hts,ss,sts,xs  = get_props(model,flash)
                    holdups[ih].mw = mwts
                    holdups[ih].v  = vts
                    holdups[ih].h  = hts
                    holdups[ih].s  = sts
                    moles_fluid = (1.0 - β_solid)*moles_total
                    holdups[ih].moles = moles_fluid
                    holdups[ih].volume  = moles_fluid*holdups[ih].v
                    phases = Vector{phase}(undef,0)
                    for ip = eachindex(βs)      
                        push!(phases, phase(βs[ip],βs[ip]*holdups[ih].moles*vs[ip],βs[ip]*holdups[ih].moles,mws[ip],vs[ip],hs[ip],ss[ip],xs[ip]))
                    end
                    holdups[ih].phases = phases

                    if β_solid > 0.0
                        println("dry ice detected in phase[$(ih)] = $(β_solid)")
                        v_solid = volume_dryice(ps,Ts,Ttp,hltp,sltp,gtp0,stp0,1,1)
                        s_solid = entropy_dryice(ps,Ts,Ttp,hltp,sltp,gtp0,stp0,1,1)
                        g_solid = gibbs_dryice(ps,Ts,Ts,hltp,sltp,gtp0,stp0,1,1)
                        h_solid = g_solid + Ts*s_solid
                        moles_solid += β_solid*moles_total
                        volume_solid += β_solid*moles_total*v_solid
                        hm_solid += β_solid*moles_total*h_solid
                    #    println("moles_solid $(ih)] = $(β_solid*moles_total)")
                    #    println("h_solid $(ih)] = $(h_solid)")
                    else
                        println("no dry ice detected in phase[$(ih)] = $(β_solid)")
                    end

                else

                    holdups[ih].p = ps
                    holdups[ih].T = Ts
                    flash = Clapeyron.tp_flash_impl(model,holdups[ih].p,holdups[ih].T,holdups[ih].z, HELDTPFlash(verbose = false))
                    βs,mws,mwts,vs,vts,hs,hts,ss,sts,xs  = get_props(model,flash)
                    holdups[ih].mw = mwts
                    holdups[ih].v  = vts
                    holdups[ih].h  = hts
                    holdups[ih].s  = sts
                    moles_fluid = moles_total
                    holdups[ih].moles = moles_fluid
                    holdups[ih].volume  = moles_fluid*holdups[ih].v
                    phases = Vector{phase}(undef,0)
                    for ip = eachindex(βs)      
                        push!(phases, phase(βs[ip],βs[ip]*holdups[ih].moles*vs[ip],βs[ip]*holdups[ih].moles,mws[ip],vs[ip],hs[ip],ss[ip],xs[ip]))
                    end
                    holdups[ih].phases = phases
                    
                end

        #    end

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
        #    println("holdups[$(ih)] = $(holdups[ih])")

        end

        # we should have liquid formed for condensation
        # condensation also requires a cold wall Tvap - Tw > 0.0
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

        # we should have vapour formed for boiling
        # boiling also requires a hot wall Tw - Tliq > 0.0
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

        
        # add each phases solid to the overall solid phase
    #    println("h_initial = $(holdup_solid.h)")
    #    println("moles_initial = $(holdup_solid.moles)")
    #    println("moles_solid = $(moles_solid)")
    #    println("hm_solid = $(hm_solid)")
    #    println("ΔQ_solid = $(ΔQ_solid)")
        moles_initial = holdup_solid.moles 
        holdup_solid.moles += moles_solid
        holdup_solid.volume += volume_solid
        if holdup_solid.moles > 0.0
            holdup_solid.h = (holdup_solid.h*moles_initial + hm_solid)/holdup_solid.moles + ΔQ_solid
        end
        holdup_solid.p = ps

        # we need a phflash to find mixing temperature for now assume just avergae of fluids
        # the flash may produce a fluid this needs to be added to the fluid holdups

        Tsv = (holdups[1].T + holdups[2].T)/2.0

        βsolid = 1.0
        vs = 0.0
        hs = 0.0
        ss = 0.0
        vv = 0.0
        hv = 0.0
        sv = 0.0

        if holdup_solid.moles > 0.0

        # This assumes the solid form below the triple point so we have solid/vapour equilibria
        # need to check vv = volume(modelCO2,ps,T,[1]) retorns the vapour root.

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
        #    println("Gibbs Tsv = $(Tsv), error = $(error)")
            dTsv = 0.001
            while (error > sqrt(eps(Float64)))
                df_dT = (fnc_gibbs(Tsv+dTsv) - f)/dTsv
                Tsv = Tsv - f/df_dT
                f = fnc_gibbs(Tsv)
                error = abs(f)
        #        println("Gibbs Tsv = $(Tsv), error = $(error)")
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

        #    println("vv = $(vv) m3/mol")
        #    println("gv = $(gv) J/mol")
        #    println("μv = $(μv) J/mol")
        #    println("gs = $(gs) J/mol")

        #    println("hv = $(hv) J/mol")
        #    println("hs = $(hs) J/mol")

        #    println("holdup_solid.h = $(holdup_solid.h) J/mol")

            if holdup_solid.h > hv

                βsolid = 0.0

                function fnc_hv(T)
                    vv = volume(modelCO2,ps,T,[1])
                    hv = Clapeyron.VT_enthalpy(modelCO2, vv, T, [1])
                    return hv - holdup_solid.h
                end

                f = fnc_hv(Tsv)
                error = abs(f)
            #    println("hv Tsv = $(Tsv), error = $(error)")
                dTsv = 0.001
                while (error > sqrt(eps(Float64)))
                    df_dT = (fnfnc_hvc(Tsv+dTsv) - f)/dTsv
                    Tsv = Tsv - f/df_dT
                    f = fnc_hv(Tsv)
                    error = abs(f)
                #    println("hv Tsv = $(Tsv), error = $(error)")
                end
                println("Superheated vapour")

                vv = volume(modelCO2,ps,Tsv,[1])
                hv = Clapeyron.VT_enthalpy(modelCO2, vv, Tsv, [1])
                sv = Clapeyron.VT_entropy(modelCO2,vv,Tsv,[1])

            elseif holdup_solid.h < hs

                βsolid = 1.0

                function fnc_hs(T)
                    gs = gibbs_dryice(ps,T,Ttp,hltp,sltp,gtp0,stp0,1,1)
                    ss = entropy_dryice(ps,T,Ttp,hltp,sltp,gtp0,stp0,1,1)
                    hs = gs + T*ss
                    return hs - holdup_solid.h
                end

                f = fnc_hs(Tsv)
                error = abs(f)
            #    println("hs Tsv = $(Tsv), error = $(error)")
                dTsv = 0.001
                while (error > sqrt(eps(Float64)))
                    df_dT = (fnc_hs(Tsv+dTsv) - f)/dTsv
                    Tsv = Tsv - f/df_dT
                    f = fnc_hs(Tsv)
                    error = abs(f)
                #    println("hs Tsv = $(Tsv), error = $(error)")
                end
                println("Subcooled solid")

                gs = gibbs_dryice(ps,Tsv,Ttp,hltp,sltp,gtp0,stp0,1,1)
                ss = entropy_dryice(ps,Tsv,Ttp,hltp,sltp,gtp0,stp0,1,1)
                hs = gs + Tsv*ss
                vs = volume_dryice(ps,Tsv,Ttp,hltp,sltp,gtp0,stp0,1,1)

            else

                βsolid = (holdup_solid.h - hv)/(hs - hv)
                println("Equilibrium solid/vapour")

            end

            println("βsolid = $(βsolid)")

            holdup_solid.T = Tsv


            # once we know βsolid we can separate the vapour formed and add to holdup[1]
            # do this here

            if βsolid < 1.0

                moles_total = holdup_solid.moles
                cm = Vector{Float64}(undef,0)
                mt1 = holdups[1].moles
                β1 = holdups[1].phases[1].β
                z1 = holdups[1].phases[1].z
                h1 = holdups[1].phases[1].h
                mtsv = (1.0 - βsolid)*moles_total
                zsv = holdup_solid.z
                for ic = eachindex(holdups[1].z)
                        push!(cm,mt1*z1[ic] + mtsv*zsv[ic])
                end
                h1m = mt1*h1 + mtsv*hv
                mt1 = sum(cm)
                if mt1 > 0.0
                    for ic = eachindex(cm)
                        z1[ic] = cm[ic]/mt1
                    end
                    h1 = h1m/mt1
                    holdups[1].moles = mt1
                    holdups[1].h = h1
                    holdups[1].z = z1
                end

                holdup_solid.moles = βsolid*moles_total
                holdup_solid.volume = holdup_solid.moles*vs
                holdup_solid.v = vs
                holdup_solid.h = hs
                holdup_solid.s = ss

                phases = Vector{phase}(undef,0)
                push!(phases, phase(1.0,holdup_solid.volume,holdup_solid.moles,mwCO2,holdup_solid.v,holdup_solid.h,holdup_solid.s,[1.0,0.0]))
                holdup_solid.phases = phases

            else

                holdup_solid.volume = holdup_solid.moles*vs
                holdup_solid.v = vs
                holdup_solid.h = hs
                holdup_solid.s = ss

                phases = Vector{phase}(undef,0)
                push!(phases, phase(1.0,holdup_solid.volume,holdup_solid.moles,mwCO2,holdup_solid.v,holdup_solid.h,holdup_solid.s,[1.0,0.0]))
                holdup_solid.phases = phases

            end

        end

        println("holdup_solid = $(holdup_solid)")

        # recombine vapour with liquid phase vapour and liquid with vapour phase liquid
        cm = Vector{Float64}(undef,0)
        mt1 = holdups[1].moles
        β1 = holdups[1].phases[1].β
        z1 = holdups[1].phases[1].z
        h1 = holdups[1].phases[1].h
        if length(holdups[2].phases) > 1
            mt2 = holdups[2].moles
            β2 = holdups[2].phases[1].β
            z2 = holdups[2].phases[1].z
            h2 = holdups[2].phases[1].h
            for ic = eachindex(holdups[1].z)
                push!(cm,β1*mt1*z1[ic] + β2*mt2*z2[ic])
            end
            h1m = β1*mt1*h1 + β2*mt2*h2
        else
            phase_type = Clapeyron.VT_identify_phase(model, holdups[2].v, holdups[2].T, holdups[2].z)
            if phase_type == :vapour
                mt2 = holdups[2].moles
                β2 = holdups[2].phases[1].β
                z2 = holdups[2].phases[1].z
                h2 = holdups[2].phases[1].h
                for ic = eachindex(holdups[1].z)
                    push!(cm,β1*mt1*z1[ic] + β2*mt2*z2[ic])
                end
                h1m = β1*mt1*h1 + β2*mt2*h2
            else
                for ic = eachindex(holdups[1].z)
                    push!(cm,β1*mt1*z1[ic])
                end
                h1m = β1*mt1*h1
            end
        end
        mt1 = sum(cm)
        if mt1 > 0.0
        #    z1 = Vector{Float64}(undef,0)
            for ic = eachindex(cm)
            #    push!(z1,cm[ic]/mt1)
                z1[ic] = cm[ic]/mt1
            end
            h1 = h1m/mt1
            holdups[1].moles = mt1
            holdups[1].h = h1
            holdups[1].z = z1
        #    println("cm1 = $(cm) moles")
        #    println("mt1 = $(mt1) moles")
        #    println("z1 = $(z1)")
        #    println("h1 = $(h1) J/mol")
        else
            holdups[1].moles = mt1
        end

        cm = Vector{Float64}(undef,0)
        mt2 = holdups[2].moles
        β2 = holdups[2].phases[1].β
        z2 = holdups[2].phases[1].z
        h2 = holdups[2].phases[1].h
        if length(holdups[2].phases) > 1
            β2 = holdups[2].phases[2].β
            z2 = holdups[2].phases[2].z
            h2 = holdups[2].phases[2].h
        end
        if length(holdups[1].phases) > 1
            mt1 = holdups[1].moles
            β1 = holdups[1].phases[2].β
            z1 = holdups[1].phases[2].z
            h1 = holdups[1].phases[2].h
            for ic = eachindex(holdups[2].z)
                push!(cm,β2*mt2*z2[ic] + β1*mt1*z1[ic])
            end
            h2m = β2*mt2*h2 + β1*mt1*h1
        else
            phase_type = Clapeyron.VT_identify_phase(model, holdups[2].v, holdups[2].T, holdups[2].z)
            if phase_type == :vapour
                for ic = eachindex(holdups[2].z)
                    push!(cm,0.0)
                end
            else
                for ic = eachindex(holdups[2].z)
                    push!(cm,β2*mt2*z2[ic])
                end
                h2m = β2*mt2*h2
            end
        end
        mt2 = sum(cm)
        if mt2 > 0.0
        #    z2 = Vector{Float64}(undef,0)
            for ic = eachindex(cm)
            #    push!(z2,cm[ic]/mt2)
                z2[ic] = cm[ic]/mt2
            end
            h2 = h2m/mt2
            holdups[2].moles = mt2
            holdups[2].h = h2
            holdups[2].z = z2
        #    println("cm2 = $(cm) moles")
        #    println("mt2 = $(mt2) moles")
        #    println("z2 = $(z2)")
        #    println("h2 = $(h2) J/mol")
        else
            holdups[2].moles = mt2
        end

        # Ok we have all we need to do a phflash and reset the holpups with the new conditions after separation and mixing

                for ih = eachindex(holdups)

                    hs = holdups[ih].h
                    Ts = holdups[ih].T
                    ps = holdups[ih].p
                    zs = holdups[ih].z

                    function Qh(model,p,T,x,hspec)
                        fr = Clapeyron.tp_flash_impl(model,p,T,x, HELDTPFlash(verbose = false))
                    #    println("fr = $(fr)")
                        g = fr.data.g
                        return (g - hspec/Clapeyron.R̄/T)
                    end

                #    function phflashmin(T)
                #        Q,β = Qh(model,ps,T,zs,hs)
                #        return -Q,length(β)
                #    end

                    function phflashmin(T)
                        Q = Qh(model,ps,T,zs,hs)
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

        tank_volume_new = holdups[1].volume + holdups[2].volume + holdup_solid.volume
        Δv = tank_volume_new - tank_volume
    #    ΔL = Δv/(pi*tank_radius^2)
    #    tank_length_new = tank_length + ΔL

    #    xsa_vap = holdups[1].volume/tank_length_new
    #    xsa_liq = holdups[2].volume/tank_length_new
        
    #    Δv_vap_min = xsa_vap*ΔL
    #    Δv_liq_min = xsa_liq*ΔL

        # minimum level that all vapour or liquid phase is removed in a time step.
        # this ramps down to zero at zero liquid level
        # represents the inefficiency of remove a phase at low low or high high levels

    #    println("Δv_vap_min = $(Δv_vap_min)")
    #    println("Δv_liq_min = $(Δv_liq_min)")
    #    println("Δv_vap_min + Δv_liq_min = $(Δv_vap_min+Δv_liq_min)")
        println("Δv = $(Δv)")
    #    holdups[1].volume -= Δv
    #    Δmoles = holdups[1].moles - holdups[1].volume/holdups[1].v
    #    holdups[1].moles = holdups[1].volume/holdups[1].v

        pcrit_vap = pcrit_ratio_vap*holdups[1].p
        flowout_vap2, pcrit_vap, Tcrit_vap, Tout_vap = orificeFlow(model,holdups[1].p,holdups[1].T,holdups[1].z,pcrit_vap,Tcrit_vap,pout_vap,Tout_vap,d1,d2,Cdliq,Cdvap,false)
        
        pcrit_liq = pcrit_ratio_liq*holdups[2].p
        flowout_liq2, pcrit_liq, Tcrit_liq, Tout_liq = orificeFlow(model,holdups[2].p,holdups[2].T,holdups[2].z,pcrit_liq,Tcrit_liq,pout_liq,Tout_liq,d1,d2,Cdliq,Cdvap,false)

        flowout_vap2 /= number_tanks_parallel
        flowout_liq2 /= number_tanks_parallel
        
        pcrit_ratio_vap = pcrit_vap/holdups[1].p
        pcrit_ratio_liq = pcrit_liq/holdups[2].p

        println("pcrit_ratio vap = $(pcrit_ratio_vap)")
        println("pcrit_ratio liq = $(pcrit_ratio_liq)")

        # use trapazodial rule for flows and integration for time
        flowout_vap = 0.5*(flowout_vap1 + flowout_vap2)
        flowout_liq = 0.5*(flowout_liq1 + flowout_liq2)

        flowout_vap1 = flowout_vap2
        flowout_liq1 = flowout_liq2

        println("flowout vapour = $(flowout_vap*holdups[1].mw) kg/s")
        println("flowout liquid = $(flowout_liq*holdups[2].mw) kg/s")
        
        println("flowout vapour = $(flowout_vap*holdups[1].v) m3/s")
        println("flowout liquid = $(flowout_liq*holdups[2].v) m3/s")

        if topflow
            level_noz = levelPercent/100.0
            volume_fraction_noz = alpha_lag*max(min((1.0 - level_noz)/level_min,1.0),0.0) + (1.0 - alpha_lag)*volume_fraction_noz
            println("volume_fraction_noz = $(volume_fraction_noz)")
            Δv_vap = min(Δv,holdups[1].volume*volume_fraction_noz)
            Δv_liq = Δv - Δv_vap
            Δt_liq = Δv_liq/(flowout_liq*holdups[2].v)
            Δt_vap = Δv_vap/(flowout_vap*holdups[1].v)
            holdups[2].volume -= Δv_liq
            Δmoles_liq = max(holdups[2].moles - holdups[2].volume/holdups[2].v,0.0)
            holdups[2].moles -= Δmoles_liq
            holdups[1].volume -= Δv_vap
            Δmoles_vap = max(holdups[1].moles - holdups[1].volume/holdups[1].v,0.0)
            holdups[1].moles -= Δmoles_vap
        else
            level_noz = (levelPercent - solidPercent)/100.0
            volume_fraction_noz = alpha_lag*max(min(level_noz/level_min,1.0),0.0) + (1.0 - alpha_lag)*volume_fraction_noz
            println("volume_fraction_noz = $(volume_fraction_noz)")
            Δv_liq = min(Δv,holdups[2].volume*volume_fraction_noz)
            Δv_vap = Δv - Δv_liq
            Δt_liq = Δv_liq/(flowout_liq*holdups[2].v)
            Δt_vap = Δv_vap/(flowout_vap*holdups[1].v)
            holdups[2].volume -= Δv_liq
            Δmoles_liq = max(holdups[2].moles - holdups[2].volume/holdups[2].v,0.0)
            holdups[2].moles -= Δmoles_liq
            holdups[1].volume -= Δv_vap
            Δmoles_vap = max(holdups[1].moles - holdups[1].volume/holdups[1].v,0.0)
            holdups[1].moles -= Δmoles_vap
        end

        println("holdups[$(1)] = $(holdups[1])")
        println("holdups[$(2)] = $(holdups[2])")

        println("holdup volume = $(holdups[1].volume + holdups[2].volume + holdup_solid.volume)")

        println("Δmoles_vap = $(Δmoles_vap)")
        println("Δmoles_liq = $(Δmoles_liq)")

    #    println("orifice flow = $(flowout) mol/s")
    #    println("orifice critical pressure = $(pcrit/1e5) bara")
    #    println("orifice isentropic critical Temperature = $(Tcrit - 273.15) deg C")
    #    println("orifice critical ratio pf/pin = $(pcrit/pin)")
    #    println("orifice outlet pressure = $(pout/1e5) bara")
    #    println("orifice isentropic outlet Temperature = $(Tout - 273.15) deg C")

        # to find actual valve outlet condition we need an isenthalpy flash from pin to pout.
        
    #    Δt_vap = Δmoles_vap/flowout_vap
    #    Δt_liq = Δmoles_liq/flowout_liq

        if Δmoles_vap > 0.0
            fv_plot = flowout_vap*holdups[1].mw
        else
            fv_plot = 0.0
        end

        if Δmoles_liq > 0.0
            fl_plot = flowout_liq*holdups[2].mw
        else
            fl_plot = 0.0
        end

        push!(fv_Plot, fv_plot)
        push!(fl_Plot, fl_plot)

        Δt = Δt_vap + Δt_liq

        println("Δt_vap = $(Δt_vap)")
        println("Δt_liq = $(Δt_liq)")

        time += Δt
        println("time = $(time/3600)")
        push!(time_Plot, time/3600.0)

        # new level
        for ih = eachindex(holdups)
            tank_xsa_holdups[ih] = holdups[ih].volume/tank_length
        end

        tank_xsa_holdup_solid = holdup_solid.volume/tank_length

        volumeFraction = tank_xsa_holdup_solid/(pi*tank_radius^2)

        if tank_xsa_holdup_solid > 10.0*eps(Float64)
             if theta_solid == 0.0
                theta_solid = sqrt(0.9*volumeFraction)*360/180*pi
            end
            error = 1.0
            while (error > 0.0001)
                f =  (tank_radius^2)*(theta_solid - sin(theta_solid)) - 2.0*tank_xsa_holdup_solid
                df_dtheta = (tank_radius^2)*(1.0 - cos(theta_solid))
                theta_solid = theta_solid - f/df_dtheta
                error = abs(f)
            #    println("theta_solid = $(theta_solid), error = $(error)")
            end
        else
            theta_solid = 0.0
        end
 
        # height based on tank diamter
        tank_diameter = 2.0*tank_radius
        solidPercent = (1.0 - cos(theta_solid/2.0))/2.0*100.0
        height_solid = tank_diameter*solidPercent/100.0
        println("Solid Percent = $(solidPercent)")
   
        cord_solid = 2.0*tank_radius*sin(theta_solid/2.0)
        interface_area_solid = cord_solid*tank_length

        volumeFraction = (tank_xsa_holdups[2] + tank_xsa_holdup_solid)/(pi*tank_radius^2)

        if tank_xsa_holdups[2] + tank_xsa_holdup_solid > 10.0*eps(Float64)
            if theta == 0.0
                theta = sqrt(0.9*volumeFraction)*360/180*pi
            end
            error = 1.0
            while (error > 0.0001)
                f =  (tank_radius^2)*(theta - sin(theta)) - 2.0*(tank_xsa_holdups[2] + tank_xsa_holdup_solid)
                df_dtheta = (tank_radius^2)*(1.0 - cos(theta))
                theta = theta - f/df_dtheta
                error = abs(f)
            #    println("theta = $(theta), error = $(error)")
            end
        else
            theta = 0.0
        end

        # height based on tank diamter
        tank_diameter = 2.0*tank_radius
        levelPercent = (1.0 - cos(theta/2.0))/2.0*100.0
        height[2] = tank_diameter*levelPercent/100.0
        height[1] = tank_diameter*(1.0 - levelPercent/100.0)
        println("Level Percent = $(levelPercent)")

        cord = 2.0*tank_radius*sin(theta/2.0)
        interface_area = cord*tank_length
        #println("interface_area = $(interface_area) m2")

        # used for position along wall as an arc length used to define position of interface
        #arc_length = pi*tank_radius
        #arc_length_liquid = theta/2.0*tank_radius
        #println("arc_length_liquid/arc_length = $(arc_length_liquid/arc_length)")

        # all areas are based on total area for calculating heat from fluids. Only 50% goes to wall as this is 
        # modelled symmetrically
        #wall_vapour_area = 2.0*(arc_length - arc_length_liquid)*tank_length
        #println("wall_vapour_area = $(wall_vapour_area) m2")

        #wall_liquid_area = 2.0*arc_length_liquid*tank_length
        #println("wall_liquid_area = $(wall_liquid_area) m2")

        # need htc for vapour wall and liquid wall based on transport properties, Re, Pr, Gr etc

        dia_min = 0.010
        dia_vap = max(sqrt(4.0/pi*tank_xsa_holdups[1]), dia_min)
        if tank_xsa_holdups[1] > pi/4.0*dia_min*dia_min
            vel_vap = (Δmoles_vap/Δt)*holdups[1].v/tank_xsa_holdups[1]
        else
            vel_vap = (Δmoles_vap/Δt)*holdups[1].v/(pi/4.0*dia_min*dia_min)
        end
        rho_vap = 1.0/holdups[1].v # m3/mol
        mu_vap = SuperTRAPP_mu(model,holdups[1].T,rho_vap,holdups[1].z) # Ps.s
        lambda_vap = SuperTRAPP_lambda(model,holdups[1].T,rho_vap,holdups[1].z) # W/m/K
        cp_vap = Clapeyron.VT_isobaric_heat_capacity(model,holdups[1].v,holdups[1].T,holdups[1].z)/holdups[1].mw # J/Kg/K

    #    println("vel_vap = $(vel_vap) m/s")
    #    println("mu_vap = $(mu_vap) Pa.s")
    #    println("lambda_vap = $(lambda_vap) W/m/K")
    #    println("cp_vap = $(cp_vap) J/Kg/K")

        function alpha_vap(Tw, T, u, d, rho, mu, cp , lambda)
        
            g = 9.8065
            Pr = cp*mu/lambda
            Re = rho*u*d/mu

        #    println("Pr = $(Pr)")
        #    println("Re = $(Re)")

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

        #    println("Gr = $(Gr)")

            NuNC = 0.825 + 0.387*(Gr*Pr / (1 + (0.492/Pr)^(9/16))^(16/9))^(1/6)
            NuNC = NuNC * NuNC

            Nu = sqrt(NuNC * NuNC + NuLam * NuLam + NuTurb * NuTurb)

           return lambda*Nu/d

        end

        htc_vap = alpha_vap(Twall_vap+273.15, holdups[1].T, vel_vap, dia_vap, rho_vap*holdups[1].mw, mu_vap, cp_vap ,lambda_vap)

        if holdups[1].volume < pi/4.0*dia_min*dia_min*tank_length
            htc_vap = 0.0
        end
        
        println("htc_vap = $(htc_vap)")

        dia_liq = max(sqrt(4.0/pi*tank_xsa_holdups[2]), dia_min)
        if tank_xsa_holdups[2] > pi/4.0*dia_min*dia_min
            vel_liq = (Δmoles_liq/Δt)*holdups[2].v/tank_xsa_holdups[2]
        else
            vel_liq = (Δmoles_liq/Δt)*holdups[2].v/(pi/4.0*dia_min*dia_min)
        end
        rho_liq = 1.0/holdups[2].v
        mu_liq = SuperTRAPP_mu(model,holdups[2].T,rho_liq,holdups[2].z)
        lambda_liq = SuperTRAPP_lambda(model,holdups[2].T,rho_liq,holdups[2].z)
        cp_liq = Clapeyron.VT_isobaric_heat_capacity(model,holdups[2].v,holdups[2].T,holdups[2].z)/holdups[2].mw

    #    println("vel_liq = $(vel_liq) m/s")
    #    println("mu_liq = $(mu_liq) Pa.s")
    #    println("lambda_liq = $(lambda_liq) W/m/K")
    #    println("cp_liq = $(cp_liq) J/Kg/K")

        function alpha_liq(Tw, T, u, d, rho, mu, cp , lambda, BOILING)
            
            g = 9.8065
            Pr = cp*mu/lambda
            Re = rho*u*d/mu

        #    println("Pr = $(Pr)")
        #    println("Re = $(Re)")

            NuLam = 0.0
            NuTurb = 0.0
            if Re > 1.0

                NuLam = 3.657
                f = 1.82*log10(Re) - 1.64
                f = 1.0/f/f
                NuTurb = f/8.0*Re*Pr/(1.0 + 12.7*sqrt(f/8.0)*(Pr^(2/3) - 1.0))

            end

            deltaT0 = 20000.0 / 3500.0
            deltaT = Tw - T
            if deltaT < 0.0 
                deltaT = 0.0
            end

            NuBOIL = d / lambda * 3500.0 * BOILING * (deltaT / deltaT0)^0.25

            beta = 2.0/(Tw + T)
            Gr = d*d*d*g*rho*rho*beta*abs(Tw - T)/mu/mu

        #    println("Gr = $(Gr)")

            NuNC = 0.825 + 0.387*(Gr*Pr / (1 + (0.492/Pr)^(9/16))^(16/9))^(1 /6)
            NuNC = NuNC * NuNC

            Nu = sqrt(NuNC * NuNC + NuLam * NuLam + NuTurb * NuTurb + NuBOIL * NuBOIL)

            return lambda*Nu/d

        end

        htc_liq = alpha_liq(Twall_liq+273.15, holdups[2].T, vel_liq, dia_liq, rho_liq*holdups[2].mw, mu_liq, cp_liq ,lambda_liq, BOILING)
    
        if holdups[2].volume < pi/4.0*dia_min*dia_min*tank_length
            htc_liq = 0.0
        end

        println("htc_liq = $(htc_liq)")

        dia_solid = max(sqrt(4.0/pi*tank_xsa_holdup_solid), dia_min)
        Nu_solid = 1.0
        pratio = ps/200/1e6
        A0 =  12.1860*pratio*pratio - 75.579*pratio + 166.6
        A1 = -(9.2052*pratio*pratio - 46.508*pratio + 98.634)
        A2 =  2.2499*pratio*pratio - 9.5547*pratio + 19.565
        A3 = -(0.1788*pratio*pratio - 0.6566*pratio + 1.3047)
        lambda_solid = exp(A0 + A1*log(Twall_solid+273.15) + A2*log(Twall_solid+273.15)^2 + A3*log(Twall_solid+273.15)^3) # W/m/K
        println("lambda_solid = $(lambda_solid)")
        htc_solid = Nu_solid*lambda_solid/dia_solid

        if holdup_solid.volume < pi/4.0*dia_min*dia_min*tank_length
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

        #    println("Pr = $(Pr)")
        #    println("Re = $(Re)")

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

        #    println("Gr = $(Gr)")

            NuNC = 0.825 + 0.387*(Gr*Pr / (1 + (0.492/Pr)^(9/16))^(16/9))^(1/6)
            NuNC = NuNC * NuNC

            Nu = sqrt(NuNC * NuNC + NuLam * NuLam + NuTurb * NuTurb)

           return lambda*Nu/d

        end

        if Δmoles_vap > 0.0
            vel_noz_vap = flowout_vap*holdups[1].v/(pi/4*d2^2)
        else
            vel_noz_vap = 0.0
        end
        htc_noz_vap = alpha_noz(Twall_noz+273.15, holdups[1].T, vel_noz_vap, d2, rho_vap*holdups[1].mw, mu_vap, cp_vap ,lambda_vap)
        
        println("htc_noz_vap = $(htc_noz_vap)")

        if Δmoles_liq > 0.0
            vel_noz_liq = flowout_liq*holdups[2].v/(pi/4*d2^2)
        else
            vel_noz_liq = 0.0
        end
        htc_noz_liq = alpha_noz(Twall_noz+273.15, holdups[2].T, vel_noz_liq, d2, rho_liq*holdups[2].mw, mu_liq, cp_liq ,lambda_liq)
        
        println("htc_noz_liq = $(htc_noz_liq)")

        htc_ambient_noz = alpha_ambient(Twall_noz+273.15, T_ambient+273.15, u_ambient, d2+2.0*overall_thickness, rho_ambient*mwAir, mu_ambient, cp_ambient , lambda_ambient)

        htc_insulation_noz = lambda_insulation/((d2/2.0)*log((d2/2.0+insulation_thickness)/(d2/2.0)))

        htc_ambient_overall = htc_ambient_noz*htc_insulation_noz/(htc_ambient_noz + htc_insulation_noz)

        println("htc_ambient_overall = $(htc_ambient_overall)")

        Twall_noz = (d2*Δt_vap*htc_noz_vap*(holdups[1].T - 273.15) + d2*Δt_liq*htc_noz_liq*(holdups[2].T - 273.15) + d2*Δt*htc_ambient_overall*T_ambient + cp_steel*rho_steel*((d2 + 2.0*wall_thickness)^2 - d2^2)/4*Twall_noz)/(cp_steel*rho_steel*((d2 + 2.0*wall_thickness)^2 - d2^2)/4 + d2*Δt_vap*htc_noz_vap + d2*Δt_liq*htc_noz_liq + d2*Δt*htc_ambient_overall)
        
        push!(Twallnoz_Plot, Twall_noz)

        println("Twall_noz = $(Twall_noz)")

        # then we use Ferrite to model finite element wall and do a transient step for the Δt

        # Start wall calc using Ferrite
        htcf1 = htc_solid
        Tf1 = holdup_solid.T-273.15
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

        if wtite_VTK
            if iVTK == nVTK
                VTKGridFile("transient-heat-$istep", dh) do vtk
                    write_solution(vtk, dh, u)
                    pvd[time] = vtk
                end
                iVTK = 0
            end
        end

        #At the end of the time loop, we set the previous solution to the current one and go to the next time step.
        uₙ .= u

        points = [Vec(((tank_radius+overall_thickness-0.0001)*cos(theta), (tank_radius+overall_thickness-0.0001)*sin(theta))) for theta in range(-pi/2, pi/2, length = 101)];
        ph = PointEvalHandler(grid, points);
        u_points = evaluate_at_points(ph, dh, uₙ, :u);

        npoints = length(u_points)
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

        println("Twall_outer = $(Twall_outer)")

        q_gp = compute_heat_fluxes(cellvalues, dh, uₙ);
        projector = L2Projector(ip, grid);
        q_projected = project(projector, q_gp, qr);
        points = [Vec((tank_radius*cos(theta), tank_radius*sin(theta))) for theta in range(-pi/2, pi/2, length = 101)];
        ph = PointEvalHandler(grid, points);
        q_points = evaluate_at_points(ph, projector, q_projected);
        u_points = evaluate_at_points(ph, dh, uₙ, :u);

        # this is inner wall temperature
        npoints = length(u_points)
        Twallmin = u_points[1]
        ipoints_min = 1
        for ipoints = 2:npoints
            if u_points[ipoints] < Twallmin
                Twallmin = u_points[ipoints]
                ipoints_min = ipoints
            end
        end

        println("Twallmin = $(Twallmin)")

        push!(Twallmin_Plot, Twallmin)

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

        println("Twall_vap = $(Twall_vap)")
        println("Twall_liq = $(Twall_liq)")
        println("Twall_solid = $(Twall_solid)")

        perimiter = [tank_radius*theta for theta in range(0.0, pi, length = 101)];

    #    np = length(q_points)
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
        Q_wallliq *= 2.0*Δt
        Q_wallvap *= 2.0*Δt

        println("Q_wallvap = $(Q_wallvap)")
        println("Q_wallliq = $(Q_wallliq)")
        println("Q_wallsolid = $(Q_wallsolid)")

        # End wall calc

        # find T at outside of steel

        points = [Vec(((tank_radius+wall_thickness)*cos(theta), (tank_radius+wall_thickness)*sin(theta))) for theta in range(-pi/2, pi/2, length = 101)];
        ph = PointEvalHandler(grid, points);
        u_points = evaluate_at_points(ph, dh, uₙ, :u);

        Twallmin_out = u_points[ipoints_min]
        Twallmin_av = 0.5*(Twallmin_out + Twallmin)

        println("Twallmin outer = $(Twallmin_out)")
        println("Twallmin average = $(Twallmin_av)")

        push!(Twallmin_av_Plot, Twallmin_av)
        push!(Twall_MMDT_Plot, -46.0)

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
            ΔQ[1] = Q_liqvap/holdups[1].moles
            Q_vapliq = Δt*htc_vap_liq*interface_area*(holdups[1].T - holdups[2].T)
            ΔQ[2] = Q_vapliq/holdups[2].moles
        else
            ΔQ[1] = 0.0
            ΔQ[2] = 0.0
        end

        if holdups[1].moles > 0.0
            ΔQ[1] += Q_wallvap/holdups[1].moles
        end

        if holdups[2].moles > 0.0
            ΔQ[2] += Q_wallliq/holdups[2].moles
        end

        if holdup_solid.moles > 0.0 && holdups[2].moles > 0.0
            Q_liqsolid = Δt*htc_liq_solid*interface_area_solid*(holdups[2].T - holdup_solid.T)
            ΔQ_solid = Q_liqsolid/holdup_solid.moles
            Q_solidliq = Δt*htc_liq_solid*interface_area_solid*(holdup_solid.T - holdups[2].T)
            ΔQ[2] += Q_solidliq/holdups[2].moles
        else
            ΔQ_solid = 0.0
        end

        if holdup_solid.moles > 0.0
            ΔQ_solid += Q_wallsolid/holdup_solid.moles
        end

        println("ΔQ[1] = $(ΔQ[1]) J/mol")
        println("ΔQ[2] = $(ΔQ[2]) J/mol")
        println("ΔQ_solid = $(ΔQ_solid) J/mol")

        push!(T1_Plot, holdups[1].T-273.15)
        push!(T2_Plot, holdups[2].T-273.15)

        push!(p1_Plot, holdups[1].p/1e5)
        push!(p2_Plot, holdups[2].p/1e5)
        push!(ps_Plot, holdup_solid.p/1e5)

        push!(m1_Plot, holdups[1].moles*holdups[1].mw/1000)
        push!(m2_Plot, holdups[2].moles*holdups[2].mw/1000)
        push!(ms_Plot, holdup_solid.moles*holdup_solid.mw/1000)

        println("solid mass = $(holdup_solid.moles*holdup_solid.mw) kg")

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
    label=["vapour" "liquid" "minimum inner wall" "minimum average wall" "nozzle" "LTCS MMDT"], 
    [time_Plot,time_Plot,time_Plot,time_Plot,time_Plot,time_Plot],
    [T1_Plot,T2_Plot,Twallmin_Plot,Twallmin_av_Plot,Twallnoz_Plot,Twall_MMDT_Plot],
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

end

delta_P = 0.125e5
pe = 1.0e5
Δtmax = 5.0
nstep = Int((ps - pe)/delta_P)
println("nstep = $(nstep)")
nVTK = 2
BlowDown(model, ps, Ts, zs, tank_volume, tank_radius, tank_length, level, pe, Δtmax, nstep, nVTK)