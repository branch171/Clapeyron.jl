
#=
CEmerald::SuperTRAPP::SuperTRAPP(Basis para2) : basis(para2)
{
	E = {
		-14.113294896,
		968.22940153,
		13.686545032,
		-12511.628378,
		0.0168910864,
		43.527109444,
		7659.4543472
	};

	C = {
		15.2583985944,
		5.29917319127,
		-3.05330414748,
		0.450477583739,
		1.03144050679,
		-0.185480417707
	};

	paramsref.cname = L"Propane";
	paramsref.type = L"base";
	paramsref.Mw = 44.0956;
	paramsref.CN = 0.0;
	paramsref.SG = 0.0;
	paramsref.Tb = 0.0;
	paramsref.Tc = 369.85;
	paramsref.Pc = 42.4766;
	paramsref.Zc = 0.280363;
	paramsref.Omega = 0.152;
	paramsref.m = 2.002;
	paramsref.sigma = 3.6184;
	paramsref.eta = 208.11;
	paramsref.etaAB = 0.0;
	paramsref.kappaAB = 0.0;
	paramsref.S = { 0.0, 0.0 };
	paramsref.CpIdeal = { 4.0,	31.0,	290.0,	0.2675,	3.124, -4.826, 0.0 };

	//CPA	
	paramsref.cpa_parameters.C = 0.0;
	paramsref.cpa_parameters.B = 0.0;
	paramsref.cpa_parameters.m = 0.0;
	paramsref.cpa_parameters.association_parameters.type = L"none";
	paramsref.cpa_parameters.association_parameters.kappaAB = 0.0;
	paramsref.cpa_parameters.association_parameters.etaAB = 0.0;
	paramsref.cpa_parameters.association_parameters.association_structure = { 0,  0, 0 };
	paramsref.cpa_parameters.association_parameters.doners_per_site = { 0,  0, 0 };

	//PCSAFT	
	paramsref.pcsaft_parameters.m = 2.002;
	paramsref.pcsaft_parameters.sigma = 3.6184;
	paramsref.pcsaft_parameters.eta = 208.11;
	paramsref.pcsaft_parameters.association_parameters.type = L"none";
	paramsref.pcsaft_parameters.association_parameters.kappaAB = 0.0;
	paramsref.pcsaft_parameters.association_parameters.etaAB = 0.0;
	paramsref.pcsaft_parameters.association_parameters.association_structure = { 0,  0, 0 };
	paramsref.pcsaft_parameters.association_parameters.doners_per_site = { 0,  0, 0 };

	//SAFTVRMie	
	paramsref.saftvrmie_parameters.m = 2.002;
	paramsref.saftvrmie_parameters.sigma = 3.6184;
	paramsref.saftvrmie_parameters.eta = 208.11;
	paramsref.saftvrmie_parameters.La = 6.0;
	paramsref.saftvrmie_parameters.Lr = 12.0;
	paramsref.saftvrmie_parameters.association_parameters.type = L"none";
	paramsref.saftvrmie_parameters.association_parameters.kappaAB = 0.0;
	paramsref.saftvrmie_parameters.association_parameters.etaAB = 0.0;
	paramsref.saftvrmie_parameters.association_parameters.association_structure = { 0,  0, 0 };
	paramsref.saftvrmie_parameters.association_parameters.doners_per_site = { 0,  0, 0 };
	
}
double CEmerald::SuperTRAPP::mu(const double T, const double rho, const vector<double> z) const
{
	size_t nc;
	double third = 1.0 / 3.0, temp, temp2, Mwij, Mwik1, Mwik2, sij, deltaij, deltajk;
	double Mw_mix, neta = 0.0, MwRef, ZcRef, TcRef, PcRef, VcRef, rhocRef, OmegaRef;

	vector<Component> params;

	params = basis.component;

	nc = z.size();

	Mw_mix = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		Mw_mix += z[i] * params[i].Mw;
	}

	vector<double > Vc(nc);

	for (size_t i = 0; i < nc; i++)
	{
		Vc[i] = params[i].Zc * RGas * params[i].Tc / params[i].Pc / 10.0;
	}

	MwRef = paramsref.Mw;
	ZcRef = paramsref.Zc;
	TcRef = paramsref.Tc;
	PcRef = paramsref.Pc;
	VcRef = ZcRef * RGas * TcRef / PcRef / 10.0;
	rhocRef = 1000.0 / VcRef;
	OmegaRef = paramsref.Omega;

	vector<double > f(nc);
	for (size_t i = 0; i < nc; i++)
	{
		f[i] = params[i].Tc / TcRef * (1.0 + (params[i].Omega - OmegaRef) * (0.05202976 - 0.7498189 * log(T / params[i].Tc)));
	}

	vector<double > h(nc);
	for (size_t i = 0; i < nc; i++)
	{
		h[i] = Vc[i] / VcRef * ZcRef / params[i].Zc * (1.0 + (params[i].Omega - OmegaRef) * (0.1435971 - 0.2821562 * log(T / params[i].Tc)));
	}

	double hm = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			temp = (pow(h[i], third) + pow(h[j], third));
			temp = temp * temp * temp / 8.0;
			hm += z[i] * z[j] * temp;
		}
	}

	double fm = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			temp = (pow(h[i], third) + pow(h[j], third));
			temp = temp * temp * temp / 8.0;
			fm += z[i] * z[j] * temp * sqrt(f[i] * f[j]);
		}
	}
	fm = fm / hm;

	double To, rhoo, zeta;

	To = T / fm;
	rhoo = rho * hm;

	double G1, G2, G3, DnetaR;

	G1 = exp(E[0] + E[1] / To);
	temp = pow(To, 1.5);
	G2 = E[2] + E[3] / temp;
	G3 = E[4] + E[5] / To + E[6] / To / To;

	DnetaR = G1 * (exp(pow(rhoo, 0.1) * G2 + sqrt(rhoo) * (rhoo / rhocRef - 1.0) * G3) - 1.0);

	vector<double > s(nc);
	for (size_t i = 0; i < nc; i++)
	{
		s[i] = 4.771 * pow(h[i], third);
	}

	double ss2 = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		ss2 += z[i] * s[i] * s[i];
	}

	double ss3 = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		ss3 += z[i] * s[i] * s[i] * s[i];
	}

	zeta = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		zeta += z[i] * s[i] * s[i] * s[i];
	}
	zeta = 6.023e-4 * pi / 6.0 * rho * zeta;

	double Fnetam = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			temp = (pow(h[i], third) + pow(h[j], third));
			temp = temp * temp * temp / 8.0;
			temp = pow(temp, 4.0 * third);
			Mwij = 2.0 * params[i].Mw * params[j].Mw / (params[i].Mw + params[j].Mw);
			temp2 = sqrt(f[i] * f[j]);
			Fnetam += z[i] * z[j] * temp * sqrt(temp2 * Mwij);
		}
	}
	Fnetam = Fnetam / sqrt(MwRef) / hm / hm;

	vector<vector<double >> g;
	g.resize(nc);
	for (size_t i = 0; i < nc; i++)
	{
		g[i].resize(nc);
	}

	temp = 1.0 - zeta;
	double theta;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			theta = s[i] * s[j] / (s[i] + s[j]) * ss2 / ss3;
			g[i][j] = 1.0 / temp + 3.0 * zeta / temp / temp * theta + 2.0 * zeta * zeta / temp / temp / temp * theta * theta;
		}
	}

	vector<double > Y(nc);
	for (size_t i = 0; i < nc; i++)
	{
		temp = 0.0;
		for (size_t j = 0; j < nc; j++)
		{
			Mwij = params[j].Mw / (params[i].Mw + params[j].Mw);
			sij = (s[i] + s[j]) / 2.0;
			temp += z[j] * Mwij * sij * sij * sij * g[i][j];
		}
		Y[i] = z[i] * (1.0 + 8.0 * pi / 15.0 * 6.023e-4 * rho * temp);
	}

	vector<vector<double >> neta0;
	neta0.resize(nc);
	for (size_t i = 0; i < nc; i++)
	{
		neta0[i].resize(nc);
	}

	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			Mwij = 2.0 * params[i].Mw * params[j].Mw / (params[i].Mw + params[j].Mw);
			sij = (s[i] + s[j]) / 2.0;
			neta0[i][j] = 26.692 * sqrt(Mwij * T) / sij / sij;
		}
	}

	vector<vector<double >> B;
	B.resize(nc);
	for (size_t i = 0; i < nc; i++)
	{
		B[i].resize(nc);
	}

	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			temp = 0.0;
			deltaij = 0.0;
			if (i == j) deltaij = 1.0;
			for (size_t k = 0; k < nc; k++)
			{
				Mwik1 = params[k].Mw / (params[i].Mw + params[k].Mw);
				Mwik2 = params[i].Mw / params[k].Mw;
				deltajk = 0.0;
				if (j == k) deltajk = 1.0;
				temp += z[i] * z[k] * g[i][k] / neta0[i][k] * Mwik1 * Mwik1 * ((1.0 + 5.0 / 3.0 * Mwik2) * deltaij - 2.0 / 3.0 * Mwik2 * deltajk);
			}
			B[i][j] = 2.0 * temp;
		}
	}

	vector<double > beta(nc);
	LUdcmpv alu(B);
	alu.solve(Y, beta);

	temp = 2.0 * pi / 3.0 * 6.023e-4;
	temp *= temp;
	temp *= 48.0 / 25.0 / pi * rho * rho;

	double neta_Enskog_m = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		neta_Enskog_m += Y[i] * beta[i];
	}

	temp2 = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			sij = (s[i] + s[j]) / 2.0;
			temp2 += z[i] * z[j] * sij * sij * sij * sij * sij * sij * neta0[i][j] * g[i][j];
		}
	}

	neta_Enskog_m += temp * temp2;

	double sx = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			sij = (s[i] + s[j]) / 2.0;
			sx += z[i] * z[j] * sij * sij * sij;
		}
	}
	sx = pow(sx, third);

	double Mwx = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			Mwij = 2.0 * params[i].Mw * params[j].Mw / (params[i].Mw + params[j].Mw);
			sij = (s[i] + s[j]) / 2.0;
			Mwx += z[i] * z[j] * sqrt(Mwij) * sij * sij * sij * sij;
		}
	}
	Mwx = Mwx / sx / sx / sx / sx;
	Mwx = Mwx * Mwx;

	double zetax = 6.023e-4 * pi / 6.0 * rho * sx * sx * sx;

	temp = (1.0 - zetax);
	double gxx = 1.0 / temp + 3.0 * zetax / temp / temp * 0.5 + 2.0 * zetax * zetax / temp / temp / temp * 0.25;

	double Yx = 1.0 + 8.0 * pi / 15.0 * 6.023e-4 * rho * sx * sx * sx * gxx / 2.0;

	double neta0x = 26.692 * sqrt(Mwx * T) / sx / sx;

	double Bxx = gxx / neta0x;

	double betax = Yx / Bxx;

	temp = 2.0 * pi / 3.0 * 6.023e-4;
	temp *= temp;
	temp *= 48.0 / 25.0 / pi * rho * rho;

	double neta_Enskog_x = betax * Yx + temp * sx * sx * sx * sx * sx * sx * neta0x * gxx;

	double Dneta_Enskog = 0.1 * (neta_Enskog_m - neta_Enskog_x);

	for (size_t i = 0; i < nc; i++)
	{
		s[i] = 0.809 * pow(Vc[i], third);
	}

	double sm = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			sij = sqrt(s[i] * s[j]);
			sm += z[i] * z[j] * sij * sij * sij;
		}
	}
	sm = pow(sm, third);

	vector<double > eta(nc);
	for (size_t i = 0; i < nc; i++)
	{
		eta[i] = params[i].Tc / 1.2593;
	}

	double etaij, etam = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			sij = sqrt(s[i] * s[j]);
			etaij = sqrt(eta[i] * eta[j]);
			etam += z[i] * z[j] * etaij * sij * sij * sij;
		}
	}
	etam = etam / sm / sm / sm;

	double Mwm = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			sij = sqrt(s[i] * s[j]);
			etaij = sqrt(eta[i] * eta[j]);
			Mwij = 2.0 * params[i].Mw * params[j].Mw / (params[i].Mw + params[j].Mw);
			Mwm += z[i] * z[j] * etaij * sij * sij * sqrt(Mwij);
		}
	}
	Mwm = Mwm / etam / sm / sm;
	Mwm = Mwm * Mwm;

	double Omegaij, Omegam = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			sij = sqrt(s[i] * s[j]);
			Omegaij = (params[i].Omega + params[i].Omega) / 2.0;
			Omegam += z[i] * z[j] * Omegaij * sij * sij * sij;
		}
	}
	Omegam = Omegam / sm / sm / sm;

	double Fc = 1.0 - 0.2756 * Omegam;

	double FF = 1.16145 * pow(T / etam, -0.14874) + 0.52487 * exp(-0.77320 * T / etam) + 2.16178 * exp(-2.43787 * T / etam);

	double neta0m = 2.6692 * Fc * sqrt(Mwm * T) / sm / sm / FF;

	neta = (neta0m + Fnetam * DnetaR + Dneta_Enskog) * 1.0e-6;

	return neta; // Pa.s
}
double CEmerald::SuperTRAPP::lambda(const double T, const double rho, const vector<double> z) const
{
	size_t nc;
	double third = 1.0 / 3.0, temp, temp2, invMwij;
	double Mw_mix, lambda = 0.0, MwRef, ZcRef, TcRef, PcRef, VcRef, rhocRef, OmegaRef;

	vector<Component> params;

	params = basis.component;

	nc = z.size();

	Mw_mix = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		Mw_mix += z[i] * params[i].Mw;
	}

	vector<double > Vc(nc);

	for (size_t i = 0; i < nc; i++)
	{
		Vc[i] = params[i].Zc * RGas * params[i].Tc / params[i].Pc / 10.0;
	}

	MwRef = paramsref.Mw;
	ZcRef = paramsref.Zc;
	TcRef = paramsref.Tc;
	PcRef = paramsref.Pc;
	VcRef = ZcRef * RGas * TcRef / PcRef / 10.0;
	rhocRef = 1000.0 / VcRef;
	OmegaRef = paramsref.Omega;

	vector<double > f(nc);
	for (size_t i = 0; i < nc; i++)
	{
		f[i] = params[i].Tc / TcRef * (1.0 + (params[i].Omega - OmegaRef) * (0.05202976 - 0.7498189 * log(T / params[i].Tc)));
	}

	vector<double > h(nc);
	for (size_t i = 0; i < nc; i++)
	{
		h[i] = Vc[i] / VcRef * ZcRef / params[i].Zc * (1.0 + (params[i].Omega - OmegaRef) * (0.1435971 - 0.2821562 * log(T / params[i].Tc)));
	}

	double hm = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			temp = (pow(h[i], third) + pow(h[j], third));
			temp = temp * temp * temp / 8.0;
			hm += z[i] * z[j] * temp;
		}
	}

	double fm = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			temp = (pow(h[i], third) + pow(h[j], third));
			temp = temp * temp * temp / 8.0;
			fm += z[i] * z[j] * temp * sqrt(f[i] * f[j]);
		}
	}
	fm = fm / hm;

	double To, rhoo;

	To = T / fm;
	rhoo = rho * hm;

	double Flambdam = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			temp = (pow(h[i], third) + pow(h[j], third));
			temp = temp * temp * temp / 8.0;
			temp = pow(temp, -4.0 * third);
			invMwij = (1.0 / params[i].Mw + 1.0 / params[j].Mw) / 2.0;
			temp2 = sqrt(f[i] * f[j]);
			Flambdam += z[i] * z[j] * temp * sqrt(temp2 * invMwij);
		}
	}
	Flambdam = Flambdam * sqrt(MwRef) * pow(hm, 2.0 * third);

	double Omegam = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		Omegam += z[i] * params[i].Omega;
	}

	double Xm = sqrt((1.0 + 2.186634 * (Omegam - OmegaRef))) / (1.0 - 0.5050059 * (Omegam - OmegaRef));

	double Tr = To / TcRef;
	double rhor = rhoo / rhocRef;

	double DlambdaR = C[0] * rhor + C[1] * rhor * rhor * rhor + (C[2] + C[3] / Tr) * rhor * rhor * rhor * rhor + (C[4] + C[5] / Tr) * rhor * rhor * rhor * rhor * rhor;

	vector<double > s(nc);
	for (size_t i = 0; i < nc; i++)
	{
		s[i] = 0.809 * pow(Vc[i], third);
	}

	vector<double > eta(nc);
	for (size_t i = 0; i < nc; i++)
	{
		eta[i] = params[i].Tc / 1.2593;
	}

	double y, F;
	vector<double > Cp(nc);
	for (size_t i = 0; i < nc; i++)
	{
		y = T / (T + params[i].CpIdeal[2]);
		F = params[i].CpIdeal[3] + params[i].CpIdeal[4] * y + params[i].CpIdeal[5] * y * y + params[i].CpIdeal[6] * y * y * y;
		Cp[i] = (params[i].CpIdeal[0] + y * y * (params[i].CpIdeal[1] - params[i].CpIdeal[0]) * (1.0 + (y - 1.0) * F)) * RGas / 100.0;
	}

	double Fc, FF;
	vector<double > neta0x(nc);
	for (size_t i = 0; i < nc; i++)
	{
		Fc = 1.0 - 0.2756 * params[i].Omega;
		FF = 1.16145 * pow(T / eta[i], -0.14874) + 0.52487 * exp(-0.77320 * T / eta[i]) + 2.16178 * exp(-2.43787 * T / eta[i]);
		neta0x[i] = 2.6692 * Fc * sqrt(params[i].Mw * T) / s[i] / s[i] / FF;
	}

	vector<double > lambda0x(nc);
	for (size_t i = 0; i < nc; i++)
	{
		lambda0x[i] = 15.0 / 4.0 * RGas / 100.0 * neta0x[i] * 1.0e-6 / params[i].Mw * 1000.0;
	}

	vector<double > lambdaintx(nc);
	for (size_t i = 0; i < nc; i++)
	{
		lambdaintx[i] = 1.32 * neta0x[i] * 1.0e-6 / params[i].Mw * (Cp[i] - 2.5 * RGas / 100.0) * 1000.0;
	}

	double lambda0ij, lambda0m = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			lambda0ij = 2.0 * lambda0x[i] * lambda0x[j] / (lambda0x[i] + lambda0x[j]);
			lambda0m += z[i] * z[j] * lambda0ij;
		}
	}

	double lambdaintij, lambdaintm = 0.0;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < nc; j++)
		{
			lambdaintij = 2.0 * lambdaintx[i] * lambdaintx[j] / (lambdaintx[i] + lambdaintx[j]);
			lambdaintm += z[i] * z[j] * lambdaintij;
		}
	}

	lambda = lambda0m + lambdaintm + Flambdam * Xm * DlambdaR / 1000.0;

	return lambda; // W/m/K
}

=#

using Clapeyron

# Propane
MwRef = 44.0956
ZcRef = 0.280363
TcRef = 369.85
PcRef = 42.4766

fluid = ["carbon dioxide","nitrogen"]
nfluid = length(fluid)
model = GERG2008(fluid)

mw = model.params.Mw.values

println("mw = $(mw)")

Tc = model.params.Tc.values

println("Tc = $(Tc)")

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

xCO2 = 0.997275
x = [xCO2,1.0-xCO2]

pbase = 27.0e5
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

Tvap = Tin
rhovap = 1.0/vf[1]
#println("rhovap = $(rhovap)")
zvap = xf[1]
mu_vap = SuperTRAPP_mu(model, Tvap, rhovap, zvap)
lambda_vap = SuperTRAPP_lambda(model, Tvap, rhovap, zvap)

Tliq = Tin
rholiq = 1.0/vf[2]
#println("rholiq = $(rholiq)")
zliq = xf[2]
mu_liq = SuperTRAPP_mu(model, Tliq, rholiq, zliq)
lambda_liq = SuperTRAPP_lambda(model, Tliq, rholiq, zliq)

println("mu_vap, mu_liq = $(mu_vap) Pa.s, $(mu_liq) Pa.s")

println("lambda_vap, lambda_liq = $(lambda_vap) W/m/K, $(lambda_liq) W/m/K")