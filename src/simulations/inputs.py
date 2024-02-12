
import vice
import math as m

METDEPYIELDS = False
R_ETA = 3
YIELDFACTOR = 1
# XH_CGM = -float("inf")
XH_CGM = -0.7

def molecular_tau_star(time, value_today = 2):
	return value_today * ((0.5 + time) / 13.7)**0.5

def sfe_function(time, sigma_sfr):
	mol = molecular_tau_star(time)
	N = plaw_index(time, sigma_sfr)
	return mol * (sigma_sfr * mol / 1e8)**(1 / N - 1)

def plaw_index(time, sigma_sfr):
	mol = molecular_tau_star(time)
	if sigma_sfr > 1e8 / mol:
		return 1
	else:
		return 1.5

def eta_function(radius, time, scale = R_ETA, rsun = 8):
	# etasun = YIELDFACTOR - 0.6
	etasun = 0.6
	return etasun * m.exp((radius - rsun) / scale)
	# return 0
	# return 1
	# return 0.3

def beta_phi_in(radius, time):
	return 0.7

def beta_phi_out(radius, time):
	return 0
	# return 1


##### things that shouldn't need modified below this line
baseline_cc = {
	"o": 0.00713,
	"fe": 0.000473
}
baseline_ia = {
	"o": 0,
	"fe": 0.00077
}

if METDEPYIELDS:
	for elem in ["o", "fe"]:
		def metdepcc(z, maxincrease = 3, zsun = 0.014, element = elem):
			base = PREFACTOR * baseline_cc[element]
			return base * min((z / zsun)**(-1 / 2), maxincrease)
		def metdepia(z, maxincrease = 3, zsun = 0.014, element = elem):
			base = PREFACTOR * baseline_ia[element]
			return base * min((z / zsun)**(-1/ 2), maxincrease)
		vice.yields.ccsne.settings[elem] = metdepcc
		vice.yields.sneia.settings[elem] = metdepia
else:
	vice.yields.ccsne.settings['o'] = YIELDFACTOR * baseline_cc["o"]
	vice.yields.sneia.settings['o'] = YIELDFACTOR * baseline_ia["o"]
	vice.yields.ccsne.settings['fe'] = YIELDFACTOR * baseline_cc["fe"]
	vice.yields.sneia.settings['fe'] = YIELDFACTOR * baseline_ia["fe"]


