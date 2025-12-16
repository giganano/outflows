
import vice

YIELDSOLAR = 2
FE_CC_FRAC = 0.35
XH_CGM = -float("inf")
ETA_VARY = False
METDEPYIELDS = False

vice.yields.ccsne.settings["o"] = YIELDSOLAR * vice.solar_z["o"]
vice.yields.sneia.settings["o"] = 0
vice.yields.ccsne.settings["fe"] = FE_CC_FRAC * YIELDSOLAR * vice.solar_z["fe"]
vice.yields.sneia.settings["fe"] = (1 - FE_CC_FRAC) * YIELDSOLAR * vice.solar_z["fe"]

class metdepyield:

	def __init__(self, baseline, maxincrease = 3, plawindex = -0.5,
		zsun = 0.014):
		self.baseline = baseline
		self.maxincrease = maxincrease
		self.plawindex = plawindex
		self.zsun = zsun

	def __call__(self, z):
		if z:
			prefactor = min((z / self.zsun)**self.plawindex, self.maxincrease)
		else:
			prefactor = self.maxincrease
		return self.baseline * prefactor

if METDEPYIELDS:
	for elem in ["o", "fe"]:
		vice.yields.ccsne.settings[elem] = metdepyield(
			vice.yields.ccsne.settings[elem])
		vice.yields.sneia.settings[elem] = metdepyield(
			vice.yields.sneia.settings[elem])
else: pass


# from vice.yields.presets import JW20
# vice.yields.sneia.settings['fe'] *= 10**0.1

# for elem in ["o", "fe"]:
# 	vice.yields.ccsne.settings[elem] /= 3
# 	vice.yields.sneia.settings[elem] /= 3

# vice.yields.ccsne.settings['o'] = 0.015
# vice.yields.sneia.settings['o'] = 0
# vice.yields.ccsne.settings['fe'] = 0.0012
# vice.yields.sneia.settings['fe'] = 0.0021

# vice.yields.ccsne.settings['o'] = 0.01
# vice.yields.sneia.settings['o'] = 0
# vice.yields.ccsne.settings['fe'] = 8e-4
# vice.yields.sneia.settings['fe'] = 0.0014

# vice.yields.ccsne.settings['o'] = 0.005
# vice.yields.sneia.settings['o'] = 0
# vice.yields.ccsne.settings['fe'] = 4e-4
# vice.yields.sneia.settings['fe'] = 0.0007
