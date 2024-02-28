
import vice

YIELDSOLAR = 1
FE_CC_FRAC = 0.35

vice.yields.ccsne.settings["o"] = YIELDSOLAR * vice.solar_z["o"]
vice.yields.sneia.settings["o"] = 0
vice.yields.ccsne.settings["fe"] = FE_CC_FRAC * YIELDSOLAR * vice.solar_z["fe"]
vice.yields.sneia.settings["fe"] = (1 - FE_CC_FRAC) * YIELDSOLAR * vice.solar_z["fe"]




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
