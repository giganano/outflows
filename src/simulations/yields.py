
import vice
# from vice.yields.presets import JW20
# vice.yields.sneia.settings['fe'] *= 10**0.1
vice.yields.ccsne.settings['o'] = 0.01
vice.yields.sneia.settings['o'] = 0
vice.yields.ccsne.settings['fe'] = 8e-4
vice.yields.sneia.settings['fe'] = 0.0014

# for elem in ["o", "fe"]:
# 	vice.yields.ccsne.settings[elem] *= 2 / 3
# 	vice.yields.sneia.settings[elem] *= 2 / 3


