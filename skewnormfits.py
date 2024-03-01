r"""
Fit skewnormals to the metallicity distributions predicted by VICE models
in bins of age and radius.
"""

from src.stats import skewnormal_mode_sample
import numpy as np
import vice
import sys

RADIAL_BINS = list(range(16))
AGE_BINS = list(range(11))
ZONE_WIDTH = 0.1

output = vice.output(sys.argv[1])
extra = np.genfromtxt("%s_analogdata.out" % (output.name))
output.stars["absz"] = [abs(_) for _ in extra[:, -1][:output.stars.size[0]]]
stars = output.stars.filter("absz", "<=", 0.5)

with open(sys.argv[2], "w") as out:
	out.write("# age_min [Gyr]    age_max [Gyr]    Rgal_min [kpc]    ")
	out.write("Rgal_max [kpc]    mode([O/H])    mode([Fe/H])\n")
	for i in range(len(AGE_BINS) - 1):
		sub = stars.filter(
			"age", ">=", AGE_BINS[i]).filter(
			"age", "<=", AGE_BINS[i + 1]).filter(
			"mass", ">=", 1)
		for j in range(len(RADIAL_BINS) - 1):
			subsub = sub.filter(
				"zone_final", ">=", int(RADIAL_BINS[j] / ZONE_WIDTH)).filter(
				"zone_final", "<=", int(RADIAL_BINS[j + 1] / ZONE_WIDTH) - 1)
			out.write("%.3e\t%.3e\t" % (AGE_BINS[i], AGE_BINS[i + 1]))
			out.write("%.3e\t%.3e\t" % (RADIAL_BINS[j], RADIAL_BINS[j + 1]))
			try:
				oh = skewnormal_mode_sample(subsub["[o/h]"],
					weights = subsub["mass"])
			except RuntimeError:
				oh = float("nan")
			try:
				feh = skewnormal_mode_sample(subsub["[fe/h]"],
					weights = subsub["mass"])
			except RuntimeError:
				feh = float("nan")
			out.write("%.3e\t%.3e\n" % (oh, feh))
			sys.stdout.write("\rAge = %d - %d Gyr ; R = %d - %d kpc    " % (
				AGE_BINS[i], AGE_BINS[i + 1], RADIAL_BINS[j], RADIAL_BINS[j + 1]))
	sys.stdout.write("\n")
	out.close()


