r"""
Determine the median stellar age in bins of radius predicted by a model.
"""

import numpy as np
import vice
import sys

RADIAL_BINS = list(range(16))
ZONE_WIDTH = 0.1

output = vice.output(sys.argv[1])
extra = np.genfromtxt("%s_analogdata.out" % (output.name))
output.stars["absz"] = [abs(_) for _ in extra[:, -1][:output.stars.size[0]]]
stars = output.stars.filter("absz", "<=", 0.5)

def median_age(stars):
	indices = np.argsort(stars["age"])
	norm = sum(stars["mass"])
	s = 0
	for idx in indices:
		s += stars["mass"][idx] / norm
		if s > 0.5: return stars["age"][idx]
	return float("nan")

with open(sys.argv[2], "w") as out:
	out.write("# Rgal_min [kpc]    Rgal_max [kpc]   median age [Gyr]\n")
	for i in range(len(RADIAL_BINS) - 1):
		sys.stdout.write("\rR = %d - %d kpc    " % (RADIAL_BINS[i],
			RADIAL_BINS[i + 1]))
		inner = int(RADIAL_BINS[i] / ZONE_WIDTH)
		outer = int(RADIAL_BINS[i + 1] / ZONE_WIDTH)
		stars = output.stars.filter(
			"zone_final", ">=", inner).filter(
			"zone_final", "<=", outer).filter(
			"absz", "<=", 0.5).filter(
			"mass", ">=", 1)
		out.write("%.3e\t%.3e\t%.3e\n" % (RADIAL_BINS[i], RADIAL_BINS[i + 1],
			median_age(stars)))
	sys.stdout.write("\n")
	out.close()

