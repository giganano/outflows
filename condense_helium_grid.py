
import numpy as np
import vice
import sys

output_file = "helium_grid_presentday.npy"
output_directory = "/Volumes/Elements/helium-grid"
zone_width = 0.1
lookback = 0

primordial_he_ratios = np.arange(0.1e-4, 2.001e-4, 0.05e-4) # by number -- 39 choices
yieldsolars = np.arange(0.7, 2.6, 0.1) # 19 choices
rgal = [zone_width * (i + 0.5) for i in range(int(20 / zone_width))]

results = np.zeros((
	len(primordial_he_ratios),
	len(yieldsolars),
	len(rgal),
	))

for i in range(len(primordial_he_ratios)):
	for j in range(len(yieldsolars)):
		subdir = ("%.2e" % (primordial_he_ratios[i])).replace('.', 'p')
		subsubdir = ("%.1f" % (yieldsolars[j])).replace('.', 'p')
		name = "%s/%s/%s/output.vice" % (output_directory, subdir, subsubdir)
		sys.stdout.write("\r%s             " % (name))
		output = vice.output(name)
		diff = [abs(_ - lookback) for _ in output.zones["zone0"].history["lookback"]]
		idx = diff.index(min(diff))
		for k in range(len(rgal)):
			he3 = output.zones["zone%d" % (k)].history["mass(au)"][idx]
			he4 = output.zones["zone%d" % (k)].history["mass(he)"][idx]
			results[i][j][k] = he3 / he4 * 4 / 3

sys.stdout.write("\n")
np.save(output_file, results)

