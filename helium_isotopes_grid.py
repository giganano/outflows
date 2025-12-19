
import numpy as np
import sys
import os

OUTPUT_DIRECTORY = "./outputs-helium-grid"
# OUTPUT_DIRECTORY = "/Volumes/Elements/helium-grid"

Yp = 0.24719 + 2.341e-05
primordial_he_ratios = np.arange(0.1e-4, 2.001e-4, 0.05e-4) # by number -- 39 choices
yieldsolars = np.arange(0.5, 2.6, 0.1) # 21 choices

if os.path.exists(OUTPUT_DIRECTORY): os.system("rm -rf %s" % (OUTPUT_DIRECTORY))
os.system("mkdir %s" % (OUTPUT_DIRECTORY))

for i in range(len(primordial_he_ratios)):
	subdir = ("%.2e" % (primordial_he_ratios[i])).replace('.', 'p')
	os.system("mkdir %s/%s" % (OUTPUT_DIRECTORY, subdir))
	for j in range(len(yieldsolars)):
		subsubdir = ("%.1f" % (yieldsolars[j])).replace('.', 'p')
		os.system("mkdir %s/%s/%s" % (OUTPUT_DIRECTORY, subdir, subsubdir))

for i in range(len(primordial_he_ratios)):
	for j in range(len(yieldsolars)):
		subdir = ("%.2e" % (primordial_he_ratios[i])).replace('.', 'p')
		subsubdir = ("%.1f" % (yieldsolars[j])).replace('.', 'p')
		with open("./src/simulations/nonmetal_inputs.py", 'w') as f:
			primordial_he_ratio_bymass = primordial_he_ratios[i] * 3 / 4
			primordial_he4 = Yp / (1 + primordial_he_ratio_bymass)
			f.write("PRIMORDIAL_HE4 = %.5e\n" % (primordial_he4))
			f.write("PRIMORDIAL_HE3 = %.5e\n" % (Yp - primordial_he4))
			f.write("YIELDSOLAR = %.5e\n" % (yieldsolars[j]))
			f.close()
		cmd = """\
python simulations.py -f --name=%s/%s/%s/output \
--dt=0.05 --nstars=1 --elements=fe_o_he_au \
""" % (OUTPUT_DIRECTORY, subdir, subsubdir)
		if len(sys.argv) > 1: cmd += "--seed=%s" % (sys.argv[1])
		os.system(cmd)

# 		os.system("""\
# python simulations.py -f --name=%s/%s/%s/output \
# --dt=0.05 --nstars=1 --elements=fe_o_he_au --seed=%s
# """ % (OUTPUT_DIRECTORY, subdir, subsubdir, sys.argv[1]))



