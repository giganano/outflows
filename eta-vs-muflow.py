#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
from plots.mpltoolkit import (named_colors, mpl_loc, fancy_legend,
	load_mpl_presets)
load_mpl_presets()
mpl.rcParams["axes.linewidth"] = 0.5
mpl.rcParams["figure.titlesize"] = 14
mpl.rcParams["axes.titlesize"] = 14
mpl.rcParams["axes.labelsize"] = 14
mpl.rcParams["xtick.labelsize"] = 14
mpl.rcParams["ytick.labelsize"] = 14
mpl.rcParams["legend.fontsize"] = 14
import numpy as np
import sys
print(sys.version_info)


def eta(r, yieldsolar = 3, gradient = -0.06):
	return yieldsolar * np.exp(-gradient * np.log(10) * (r - 8)) - 0.6

def mu_flow(r, taustar0 = 2, N = 1.5, Rg = 3.75, gradient = -0.06, vgas = -1):
	taustar = taustar0 * np.exp((N - 1) * r / Rg)
	mu = 1 / r - 1 / Rg + (np.log(10) * gradient)
	mu *= -taustar * vgas
	return mu

if __name__ == "__main__":
	fig = plt.figure(figsize = (3.5, 3.5))
	ax = fig.add_subplot(111)
	ax.set_xlabel(r"$R_\text{gal}$ [kpc]")
	# ax.set_ylabel(r"$\eta$ (solid) ; $-\mu_\text{flow}$ (dotted)")
	ax.set_ylabel(r"Coefficient Value")
	ax.set_xlim([0, 15])
	ax.set_xticks([0, 5, 10, 15])
	ax.set_ylim([-1, 10])
	ax.set_yticks([0, 5, 10])

	line1 = ax.plot([-1, -2], [-1, -2], c = named_colors()["black"],
		linestyle = "-", label = r"$\eta$")[0]
	line2 = ax.plot([-1, -2], [-1, -2], c = named_colors()["black"],
		linestyle = ":", label = r"$-\mu_\text{flow}$")[0]
	kwargs = {
		"loc": mpl_loc("lower right"),
		"handlelength": 1.5,
		"fontsize": 12
	}
	leg = ax.legend(**kwargs)
	ax.add_artist(leg)
	line1.remove()
	line2.remove()

	gradients = [-0.04, -0.06, -0.08]
	colors = ["crimson", "black", "dodgerblue"]
	xvals = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 1000)
	for i in range(len(gradients)):
		yvals = [eta(x, gradient = gradients[i]) for x in xvals]
		kwargs = {
			"c": named_colors()[colors[i]],
			"label": r"$\nabla\text{[O/H]} = %g$ kpc$^{-1}$" % (
				gradients[i])
		}
		ax.plot(xvals, yvals, **kwargs)

	vgas = [-1.2, -1.5, -1.8]
	for i in range(len(vgas)):
		yvals = [-mu_flow(x, vgas = vgas[i]) for x in xvals]
		kwargs = {
			"linestyle": ":",
			"c": named_colors()[colors[i]],
			"label": r"$\langle v_\text{g} \rangle = %g$ km/s" % (vgas[i])
		}
		ax.plot(xvals, yvals, **kwargs)

	kwargs = {
		"loc": mpl_loc("upper left"),
		"ncol": 1,
		"handlelength": 0,
		"fontsize": 12
	}
	leg = ax.legend(**kwargs)
	colors.extend(colors)
	fancy_legend(leg, colors)

	plt.tight_layout()
	for ext in ["pdf", "jpeg"]:
		kwargs = {}
		if ext == "jpeg": kwargs["dpi"] = 200
		plt.savefig("./eta-vs-muflow.%s" % (ext))
else: pass
