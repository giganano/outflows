
from plots.mpltoolkit import (mpl_loc, named_colors, markers, fancy_legend,
	load_mpl_presets)
import matplotlib.pyplot as plt
import numpy as np
import vice
import sys
load_mpl_presets()

_ZONE_WIDTH_ = 0.1
_RADII_ = [3, 5, 7, 9, 11, 13]
_LOOKBACK_TIMES_ = [10, 8, 6, 4, 2, 0]
_COLORS_ = ["crimson", "darkorange", "gold", "green", "blue", "darkviolet"]

_SECONDS_PER_GYR_ = 3.1536e16
_KPC_PER_KM_ = 3.24e-17


def main():
	axes = setup_axes()
	plot_evolution(axes)
	axes[0].set_ylim([1.e-5, axes[0].get_ylim()[1]])
	axes[1].set_ylim([axes[1].get_ylim()[0], 1])
	axes[2].set_ylim([1.e5, axes[2].get_ylim()[1]])
	axes[5].set_ylim([-1, 1])
	axes[6].set_ylim([0, axes[6].get_ylim()[1]])
	axes[3].set_ylim([-1, axes[3].get_ylim()[1]])
	axes[4].set_xlim([-1.7, 0.7])
	axes[4].set_ylim([-0.2, 0.5])
	axes[6].set_xlim([-1.5, 1])
	# axes[7].set_ylim([-1, 0])
	axes[7].set_ylim([0, 0.2])
	plt.tight_layout()
	plt.savefig(sys.argv[2])
	plt.clf()


def reference_gradient(radius, slope = -0.06):
	return slope * (radius - 8)


def plot_evolution(axes):
	output = vice.output(sys.argv[1])
	for i, radius in enumerate(_RADII_):
		zone_idx = int(radius / _ZONE_WIDTH_)
		area = np.pi * ((radius + _ZONE_WIDTH_)**2 - radius**2)
		zone = output.zones["zone%d" % (zone_idx)]
		time = zone.history["time"]
		centers = [(a + b) / 2 for a, b in zip(
			zone.mdf["bin_edge_left"], zone.mdf["bin_edge_right"])]
		sigma_sfr = [_ / area for _ in  zone.history["sfr"]]
		sigma_in = [_ / area for _ in zone.history["ifr"]]
		sigma_gas = [_ / area for _ in zone.history["mgas"]]
		oh = zone.history["[o/h]"]
		feh = zone.history["[fe/h]"]
		ofe = zone.history["[o/fe]"]
		kwargs = {
			"c": named_colors()[_COLORS_[i]],
		"label": "%d kpc" % (radius)
		}
		axes[0].plot(time, sigma_sfr, **kwargs)
		axes[1].plot(time, sigma_in, **kwargs)
		axes[2].plot(time, sigma_gas, **kwargs)
		axes[3].plot(time, oh, **kwargs)
		axes[4].plot(feh, ofe, **kwargs)
		axes[6].plot(centers, zone.mdf["dn/d[o/h]"], **kwargs)
		axes[7].plot(time[:-1], velocity_comparison(radius, zone,
			output.zones["zone%d" % (zone_idx + 1)]), **kwargs)

	for i, lookback in enumerate(_LOOKBACK_TIMES_):
		oh = []
		for j in range(int(15.5 / _ZONE_WIDTH_)):
			zone = output.zones["zone%d" % (j)]
			diff = [abs(_ - lookback) for _ in zone.history["lookback"]]
			idx = diff.index(min(diff))
			oh.append(zone.history["[o/h]"][idx])
		kwargs = {
			"c": named_colors()[_COLORS_[i]],
			"label": "%d Gyr ago" % (lookback) if lookback else "Today"
		}
		radii = [_ZONE_WIDTH_ / 2]
		while True:
			if radii[-1] + _ZONE_WIDTH_ > 15.5: break
			radii.append(radii[-1] + _ZONE_WIDTH_)
		axes[5].plot(radii, oh, **kwargs)

	kwargs = {
		"loc": mpl_loc("upper right"),
		"handlelength": 0,
		"ncol": 1
	}
	leg = axes[1].legend(**kwargs)
	fancy_legend(leg, _COLORS_)

	leg = axes[5].legend(**kwargs)
	fancy_legend(leg, _COLORS_)

	axes[5].plot([0, 15], [reference_gradient(_) for _ in [0, 15]],
		c = named_colors()["black"], linestyle = ":")


def velocity_comparison(radius, zone_output, outer_neighbor):
	values = []
	for i in range(len(zone_output.history["time"]) - 1):
		if zone_output.history["sfr"][i]:
			taustar = 1.e-9 * (zone_output.history["mgas"][i] / 
				zone_output.history["sfr"][i])
		else:
			taustar = float("inf")
		eta = zone_output.history["eta_0"][i]
		r = zone_output.history["r_eff"][i]
		mdotstar = zone_output.history["sfr"][i]
		mddotstar = (zone_output.history["sfr"][i + 1] - mdotstar) / (
			zone_output.history["time"][1] - zone_output.history["time"][0])
		if mdotstar:
			deriv = mddotstar / mdotstar
		else:
			deriv = float("inf")
		if zone_output.history["z(o)"][i]:
			yieldratio = zone_output.ccsne_yields["o"] / (taustar *
				zone_output.history["z(o)"][i])
		else:
			yieldratio = float("nan")
		a = (1 + eta - r) / taustar + deriv - yieldratio

		area = np.pi * ((radius + _ZONE_WIDTH_)**2 - radius**2)
		area_outer = np.pi * ((radius + 2 * _ZONE_WIDTH_)**2 -
			(radius + _ZONE_WIDTH_)**2)
		sigma_g = zone_output.history["mgas"][i] / area
		sigma_outer = outer_neighbor.history["mgas"][i] / area_outer
		dlnsigmag_dr = (sigma_outer - sigma_g) / (sigma_g * _ZONE_WIDTH_)
		if zone_output.history["z(o)"][i]:
			dlnz_dr = (outer_neighbor.history["z(o)"][i] -
				zone_output.history["z(o)"][i]) / (
				zone_output.history["z(o)"][i] * _ZONE_WIDTH_)
		else:
			dlnz_dr = float("nan")
		b = (1 / radius + dlnsigmag_dr + dlnz_dr)**(-1)

		vg = float(sys.argv[3]) * _KPC_PER_KM_ * _SECONDS_PER_GYR_

		if abs(vg) > abs(a * b):
			values.append(abs(vg) - abs(a * b))
		else:
			values.append(float("nan"))

	return values


def setup_axes():
	plt.clf()
	fig = plt.figure(figsize = (15, 15))
	axes = []
	for i in range(8):
		axes.append(fig.add_subplot(331 + i))
		if i < 4 or i == 7:
			axes[i].set_xlabel(r"Time [Gyr]")
			axes[i].set_xlim([0, 13.2])
		else: pass
		if i < 3: axes[i].set_yscale("log")
	axes[0].set_ylabel(r"$\dot{\Sigma}_\star$ [M$_\odot$ kpc$^{-2}$ yr$^{-1}$]")
	axes[1].set_ylabel(
		r"$\dot{\Sigma}_\text{in}$ [M$_\odot$ kpc$^{-2}$ yr$^{-1}$]")
	axes[2].set_ylabel(r"$\Sigma_g$ [M$_\odot$ kpc$^{-2}$]")
	axes[3].set_ylabel(r"[O/H]$_\text{gas}$")
	axes[4].set_xlabel(r"[Fe/H]")
	axes[4].set_ylabel(r"[O/Fe]")
	axes[5].set_xlabel(r"$R_\text{gal}$ [kpc]")
	axes[5].set_ylabel(r"[O/H]$_\text{gas}$")
	axes[6].set_xlabel("[O/H]")
	axes[6].set_ylabel(r"$dN_\star/d$[O/H]")
	axes[7].set_xlabel(r"Time [Gyr]")
	axes[7].set_ylabel(r"Flows are overpowering sink if $>$ 0")
	return axes


if __name__ == "__main__": main()

