
from scipy.optimize import bisect
from ..._globals import END_TIME, ZONE_WIDTH
from .utils import exponential, modified_exponential
from .insideout import insideout
from .normalize import normalize
from .gradient import gradient
import numpy as np
import math as m
import vice

TAUSFHMAX = 200
TAURISEMAX = 2 * END_TIME

class subequilibrium(modified_exponential):

	def __init__(self, radius, dt = 0.01, dr = 0.1):
		raw = np.genfromtxt("subequilibrium-input.out")
		radii = [_[0] for _ in raw]
		diff = [abs(_ - radius) for _ in radii]
		idx = diff.index(min(diff))
		super().__init__(radius, rise = raw[idx][1], timescale = raw[idx][2])
		self.norm *= normalize(self, gradient, radius, dt = dt, dr = dr)

def calibrate_subequilibrium(**kwargs):
	radii = [ZONE_WIDTH * (i + 0.5) for i in range(200)]
	with open("subequilibrium-input.out", "w") as out:
		out.write("# radius [kpc]    tau_rise [Gyr]    tau_fall [Gyr]\n")
		for r in radii:
			tau_sfh, tau_rise = find_tausfh_taurise(r, eta = 1)
			if np.isnan(tau_sfh): tau_sfh = TAUSFHMAX
			if np.isnan(tau_rise): tau_rise = TAURISEMAX
			out.write("%.3e\t%.3e\t%.3e\n" % (r, tau_rise, tau_sfh))
		out.close()


def harmonic(*args):
	if len(args) > 0:
		s = 1 / args[0]
		for arg in args[1:]: s -= 1 / arg
		return 1 / s
	else:
		raise TypeError("At least one argument expected.")


class risefall_zalpha:

	def __init__(self, radius, taurise = 2, tausfh = 6, yieldsolar = 1,
		gradient = -0.062, eta = 0, recycling = 0.4, taustar0 = 2, N = 1.5,
		Rg = 3.75):
		self.radius = radius
		self.taurise = taurise
		self.tausfh = tausfh
		self.yalpha = yieldsolar * vice.solar_z['o']
		self.gradient = gradient
		self.eta = eta
		# self.taustar = taustar
		self.taustar = taustar0 * np.exp((N - 1) * radius / Rg)
		self.recycling = recycling

	def __call__(self, time):
		if time == 0: return 0
		taudep = self.taustar / (1 + self.eta - self.recycling)
		term1 = 1 - np.exp(-time / harmonic(taudep, self.tausfh))
		term2 = harmonic(taudep, self.taurise, self.tausfh) / harmonic(
			taudep, self.tausfh)
		term2 *= np.exp(-time / self.taurise) - np.exp(
			-time / harmonic(taudep, self.tausfh))
		return self.equilibrium / (1 - np.exp(-time / self.taurise)) * (
			term1 - term2)

	@property
	def equilibrium(self):
		return self.yalpha / (1 + self.eta - self.recycling -
			self.taustar / self.tausfh)

	@staticmethod
	def at_present_day(radius, **kwargs):
		return risefall_zalpha(radius, **kwargs)(END_TIME)

	@staticmethod
	def gradient(radius, slope = -0.062):
		ralpha = -(np.log(10) * slope)**(-1)
		return vice.solar_z["o"] * np.exp(-(radius - 8) / ralpha)

_DEFAULT_KWARGS_ = {
	"yieldsolar": 1,
	"gradient": -0.062,
	"eta": 0,
	"recycling": 0.4,
	# "taustar": 1.5
	"taustar0": 2,
	"N": 1.5,
	"Rg": 3.75
}
_KWARGS_ = _DEFAULT_KWARGS_.copy()

def find_tausfh_driver(tausfh, radius, taurise):
	diff = risefall_zalpha.at_present_day(radius, tausfh = tausfh,
		taurise = taurise, **_KWARGS_)
	diff -= risefall_zalpha.gradient(radius, slope = _KWARGS_["gradient"])
	return diff

def find_taurise_driver(taurise, radius, tausfh):
	diff = risefall_zalpha.at_present_day(radius, tausfh = tausfh,
		taurise = taurise, **_KWARGS_)
	diff -= risefall_zalpha.gradient(radius, slope = _KWARGS_["gradient"])
	return diff

def find_tausfh_taurise(radius, **kwargs):
	for key in _DEFAULT_KWARGS_.keys():
		if key in kwargs.keys():
			_KWARGS_[key] = kwargs[key]
		else:
			_KWARGS_[key] = _DEFAULT_KWARGS_[key]
	taurise = 2
	if (find_tausfh_driver(0.1, radius, taurise) *
		find_tausfh_driver(TAUSFHMAX, radius, taurise) < 0):
		tausfh = bisect(find_tausfh_driver, 0.1, TAUSFHMAX,
			args = (radius, taurise))
		return [tausfh, taurise]
	elif (find_taurise_driver(0.1, radius, TAUSFHMAX) *
		find_taurise_driver(TAURISEMAX, radius, TAUSFHMAX) < 0):
		taurise = bisect(find_taurise_driver, 0.1, TAURISEMAX,
			args = (radius, TAUSFHMAX))
		return [TAUSFHMAX, taurise]
	else:
		return [float("nan"), float("nan")]

