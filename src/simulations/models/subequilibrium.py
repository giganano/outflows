
from .utils import exponential, modified_exponential
from .normalize import normalize
from .gradient import gradient
import numpy as np
import math as m
import vice

class subequilibrium(modified_exponential):

	def __init__(self, radius, dt = 0.01, dr = 0.1):
		raw = np.genfromtxt("subequilibrium-input.out")
		radii = [_[0] for _ in raw]
		diff = [abs(_ - radius) for _ in radii]
		idx = diff.index(min(diff))
		super().__init__(radius, timescale = raw[idx][1], rise = raw[idx][2])
		self.norm = normalize(self, gradient, radius, dt = dt, dr = dr)
		# tausfh = [_[1] for _ in raw]
		# taurise = [_[2] for _ in raw]

		# taustar0 = 2
		# N = 1.5
		# Rg = 4
		# grad = -0.06
		# ralpha = -(grad * m.log(10))**(-1)
		# tausfh = taustar0 * vice.solar_z['o'] / vice.yields.ccsne.settings['o']
		# tausfh *= m.exp((N - 1) * radius / Rg)
		# tausfh *= m.exp((radius - 8) / ralpha)
		# super().__init__(radius, timescale = tausfh)
		# self.norm *= normalize(self, gradient, radius, dt = dt, dr = dr)


