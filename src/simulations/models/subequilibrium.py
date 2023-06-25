
from .utils import exponential, modified_exponential
from .normalize import normalize
from .gradient import gradient
import math as m
import vice

class subequilibrium(exponential):

	def __init__(self, radius, dt = 0.01, dr = 0.1):
		taustar0 = 2
		N = 1.5
		Rg = 4
		grad = -0.06
		ralpha = -(grad * m.log(10))**(-1)
		tausfh = taustar0 * vice.solar_z['o'] / vice.yields.ccsne.settings['o']
		tausfh *= m.exp((N - 1) * radius / Rg)
		tausfh *= m.exp((radius - 8) / ralpha)
		super().__init__(radius, timescale = tausfh)
		self.norm *= normalize(self, gradient, radius, dt = dt, dr = dr)


