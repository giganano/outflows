
from .utils import exponential
from .normalize import normalize
from .gradient import gradient
import math as m
import vice

GRAD = -0.06
TAUSTAR0 = 2
RG = 3
N = 1.5

def tausfh(r):
	ralpha = -(GRAD * m.log(10))**(-1)
	tausfh0 = TAUSTAR0 * m.exp(8 / ralpha) * vice.solar_z['o']
	tausfh0 /= vice.yields.ccsne.settings['o']
	Rsfh = ((N - 1) / RG - 1 / ralpha)**(-1)
	return tausfh0 * m.exp(r / Rsfh)


class simple_exponential(exponential):

	def __init__(self, radius, dt = 0.01, dr = 0.1):
		tsfh = tausfh(radius)
		print(radius, tsfh)
		super().__init__(radius, timescale = tsfh)
		self.norm *= normalize(self, gradient, radius, dt = dt, dr = dr)
