
from .utils import exponential, modified_exponential
from .normalize import normalize
from .gradient import gradient
import math as m
import vice

class subequilibrium(exponential):

	def __init__(self, radius, dt = 0.01, dr = 0.1):
		taustar0 = 2
		yieldsolar = vice.yields.ccsne.settings['o'] / vice.solar_z['o']
		N = 1.5
		Rg = 4
		grad = -0.06
		super().__init__(radius,
			timescale = taustar0 / yieldsolar * m.exp(
				(N - 1) * radius / Rg + (radius - 8) * grad * m.log(10))
			)
		self.norm *= normalize(self, gradient, radius, dt = dt, dr = dr)


