
from .utils import exponential
from .normalize import normalize
from .gradient import gradient
import math as m
import vice

class simple_exponential(exponential):

	def __init__(self, radius, dt = 0.01, dr = 0.1):
		super().__init__(timescale = 1 + m.exp(radius / 5)) # ~6 Gyr at 8 kpc
		self.norm *= normalize(self, gradient, radius, dt = dt, dr = dr)
