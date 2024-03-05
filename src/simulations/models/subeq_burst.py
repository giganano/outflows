
from .subequilibrium import subequilibrium
from .lateburst import _BURST_TIME_
from .outerburst import burst_amplitude
# from .utils import gaussian
from .utils import skewnormal
from .normalize import normalize
from .gradient import gradient


# class subeq_burst(subequilibrium, gaussian):
class subeq_burst(subequilibrium, skewnormal):

	def __init__(self, radius, dt = 0.01, dr = 0.1, **kwargs):
		self._prefactor = 1
		# gaussian.__init__(self, mean = _BURST_TIME_, amplitude = 1.5)
		skewnormal.__init__(self, mean = 7, amplitude = burst_amplitude(radius),
			std = 1.5, skewness = 3)
		subequilibrium.__init__(self, radius, dt = dt, dr = dr, **kwargs)
		self._prefactor *= normalize(self, gradient ,radius, dt = dt, dr = dr)

	def __call__(self, time):
		return self._prefactor * subequilibrium.__call__(self, time) * (
			# 1 + gaussian.__call__(self, time)
			1 + skewnormal.__call__(self, time)
		)

