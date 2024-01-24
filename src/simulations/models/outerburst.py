r"""
This file declares the time-dependence of the star formation history at a
given radius in the outerburst model from Johnson et al. (2021).
"""

from .utils import modified_exponential, gaussian, skewnormal
from .lateburst import _BURST_TIME_, lateburst
from .insideout import _TAU_RISE_, insideout
from .normalize import normalize
from .gradient import gradient
import math as m
import os

_RADIUS_ = 6 # radius in kpc beyond which there is a late starburst

def burst_amplitude(radius):
	if radius > 4:
		testval = m.exp((radius - 4) / 6) - 1
		# if testval > 2: testval = 2
		if testval > 3: testval = 3
		return testval
	else:
		return 0

# class outerburst(modified_exponential, gaussian):

# 	def __init__(self, radius, dt = 0.01, dr = 0.1):
# 		modified_exponential.__init__(self,
# 			timescale = insideout.timescale(radius),
# 			rise = 2)
# 		gaussian.__init__(self, mean = 8, amplitude = burst_amplitude(radius),
# 			std = 0.75)
# 		self._prefactor = 1
# 		self._prefactor = normalize(self, gradient, radius, dt = dt, dr = dr)

# 	def __call__(self, time):
# 		return self._prefactor * modified_exponential.__call__(self, time) * (
# 			1 + gaussian.__call__(self, time))

class outerburst(modified_exponential, skewnormal):

	def __init__(self, radius, dt = 0.01, dr = 0.1):
		modified_exponential.__init__(self,
			timescale = insideout.timescale(radius),
			rise = 2)
		skewnormal.__init__(self, mean = 7, amplitude = burst_amplitude(radius),
			# std = 1, skewness = 2)
			std = 1.5, skewness = 3)
		self.floor = 3.e-3 * m.exp(-radius / 2)
		self._prefactor = 1
		self._prefactor = normalize(self, gradient, radius, dt = dt, dr = dr)

	def __call__(self, time):
		return self.floor + self._prefactor * modified_exponential.__call__(
			self, time) * (1 + skewnormal.__call__(self, time))



# class outerburst(lateburst, insideout):

# 	r"""
# 	The outer-burst SFH model from Johnson et al. (2021).

# 	Parameters
# 	----------
# 	radius : float
# 		The galactocentric radius in kpc of a given annulus in the model.
# 	dt : float [default : 0.01]
# 		The timestep size of the model in Gyr.
# 	dr : float [default : 0.1]
# 		The width of the annulus in kpc.

# 	All attributes and functionality are inherited from ``lateburst`` and
# 	``insideout`` declared in ``src/simulations/models/lateburst.py`` and
# 	``src/simulations/models/insideout.py``.
# 	"""

# 	def __init__(self, radius, dt = 0.01, dr = 0.1):
# 		self._burst = radius > _RADIUS_ # whether or not there will be a burst
# 		if self._burst:
# 			lateburst.__init__(self, radius, dt = dt, dr = dr)
# 		else:
# 			insideout.__init__(self, radius, dt = dt, dr = dr)

# 	def __call__(self, time):
# 		if self._burst:
# 			return lateburst.__call__(self, time)
# 		else:
# 			return insideout.__call__(self, time)

