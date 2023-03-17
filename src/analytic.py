r"""
This file implements analytic solutions to :math:`Z_\alpha(t)` and its first
two time-derivaties. See ../latex/onezone-mdfs.pdf.
"""

__all__ = ["Zalpha", "Zdotalpha", "Zdotdotalpha"]
import math as m


class Zalpha:

	def __init__(self, y_alpha = 0.015, taurise = 2, tausfh = 6, taustar = 5,
		eta = 2.5, recycling = 0.4):
		self.y_alpha = y_alpha
		self.taurise = taurise
		self.tausfh = tausfh
		self.taustar = taustar
		self.eta = eta
		self.recycling = recycling


	def __call__(self, time):
		if time == 0: return 0
		term1 = 1 / (1 - m.exp(-time / self.taurise))
		term2 = self.y_alpha / (1 + self.eta - self.recycling)
		return term1 * term2 * self._g(time)


	@staticmethod
	def harmonic_timescale(*args):
		s = 1 / args[0]
		for arg in args[1:]: s -= 1 / arg
		return 1 / s


	def _g(self, time):
		term1 = self.harmonic_timescale(self.taudep, self.tausfh)
		term1 /= self.taudep
		term1 *= (1 - m.exp(-time / self.harmonic_timescale(self.taudep,
			self.tausfh)))
		term2 = self.harmonic_timescale(self.taudep, self.taurise, self.tausfh)
		term2 /= self.taudep
		term2 *= m.exp(-time / self.taurise) - m.exp(
			-time / self.harmonic_timescale(self.taudep, self.tausfh))
		return term1 - term2


	def _gdot(self, time):
		term1 = m.exp(-time / self.harmonic_timescale(self.taudep, self.tausfh))
		term1 /= self.taudep
		term2 = self.harmonic_timescale(self.taudep, self.taurise, self.tausfh)
		term2 /= self.taudep
		term3 = m.exp(-time / self.taurise) / self.taurise
		term4 = m.exp(-time / self.harmonic_timescale(self.taudep, self.tausfh))
		term4 /= self.harmonic_timescale(self.taudep, self.tausfh)
		return term1 + term2 * (term3 - term4)


	def _gdotdot(self, time):
		term1 = -m.exp(-time / self.harmonic_timescale(self.taudep, self.tausfh))
		term1 /= self.taudep * self.harmonic_timescale(self.taudep, self.tausfh)
		term2 = self.harmonic_timescale(self.taudep, self.taurise, self.tausfh)
		term2 /= self.taudep
		term3 = m.exp(-time / self.taurise) / self.taurise**2
		term4 = m.exp(-time / self.harmonic_timescale(self.taudep, self.tausfh))
		term4 /= self.harmonic_timescale(self.taudep, self.tausfh)**2
		return term1 - term2 * (term3 - term4)


	@property
	def taudep(self):
		return self.taustar / (1 + self.eta - self.recycling)


class Zdotalpha(Zalpha):

	def __call__(self, time):
		if time == 0: return float("nan")
		prefactor = -m.exp(-time / self.taurise) / (self.taurise * (1 - m.exp(
			-time / self.taurise)))
		prefactor += self._gdot(time) / self._g(time)
		return prefactor * super().__call__(time)


class Zdotdotalpha(Zalpha):

	def __call__(self, time):
		if time == 0: return float("nan")
		term1 = m.exp(-time / self.taurise) + m.exp(-2 * time / self.taurise)
		term1 /= self.taurise**2 * (1 - m.exp(-time / self.taurise))**2
		term2 = 2 * m.exp(-time / self.taurise)
		term2 /= self.taurise * (1 - m.exp(-time / self.taurise))
		term2 *= self._gdot(time) / self._g(time)
		term3 = self._gdotdot(time) / self._g(time)
		prefactor = term1 - term2 + term3
		return prefactor * super().__call__(time)


































































