r"""
Handles radial gas flows in these models.
"""

import numpy as np

class radial_gas_velocity:

	def __init__(self, sigma_sfh, eta, beta_phi, sfe, N):
		self.sigma_sfh = sigma_sfh # fcn of radius and time
		self.eta = eta # function of radius
		self.beta_phi = beta_phi # function of radius
		self.sfe = sfe # function of radius and time
		self.N = N # float (power-law index)

	def __call__(self, time, dr = 0.1, dt = 0.01):
		radii = [dr * i for i in range(int(20 / dr))]
		vgas = len(radii) * [0.]
		for i in range(1, len(radii)):
			vgas[i] = vgas[i - 1]
			sfr = self.sigma_sfh(radii[i], time)
			dsfr_dt = (self.sigma_sfh(radii[i], time + dt) - sfr) / dt
			term1 = dsfr_dt / (self.N * sfr)
			term2 = (1 + self.eta(radii[i]) - 0.4) / self.sfe(radii[i], time)
			term3 = (self.beta_phi(radii[i]) - 2) / (
				radii[i] * (self.beta_phi(radii[i]) - 1))
			dsfr_dr = (self.sigma_sfh(radii[i] + dr, time) - sfr) / dr
			term4 = dsfr_dr / (N * sfr)
			dv_dr = -term1 - term2 - vgas[i] * (term3 + term4)
			vgas[i] += dv_dr * dr

