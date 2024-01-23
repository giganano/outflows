r"""
Handles radial gas flows in these models.
"""

from .._globals import MAX_SF_RADIUS
import numpy as np

class radial_gas_velocity_profile:

	def __init__(self, sigma_sfh, eta, beta_phi, sfe, N):
		self.sigma_sfh = sigma_sfh # fcn of radius and time
		self.eta = eta # function of radius
		self.beta_phi = beta_phi # function of radius
		self.sfe = sfe # function of time and sigma_sfh
		self.N = N # function of sigma_sfh and time

	def __call__(self, time, dr = 0.1, dt = 0.01):
		radii = [dr * i for i in range(int(20 / dr))]
		vgas = len(radii) * [0.]
		for i in range(1, len(radii)):
			if radii[i] <= MAX_SF_RADIUS:
				vgas[i] = vgas[i - 1] + dr * self.dvdr(time, radii[i], vgas[i - 1],
					dr = dr, dt = dt)
			else:
				vgas[i] = 0
		return [radii, vgas] # vgas in kpc/Gyr

	def dvdr(self, time, radius, vgas, dr = 0.1, dt = 0.01):
		dvdr = 0
		sfr = self.sigma_sfh(radius, time) * 1.e9 # yr^-1 -> Gyr^-1
		if sfr <= 1.e-12: sfr = 1.e-12 # prevents ZeroDivisionError
		dsfr_dt = (1.e9 * self.sigma_sfh(radius, time + dt) - sfr) / dt
		dvdr -= 1 / self.N(time, sfr) * dsfr_dt / sfr
		dvdr -= (1 + self.eta(radius) - 0.4) / self.sfe(time, sfr)
		x = (self.beta_phi(radius) - 2)
		x /= radius * (self.beta_phi(radius) - 1)
		dsfr_dr = (1.e9 * self.sigma_sfh(radius + dr, time) - sfr) / dr
		x += 1 / self.N(time, sfr) * dsfr_dr / sfr
		dvdr -= vgas * x
		return dvdr


