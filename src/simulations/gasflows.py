r"""
Handles radial gas flows in these models.
"""

from .._globals import MAX_SF_RADIUS
from vice.toolkit.interpolation import interp_scheme_1d
import numpy as np

class radial_gas_velocity_profile:

	def __init__(self, sigma_sfh, eta, sfe, N, beta_phi_in, beta_phi_out):
		self.sigma_sfh = sigma_sfh # fcn of radius and time
		self.eta = eta # function of radius
		self.sfe = sfe # function of time and sigma_sfh
		self.N = N # function of sigma_sfh and time
		self.beta_phi_in = beta_phi_in # function of radius
		self.beta_phi_out = beta_phi_out # function of radius

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
		dvdr -= (1 - 0.4) / self.sfe(time, sfr)
		dvdr += self.eta(radius) / self.sfe(time, sfr) * (
			self.beta_phi_out(radius) - self.beta_phi_in(radius)
		) / (self.beta_phi_in(radius) - 1)
		x = (self.beta_phi_in(radius) - 2)
		x /= radius * (self.beta_phi_in(radius) - 1)
		dsfr_dr = (1.e9 * self.sigma_sfh(radius + dr, time) - sfr) / dr
		x += 1 / self.N(time, sfr) * dsfr_dr / sfr
		dvdr -= vgas * x
		return dvdr


class radial_flow_driver(interp_scheme_1d):
	pass

class inward_flow_driver(radial_flow_driver):

	# return zero if the inferred velocity is positive
	def __call__(self, time):
		result = super().__call__(time)
		if result > 0:
			return 0
		else:
			return result

class outward_flow_driver(radial_flow_driver):

	# return zero if the inferred velocity is negative
	def __call__(self, time):
		result = super().__call__(time)
		if result < 0:
			return 0
		else:
			return result


