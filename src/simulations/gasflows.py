r"""
Handles radial gas flows in these models.
"""

from .._globals import MAX_SF_RADIUS
from vice.toolkit.interpolation import interp_scheme_1d
import numpy as np

class radial_gas_velocity_profile:

	def __init__(self, sigma_sfh, eta, sfe, N, beta_phi_in, beta_phi_out,
		outfilename = "gasvelocities.out"):
		self.sigma_sfh = sigma_sfh # fcn of radius and time
		self.eta = eta # function of radius and time
		self.sfe = sfe # function of time and sigma_sfh
		self.N = N # function of sigma_sfh and time
		self.beta_phi_in = beta_phi_in # function of radius and time
		self.beta_phi_out = beta_phi_out # function of radius and time
		if outfilename is not None:
			self.outfile = open(outfilename, 'w')
		else:
			self.outfile = None

	# def __call__(self, time, dr = 0.1, dt = 0.01):
	# 	dense_radii, dense_vgas = self.__call_sub__(time, dr = dr / 10,
	# 		dt = dt)
	# 	radii = dense_radii[::10]
	# 	vgas = dense_vgas[::10]
	# 	for i in range(len(radii)):
	# 		self.outfile.write("%.2e\t%.2e\t%.2e\n" % (radii[i], time, vgas[i]))
	# 	return [radii, vgas]
		
	def __call__(self, time, dr = 0.1, dt = 0.01):
		radii = [dr * i for i in range(int(20 / dr))]
		vgas = len(radii) * [0.]
		for i in range(1, len(radii)):
			if radii[i] <= MAX_SF_RADIUS:
				vgas[i] = vgas[i - 1] + dr * self.dvdr(time, radii[i - 1],
					vgas[i - 1], dr = dr, dt = dt)
			else:
				vgas[i] = 0
		if self.outfile is not None:
			for i in range(len(radii)):
				self.outfile.write("%.2e\t%.2e\t%.2e\n" % (
					radii[i], time, vgas[i]))
		else: pass
		return [radii, vgas] # vgas in kpc/Gyr

	def dvdr(self, time, radius, vgas, dr = 0.1, dt = 0.01):
		dvdr = 0
		sfr = self.sigma_sfh(radius + dr / 2, time) * 1.e9 # yr^-1 -> Gyr^-1
		if sfr <= 1.e-12: sfr = 1.e-12 # prevents ZeroDivisionError
		next_sfr = self.sigma_sfh(radius + dr / 2, time + dt) * 1.e9
		if next_sfr <= 1.e-12: next_sfr = 1.e-12
		dsfr_dt = (next_sfr - sfr) / dt
		dvdr -= dsfr_dt / (self.N(time, sfr) * sfr)
		dvdr -= (1 - 0.4) / self.sfe(time, sfr)
		dvdr += self.eta(radius, time) / self.sfe(time, sfr) * (
			self.beta_phi_in(radius, time) - self.beta_phi_out(radius, time)
		) / (1 - self.beta_phi_in(radius, time))
		next_sfr = self.sigma_sfh(radius + 3 * dr / 2, time) * 1.e9
		if next_sfr <= 1.e-12: next_sfr = 1.e-12
		dsfr_dr = (next_sfr - sfr) / dr
		if radius:
			x = (2 - self.beta_phi_in(radius, time))
			x /= radius * (1 - self.beta_phi_in(radius, time))
			x += dsfr_dr / (self.N(time, sfr) * sfr)
			dvdr -= vgas * x
		else:
			dvdr -= vgas * dsfr_dr / (self.N(time, sfr) * sfr)
			dvdr *= 1 - self.beta_phi_in(radius, time)
			dvdr /= 3 - 2 * self.beta_phi_in(radius, time)
		return dvdr


class radial_flow_driver(interp_scheme_1d):
	pass

# class inward_flow_driver(radial_flow_driver):

# 	# return zero if the inferred velocity is positive
# 	def __call__(self, time):
# 		result = super().__call__(time)
# 		if result > 0:
# 			return 0
# 		else:
# 			return result

# class outward_flow_driver(radial_flow_driver):

# 	# return zero if the inferred velocity is negative
# 	def __call__(self, time):
# 		result = super().__call__(time)
# 		if result < 0:
# 			return 0
# 		else:
# 			return result


