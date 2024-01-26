r"""
The diskmodel objects employed in the Johnson et al. (2021) study.
"""

try:
	ModuleNotFoundError
except NameError:
	ModuleNotFoundError = ImportError
try:
	import vice
except (ModuleNotFoundError, ImportError):
	raise ModuleNotFoundError("Could not import VICE.")
if vice.version[:2] < (1, 2):
	raise RuntimeError("""VICE version >= 1.2.0 is required to produce \
Johnson et al. (2021) figures. Current: %s""" % (vice.__version__))
else: pass
# from vice.yields.presets import JW20
from . import yields
from .sfe import sfe, sfe_oscil
from .gasflows import radial_gas_velocity_profile
from .gasflows import inward_flow_driver, outward_flow_driver
from vice.toolkit import hydrodisk
# from vice.toolkit.interpolation import interp_scheme_1d
# vice.yields.sneia.settings['fe'] *= 10**0.1
from .._globals import END_TIME, MAX_SF_RADIUS, ZONE_WIDTH
from . import migration
from . import models
from .models.utils import get_bin_number, interpolate, gaussian, skewnormal
from .models.gradient import gradient
import warnings
import math as m
import sys
import os

_SECONDS_PER_GYR_ = 3.1536e16
_KPC_PER_KM_ = 3.24e-17


def eta_function(radius):
	# etasun = vice.yields.ccsne.settings['o'] / vice.solar_z['o'] - 0.6
	# return etasun * m.exp(0.062 * m.log(10) * (radius - 8))
	# return 0
	return 1

def molecular_tau_star(time):
	return 2 * ((0.5 + time) / 13.7)**0.5

def sfe_function(time, sigma_sfr):
	mol = molecular_tau_star(time)
	N = plaw_index(time, sigma_sfr)
	return mol * (sigma_sfr * mol / 1e8)**(1 / N - 1)

def plaw_index(time, sigma_sfr):
	mol = molecular_tau_star(time)
	if sigma_sfr > 1e8 / mol:
		return 1
	else:
		return 1.5


class diskmodel(vice.milkyway):

	r"""
	A milkyway object tuned to the Johnson et al. (2021) models specifically.

	Parameters
	----------
	zone_width : ``float`` [default : 0.1]
		The width of each annulus in kpc.
	name : ``str`` [default : "diskmodel"]
		The name of the model; the output will be stored in a directory under
		this name with a ".vice" extension.
	spec : ``str`` [default : "static"]
		A keyword denoting the time-dependence of the star formation history.
		Allowed values:

		- "static"
		- "insideout"
		- "lateburst"
		- "outerburst"

	verbose : ``bool`` [default : True]
		Whether or not the run the models with verbose output.
	migration_mode : ``str`` [default : "diffusion"]
		A keyword denoting the time-dependence of stellar migration.
		Allowed values:

		- "diffusion"
		- "linear"
		- "sudden"
		- "post-process"

	kwargs : varying types
		Other keyword arguments to pass ``vice.milkyway``.

	Attributes and functionality are inherited from ``vice.milkyway``.
	"""

	def __init__(self, zone_width = 0.1, name = "diskmodel", spec = "static",
		verbose = True, migration_mode = "diffusion",
		radial_gas_flows = False, dt = 0.01, **kwargs):
		super().__init__(zone_width = zone_width, name = name,
			verbose = verbose, **kwargs)
		self.dt = dt
		if self.zone_width <= 0.2 and self.dt <= 0.02 and self.n_stars >= 6:
			Nstars = 3102519
		else:
			Nstars = 2 * int(MAX_SF_RADIUS / zone_width * END_TIME / self.dt *
				self.n_stars)

		# self.migration.stars = migration.diskmigration(self.annuli,
		# 	N = Nstars, mode = migration_mode,
		# 	filename = "%s_analogdata.out" % (name))
			# for elem in self.zones[0].elements:
			# 	vice.yields.ccsne.settings['o'] /= 2
			# 	vice.yields.sneia.settings['o'] /= 2

		# for i in range(self.n_zones):
			# self.zones[i].eta = 0
			# radius = ZONE_WIDTH * (i + 0.5)
			# eta = vice.yields.ccsne.settings['o'] / vice.solar_z['o']
			# eta *= m.exp(0.06 * (radius - 8) * m.log(10))
			# eta -= 0.6
			# if eta > 0:
			# 	self.zones[i].eta = eta
			# else:
			# 	self.zones[i].eta = 0
		for i in range(self.n_zones):
			self.zones[i].eta = eta_function(zone_width * (i + 0.5))

		for i in range(self.n_zones):
			rmin = zone_width * i
			area = m.pi * ((rmin + zone_width)**2 - rmin**2)
			if spec == "SFEoscil":
				spec = "insideout"
				self.zones[i].tau_star = sfe_oscil(area, mode = "sfr",
					amplitude = 0.5, period = 2)
			else:
				self.zones[i].tau_star = sfe(area, mode = "sfr")

		self.migration.stars = migration.gaussian_migration(self.annuli,
			zone_width = zone_width,
			filename = "%s_analogdata.out" % (self.name),
			post_process = self.simple)
		self.evolution = star_formation_history(spec = spec,
			zone_width = zone_width, timestep = self.zones[0].dt)
		self.mode = "sfr"

		if radial_gas_flows:
			vgas_engine = radial_gas_velocity_profile(
				self.evolution, eta_function, sfe_function, plaw_index,
				lambda r: 0.7, lambda r: 0.0)
			vgas_all = []
			times = [self.dt * i for i in range(int(END_TIME / self.dt) + 10)]
			for i in range(len(times)):
				radii, vgas = vgas_engine(i * self.dt) # in kpc/Gyr
				vgas_all.append(vgas)
			matrix_elements_inward = []
			matrix_elements_outward = []
			for i in range(self.n_zones):
				yvals = []
				vgas = [row[i] for row in vgas_all]
				for j in range(len(times)):
					# normalized to a 10 Myr time interval
					# don't turn on the flow until this many timesteps have passed
					if j > 10:
						radius = i * zone_width
						numerator = vgas[j]**2 * 0.01**2
						numerator -= 2 * radius * vgas[j] * 0.01
						denominator = 2 * radius * zone_width + zone_width**2
						areafrac = numerator / denominator
						if areafrac * self.dt / 0.01 > 1:
							warnings.warn("""\
Area fraction larger than 1. Consider comparing results with different \
timestep sizes to assess the impact of numerical artifacts.""")
							areafrac = 0.01 / self.dt - 1.e-9
						elif areafrac < 0:
							areafrac = 1e-9
						else: pass
						yvals.append(areafrac)
						# if yvals[-1] * self.dt / 0.01 > 1: os.system(
						# 	"echo \"%d %d %.5f %.2e\" >> err.out" % (
						# 		i, j, yvals[-1], vgas[j]))
					else:
						yvals.append(0)
				matrix_elements_inward.append(inward_flow_driver(times, yvals))
				matrix_elements_outward.append(outward_flow_driver(times, yvals))
			for i in range(self.n_zones):
				for j in range(self.n_zones):
					if i - 1 == j:
						# inward gas flows
						self.migration.gas[i][j] = matrix_elements_inward[i]
					elif i + 1 == j:
						# outward gas flows
						self.migration.gas[i][j] = matrix_elements_outward[i]
					else:
						self.migration.gas[i][j] = 0
					# if i - 1 == j:
					# 	self.migration.gas[i][j] = gas_matrix_elements[i]
					# else:
					# 	self.migration.gas[i][j] = 0
		else:
			for i in range(self.n_zones):
				for j in range(self.n_zones):
					self.migration.gas[i][j] = 0


		# if spec == "simple-exponential" or spec == "subequilibrium":
		# 	for i in range(self.n_zones):
		# 		self.zones[i].eta = 0
		# 		radius = ZONE_WIDTH * (i + 0.5)
		# 		# print(radius)
		# 		self.zones[i].tau_star = 2 * m.exp(radius / 6)
		# 		self.zones[i].schmidt = True
		# 		self.zones[i].MgSchmidt = (self.zones[i].tau_star *
		# 			self.zones[i].func.surface_density._evol[i].norm / m.e *
		# 			1.e9)

		# for i in range(self.n_zones):
		# 	inner = i * zone_width
		# 	self.zones[i].tau_star = 2 * m.exp((inner + zone_width / 2) / 8)
		# 	self.zones[i].schmidt = True
		# 	self.zones[i].MgSchmidt = (1.e9 * self.zones[i].tau_star *
		# 		self.zones[i].func.surface_density._evol[i].norm / m.e *
		# 		m.pi * ((inner + zone_width)**2 - inner**2))
		# 	eta = (vice.yields.ccsne.settings['o'] /
		# 		vice.solar_z['o'] * 10**(0.06 * (inner + zone_width / 2 - 8))
		# 		- 0.6 + self.zones[i].tau_star /
		# 		self.zones[i].func.surface_density._evol[i].timescale)
		# 	self.zones[i].eta = eta if eta >= 0 else 0
		# if spec == "subequilibrium":
		# 	for i in range(self.n_zones): self.zones[i].eta = 0
		# 	for elem in self.zones[0].elements:
		# 		vice.yields.ccsne.settings[elem] /= 3
		# 		vice.yields.sneia.settings[elem] /= 3
		# else: pass
		
		# CONSTANT GAS VELOCITY
		# radial_gas_velocity *= _SECONDS_PER_GYR_
		# radial_gas_velocity *= _KPC_PER_KM_ # vrad now in kpc / Gyr
		# for i in range(self.n_zones):
		# 	for j in range(self.n_zones):
		# 		if i - 1 == j:
		# 			# normalized to 10 Myr time interval
		# 			numerator = radial_gas_velocity**2 * 0.01**2
		# 			numerator -= 2 * i * zone_width * radial_gas_velocity * 0.01
		# 			denominator = zone_width**2 * (2 * i + 1)
		# 			self.migration.gas[i][j] = numerator / denominator
		# 		else:
		# 			self.migration.gas[i][j] = 0

		# # VARIABLE GAS VELOCITY
		# radial_gas_velocity *= _SECONDS_PER_GYR_
		# radial_gas_velocity *= _KPC_PER_KM_ # vrad now in kpc / Gyr
		# for i in range(self.n_zones):
		# 	for j in range(self.n_zones):
		# 		if i - 1 == j:
		# 			# normalized to 10 Myr time interval
		# 			self.migration.gas[i][j] = variable_gas_velocity_migration(
		# 				i * zone_width, radial_gas_velocity)
		# 		else:
		# 			self.migration.gas[i][j] = 0




	def run(self, *args, **kwargs):
		out = super().run(*args, **kwargs)
		self.migration.stars.close_file()
		return out


	@classmethod
	def from_config(cls, config, **kwargs):
		r"""
		Obtain a ``diskmodel`` object with the parameters encoded into a
		``config`` object.

		**Signature**: diskmodel.from_config(config, **kwargs)

		Parameters
		----------
		config : ``config``
			The ``config`` object with the parameters encoded as attributes.
			See src/simulations/config.py.
		**kwargs : varying types
			Additional keyword arguments to pass to ``diskmodel.__init__``.

		Returns
		-------
		model : ``diskmodel``
			The ``diskmodel`` object with the proper settings.
		"""
		model = cls(zone_width = config.zone_width, dt = config.timestep_size,
			**kwargs)
		model.n_stars = config.star_particle_density
		model.bins = config.bins
		model.elements = config.elements
		model.setup_nthreads = config.setup_nthreads
		model.nthreads = config.nthreads
		return model


class variable_gas_velocity_migration(skewnormal):

	def __init__(self, radius, baseline, zone_width = 0.1,
		mean = 7, amplitude = 1, std = 1.5, skewness = 3):
		self.radius = radius
		self.baseline = baseline
		self.zone_width = zone_width
		super().__init__(mean = mean, amplitude = amplitude, std = std,
			skewness = skewness)

	def __call__(self, time):
		# normalized to a 10 Myr time interval
		vgas = self.baseline * (1 + super().__call__(time))
		numerator = vgas**2 * 0.01**2
		numerator -= 2 * self.radius * vgas * 0.01
		denominator = (self.radius + self.zone_width)**2 - self.radius**2
		return numerator / denominator


class star_formation_history:

	r"""
	The star formation history (SFH) of the model galaxy. This object will be
	used as the ``evolution`` attribute of the ``diskmodel``.

	Parameters
	----------
	spec : ``str`` [default : "static"]
		A keyword denoting the time-dependence of the SFH.
	zone_width : ``float`` [default : 0.1]
		The width of each annulus in kpc.

	Calling
	-------
	- Parameters

		radius : ``float``
			Galactocentric radius in kpc.
		time : ``float``
			Simulation time in Gyr.
	"""

	def __init__(self, spec = "static", zone_width = 0.1, timestep = 0.01):
		self._radii = []
		self._evol = []
		i = 0
		max_radius = 20 # kpc, defined by ``vice.milkyway`` object.
		while (i + 1) * zone_width <= max_radius:
			self._radii.append((i + 0.5) * zone_width)
			self._evol.append({
					"static": 				models.static,
					"insideout": 			models.insideout,
					"sfroscil": 			models.SFRoscil,
					"lateburst": 			models.lateburst,
					"outerburst": 			models.outerburst,
					"subequilibrium":		models.subequilibrium,
					"simple-exponential":	models.simple_exponential
				}[spec.lower()]((i + 0.5) * zone_width, dr = zone_width,
					dt = timestep))
			i += 1

	def __call__(self, radius, time):
		# The milkyway object will always call this with a radius in the
		# self._radii array, but this ensures a continuous function of radius
		if radius > MAX_SF_RADIUS:
			return 0
		else:
			idx = get_bin_number(self._radii, radius)
			if idx != -1:
				return gradient(radius) * interpolate(self._radii[idx],
					self._evol[idx](time), self._radii[idx + 1],
					self._evol[idx + 1](time), radius)
			else:
				return gradient(radius) * interpolate(self._radii[-2],
					self._evol[-2](time), self._radii[-1], self._evol[-1](time),
					radius)

