"""
Neutral-atom Mott-Bethe patch for phenix.refine
======================================================

Purpose
-------
This script lets you run a standard ``phenix.refine`` refinement but
replaces the internal electron scattering factor calculation with the
Mott-Bethe transform of neutral X-ray scattering factors.

Phenix's built-in ``electron`` scattering table uses tabulated electron
form factors.  This wrapper instead derives electron scattering factors
on-the-fly from X-ray form factors using the Mott-Bethe relation:

    f_e(s) = (C * gamma) / s^2  *  [ Z  -  f_x(s) ]

where:
    f_e(s)  = electron scattering factor (what we want)
    f_x(s)  = neutral-atom X-ray scattering factor (from wk1995 or it1992)
    Z       = atomic number (nuclear charge)
    s       = sin(theta) / lambda  =  1 / (2*d)
    C       = Mott-Bethe constant  =  m_e * e^2 / (2 * h^2)  ~  0.023934 A
    gamma   = relativistic Lorentz factor  =  1 + eV / (m_e * c^2)

At the structure-factor level this becomes:

    F_electron(h) = (C * gamma) / s(h)^2  *  [ F_nuclear(h)  -  F_xray(h) ]

where F_nuclear is the Fourier sum over bare nuclear charges Z_i and
F_xray is the ordinary neutral X-ray structure factor.

How it works
-----------------------------
The script patches two methods on internal Phenix/CCTBX classes:

1. ``f_model.manager.compute_f_calc``  —  replaced with
   ``_patched_compute_f_calc`` so that every F_calc evaluation during
   refinement returns Mott-Bethe electron structure factors instead of
   whatever the built-in table would give.

2. ``target_result_mixin.gradients_wrt_atomic_parameters``  —  replaced
   with ``_patched_gradients_wrt_atomic_parameters`` so that the
   analytical gradients used during coordinate, B-factor, and occupancy
   refinement are consistent with the Mott-Bethe F_calc.

Because both the forward model (F_calc) and the backward pass (gradients)
are patched, ``phenix.refine`` can minimise R-factors under the Mott-Bethe
electron model without any modification to the Phenix source tree.

Usage
-----
This script is invoked by the shell launcher ``bin/phenix.refine.mb``.
See ``HOWTO.md`` for full usage instructions.
"""

from __future__ import division, print_function

import argparse
import re
import sys

from cctbx import xray
from cctbx.array_family import flex
from cctbx.eltbx import tiny_pse          # Periodic System of Elements lookup
from cctbx.eltbx import xray_scattering   # X-ray form-factor tables
from iotbx.cli_parser import run_program
from libtbx.utils import user_plus_sys_time
from mmtbx.f_model import f_model
import mmtbx.refinement.targets
from phenix.programs.phenix_refine import (
  Program,
  PhenixRefineParser,
  custom_process_arguments,
)


# ---------------------------------------------------------------------------
#  Constants
# ---------------------------------------------------------------------------

# Mott-Bethe constant C = m_e * e^2 / (2 * h^2) in Angstrom units.
# This is the prefactor that converts (Z - f_x) / s^2 into an electron
# scattering factor in Angstroms.  Value: 0.023934 A.
MOTT_BETHE_COEFFICIENT = 0.023934

# Regex to extract the 1- or 2-letter element symbol from a CCTBX
# scattering-type string.  For example "Fe2+" -> "Fe", "Ca" -> "Ca".
_ELEMENT_RE = re.compile(r"^([A-Za-z]{1,2})")


# ---------------------------------------------------------------------------
#  Configuration container
# ---------------------------------------------------------------------------

class _Config(object):
  """Holds the user-chosen settings for the Mott-Bethe patch.

  Attributes
  ----------
  xray_table : str
      Which X-ray form-factor table to use for f_x(s).
      Either ``"wk1995"`` (Waasmaier & Kirfel 1995) or ``"it1992"``
      (International Tables 1992).
  electron_voltage_kv : float
      Accelerating voltage of the electron microscope in kilovolts.
      Used to compute the relativistic Lorentz factor gamma.
  verbose : bool
      If True, print diagnostic messages to stderr during refinement.
  """

  def __init__(self, xray_table, electron_voltage_kv, verbose):
    self.xray_table = xray_table
    self.electron_voltage_kv = float(electron_voltage_kv)
    self.verbose = bool(verbose)


# ---------------------------------------------------------------------------
#  Gradient helper
# ---------------------------------------------------------------------------

class _GradientDifference(object):
  """Proxy that presents (nuclear_gradients - xray_gradients) to Phenix.

  When ``phenix.refine`` asks for gradients of the target with respect to
  atomic parameters (coordinates, B-factors, occupancies), the Mott-Bethe
  chain rule requires computing gradients through both the nuclear and the
  X-ray structure-factor branches and then subtracting:

      d_target/d_params_electron  =  d_target/d_params_nuclear
                                    - d_target/d_params_xray

  (The outer (C*gamma)/s^2 scaling is already folded into the weighted
  derivative array before the gradient engine is called, so only the
  difference is needed here. I hope this is correct.)

  This class wraps the two raw gradient objects and implements the same
  accessor interface that Phenix expects (``packed()``,
  ``d_target_d_u_cart()``, etc.), returning the difference for each.
  """

  def __init__(self, nuclear_gradients, xray_gradients):
    self.nuclear_gradients = nuclear_gradients
    self.xray_gradients = xray_gradients

  def packed(self):
    """Return packed (site) gradient vector: nuclear minus xray."""
    return self.nuclear_gradients.packed() - self.xray_gradients.packed()

  def d_target_d_u_cart(self):
    """Return anisotropic-ADP gradient tensor: nuclear minus xray."""
    return self.nuclear_gradients.d_target_d_u_cart() - self.xray_gradients.d_target_d_u_cart()

  def d_target_d_u_iso(self):
    """Return isotropic-B gradient vector: nuclear minus xray."""
    return self.nuclear_gradients.d_target_d_u_iso() - self.xray_gradients.d_target_d_u_iso()

  def d_target_d_occupancy(self):
    """Return occupancy gradient vector: nuclear minus xray."""
    return self.nuclear_gradients.d_target_d_occupancy() - self.xray_gradients.d_target_d_occupancy()


# ---------------------------------------------------------------------------
#  Module-level state for the patch
# ---------------------------------------------------------------------------

# These globals store the patch configuration and the original (unpatched)
# methods so the patch can be installed exactly once and, in principle,
# reversed.  They are set by install_patch().
_PATCH_CONFIG = None               # _Config instance once patch is active
_ORIGINAL_COMPUTE_F_CALC = None    # saved f_model.manager.compute_f_calc
_ORIGINAL_GRADIENTS = None         # saved target_result_mixin.gradients_wrt_atomic_parameters


# ---------------------------------------------------------------------------
#  Physics helpers
# ---------------------------------------------------------------------------

def relativistic_gamma(voltage_kv):
  """Compute the relativistic Lorentz factor for an electron beam.

  Parameters
  ----------
  voltage_kv : float
      Accelerating voltage in kilovolts.

  Returns
  -------
  float
      gamma = 1 + eV / (m_e * c^2), where m_e * c^2 = 510.99895 keV.
      At 200 kV this is about 1.3914.

  Notes
  -----
  The factor corrects the Mott-Bethe prefactor for the increased
  electron mass at relativistic speeds.  Without it, the scattering
  factors would be too small at high voltages. In real world data this likely
  get eaten by scaling stepss. 
  """
  # Convert kV -> eV, then divide by electron rest energy in eV.
  return 1.0 + (float(voltage_kv) * 1000.0) / 510998.95


# ---------------------------------------------------------------------------
#  Element / scattering-type helpers
# ---------------------------------------------------------------------------

def _neutral_element_label(scattering_type):
  """Convert any CCTBX scattering-type string to a neutral element symbol.

  CCTBX scattering types can be bare elements ("C", "Fe"), charged ions
  ("Fe2+", "O1-"), or isotope labels ("D", "T").  This function strips
  everything back to the neutral element because the Mott-Bethe transform
  is defined for neutral atoms only. It can in priciple handle any arbitrary
  charge state, but implementing that takes a lot more work than this toy 
  implementation.

  Parameters
  ----------
  scattering_type : str or None
      The ``scatterer.scattering_type`` from a CCTBX xray_structure.

  Returns
  -------
  str
      Capitalised neutral element symbol, e.g. "Fe", "O", "H".

  Raises
  ------
  RuntimeError
      If the element cannot be identified or is not in the periodic table.

  Examples
  --------
  >>> _neutral_element_label("Fe2+")
  'Fe'
  >>> _neutral_element_label("D")
  'H'
  """
  label = (scattering_type or "").strip()

  # Try CCTBX's own label standardiser first (handles "FE2+" -> "Fe2+", etc.)
  standardized = xray_scattering.get_standard_label(
    label=label,
    exact=False,      # allow fuzzy matching
    optional=True,    # return None instead of raising if unrecognised
  )
  if standardized:
    label = standardized

  # Extract the 1- or 2-letter element symbol from the front of the string.
  match = _ELEMENT_RE.match(label)
  if match is None:
    raise RuntimeError(
      "Could not determine neutral element label from scattering type %r" % scattering_type
    )
  element = match.group(1).capitalize()   # "fe" -> "Fe", "C" -> "C"

  # Map deuterium and tritium to hydrogen.
  if element in {"D", "T"}:
    element = "H"

  # Validate against the CCTBX periodic table.
  if tiny_pse.table(element).atomic_number() <= 0:
    raise RuntimeError("Unknown element %r from scattering type %r" % (element, scattering_type))
  return element


# ---------------------------------------------------------------------------
#  Structure-building helpers
# ---------------------------------------------------------------------------

def _discard_scattering_registry(xray_structure):
  """Force an xray_structure to rebuild its scattering-type registry.

  After we change the ``scattering_type`` on individual scatterers, the
  cached registry is stale.  This function invalidates that cache so the
  next call to ``scattering_type_registry()`` rebuilds it from scratch.

  The public API method ``discard_scattering_type_registry()`` exists in
  newer CCTBX; the fallback directly pokes the private attributes for
  older versions.
  """
  if hasattr(xray_structure, "discard_scattering_type_registry"):
    xray_structure.discard_scattering_type_registry()
  else:
    # Fallback for older CCTBX versions without the public method.
    xray_structure._scattering_type_registry = None
    xray_structure._scattering_type_registry_is_out_of_date = True


def _build_xray_structure(xray_structure, table):
  """Build a neutral-atom X-ray structure for the f_x(s) term.

  Creates a deep copy of the input structure, strips all scattering types
  to neutral element symbols (e.g. "Fe2+" -> "Fe"), and assigns the
  chosen X-ray form-factor table (wk1995 or it1992).

  This structure is used to compute ``F_xray(h)`` — the neutral X-ray
  structure factor that appears in the Mott-Bethe relation.

  Parameters
  ----------
  xray_structure : cctbx.xray.structure
      The current model from phenix.refine.
  table : str
      X-ray form-factor table name ("wk1995" or "it1992").

  Returns
  -------
  cctbx.xray.structure
      Deep copy with neutral scattering types and the specified table.
  """
  result = xray_structure.deep_copy_scatterers()
  for scatterer in result.scatterers():
    # Strip charge labels -> neutral element for X-ray form factors.
    scatterer.scattering_type = _neutral_element_label(scatterer.scattering_type)
  _discard_scattering_registry(result)
  result.scattering_type_registry(table=table)
  return result


def _build_nuclear_structure(xray_structure):
  """Build the nuclear pseudo-structure for the F_nuclear(h) term.

  The Mott-Bethe relation needs the Fourier sum over bare nuclear charges:

      F_nuclear(h) = sum_i  Z_i * exp(-B_i * s^2) * exp(2*pi*i * h . r_i)

  CCTBX computes structure factors from Gaussian form-factor models.  A
  constant scatterer (no angular dependence) is modelled by a Gaussian
  with a single coefficient equal to Z and zero exponent, so f(s) = Z for
  all s.  We register one such Gaussian per unique atomic number.

  Parameters
  ----------
  xray_structure : cctbx.xray.structure
      The current model from phenix.refine.

  Returns
  -------
  cctbx.xray.structure
      Deep copy where every atom's scattering type has been replaced by
      a custom constant-Z Gaussian, so that the resulting structure
      factor gives F_nuclear(h).
  """
  result = xray_structure.deep_copy_scatterers()
  custom_dict = {}  # maps "Z_<n>" labels to constant Gaussian(Z)

  for scatterer in result.scatterers():
    element = _neutral_element_label(scatterer.scattering_type)
    atomic_number = tiny_pse.table(element).atomic_number()

    # Label like "Z_6" for carbon (Z=6), "Z_26" for iron, etc.
    label = "Z_%d" % atomic_number
    scatterer.scattering_type = label

    if label not in custom_dict:
      # xray_scattering.gaussian(float(Z)) creates a Gaussian with a
      # single term: coefficient = Z, exponent = 0.  This gives f(s) = Z
      # for all s — a constant "point charge" form factor.
      custom_dict[label] = xray_scattering.gaussian(float(atomic_number))

  _discard_scattering_registry(result)
  # Register the custom constant-Z form factors instead of any table.
  result.scattering_type_registry(custom_dict=custom_dict)
  return result


# ---------------------------------------------------------------------------
#  Structure-factor computation helpers
# ---------------------------------------------------------------------------

def _structure_factors_from_scatterers(miller_array, xray_structure, sfg_params):
  """Compute F_calc from an xray_structure using the same FFT/direct parameters
  that phenix.refine would normally use.

  This is a thin wrapper around CCTBX's ``structure_factors_from_scatterers``
  that forwards all the structure-factor-generation parameters (algorithm,
  grid resolution, etc.) from the current refinement session.

  Parameters
  ----------
  miller_array : cctbx.miller.array
      The set of Miller indices (h,k,l) to compute F_calc for.
  xray_structure : cctbx.xray.structure
      The atomic model (with appropriate form factors registered).
  sfg_params : object
      The ``sfg_params`` from ``f_model.manager`` containing algorithm
      settings (FFT vs direct summation, grid factors, etc.).

  Returns
  -------
  cctbx.miller.array
      Complex structure factors F_calc for each reflection.
  """
  manager = miller_array.structure_factors_from_scatterers(
    xray_structure=xray_structure,
    algorithm=sfg_params.algorithm,
    cos_sin_table=sfg_params.cos_sin_table,
    grid_resolution_factor=sfg_params.grid_resolution_factor,
    quality_factor=sfg_params.quality_factor,
    u_base=sfg_params.u_base,
    b_base=sfg_params.b_base,
    wing_cutoff=sfg_params.wing_cutoff,
    exp_table_one_over_step_size=sfg_params.exp_table_one_over_step_size,
  )
  return manager.f_calc()


def _mott_bethe_scales(miller_array, voltage_kv):
  """Compute the per-reflection Mott-Bethe scaling factors.

  For each reflection with d-spacing d(hkl), the Mott-Bethe weight is:

      scale(h) = (C * gamma) / s(h)^2

  where s = sin(theta)/lambda = 1/(2*d), so s^2 = 1/(4*d^2).

  Parameters
  ----------
  miller_array : cctbx.miller.array
      Miller array whose d-spacings define the s values.
  voltage_kv : float
      Electron accelerating voltage in kV (for the gamma factor).

  Returns
  -------
  prefactor : float
      The scalar C * gamma (without the 1/s^2 part).
  scales : flex.double
      Per-reflection array of C*gamma / s^2 values.

  Notes
  -----
  Two ``max(..., 1e-12)`` guards prevent division by zero:
  - The inner one protects against d ≈ 0 (would make s^2 infinite).
  - The outer one protects against s^2 ≈ 0 (very low angle, d → ∞),
    which would make the scale blow up.  In practice, normal
    crystallographic data never has d → 0 or d → ∞, but the guards
    are cheap insurance.
  """
  prefactor = MOTT_BETHE_COEFFICIENT * relativistic_gamma(voltage_kv)
  d_spacings = miller_array.d_spacings().data()
  scales = flex.double()
  for d_value in d_spacings:
    # s = 1/(2d), so s^2 = 1/(4*d^2).
    s_sq = 1.0 / (4.0 * max(float(d_value) * float(d_value), 1.0e-12))
    # scale = (C * gamma) / s^2  =  (C * gamma) * 4 * d^2.
    scales.append(prefactor / max(s_sq, 1.0e-12))
  return prefactor, scales


def _mb_weighted_derivative_array(d_target_d_f_calc_work, voltage_kv):
  """Multiply the target derivative array by Mott-Bethe per-reflection weights.

  During gradient back-propagation, the chain rule requires propagating
  d_target/d_F_calc through the Mott-Bethe scaling.  This function
  multiplies each element of the derivative array by the corresponding
  (C*gamma)/s^2 weight, producing the weighted derivative that is then
  passed to the nuclear and X-ray gradient engines.

  Parameters
  ----------
  d_target_d_f_calc_work : cctbx.miller.array
      Complex derivatives of the target function with respect to F_calc,
      as computed by the Phenix target machinery.
  voltage_kv : float
      Electron accelerating voltage in kV.

  Returns
  -------
  prefactor : float
      The scalar C * gamma.
  weighted : cctbx.miller.array
      Copy of the input with data scaled by the per-reflection MB weights.
  """
  prefactor, scales = _mott_bethe_scales(d_target_d_f_calc_work, voltage_kv)
  weighted = d_target_d_f_calc_work.customized_copy(
    data=d_target_d_f_calc_work.data() * scales
  )
  return prefactor, weighted


# ---------------------------------------------------------------------------
#  patch replacements
# ---------------------------------------------------------------------------

def _patched_compute_f_calc(self, miller_array=None, xray_structure=None):
  """Replacement for ``f_model.manager.compute_f_calc``.

  Instead of using Phenix's built-in electron form-factor table, this
  computes electron structure factors via the Mott-Bethe relation:

      F_e(h) = [C * gamma / s(h)^2]  *  [F_nuclear(h) - F_xray(h)]

  Steps:
    1. Build a neutral X-ray structure  ->  compute F_xray(h)
    2. Build a nuclear pseudo-structure ->  compute F_nuclear(h)
    3. Compute per-reflection MB weights (C*gamma / s^2)
    4. Return (F_nuclear - F_xray) * weights

  Parameters
  ----------
  self : f_model.manager
      The Phenix f-model manager instance (passed implicitly because
      this replaces a bound method).
  miller_array : cctbx.miller.array, optional
      Miller indices to use.  Defaults to ``self.f_obs()``.
  xray_structure : cctbx.xray.structure, optional
      Atomic model to use.  Defaults to ``self.xray_structure``.

  Returns
  -------
  cctbx.miller.array
      Complex Mott-Bethe electron structure factors for each reflection.
  """
  # Use the provided arguments or fall back to the manager's own data.
  xrs = xray_structure if xray_structure is not None else self.xray_structure
  miller_array = miller_array if miller_array is not None else self.f_obs()
  if miller_array.indices().size() == 0:
    raise RuntimeError("Empty miller_array.")

  # 1. Build the two parallel structures.
  xray_structure_mb = _build_xray_structure(xrs, _PATCH_CONFIG.xray_table)
  nuclear_structure_mb = _build_nuclear_structure(xrs)

  # 2. Compute F_xray and F_nuclear using the same FFT/direct parameters
  #    that phenix.refine would normally use.
  f_x = _structure_factors_from_scatterers(
    miller_array=miller_array,
    xray_structure=xray_structure_mb,
    sfg_params=self.sfg_params,
  )
  f_n = _structure_factors_from_scatterers(
    miller_array=miller_array,
    xray_structure=nuclear_structure_mb,
    sfg_params=self.sfg_params,
  )

  # 3. Compute per-reflection Mott-Bethe weights: (C * gamma) / s^2.
  prefactor, scales = _mott_bethe_scales(f_x, _PATCH_CONFIG.electron_voltage_kv)

  if _PATCH_CONFIG.verbose:
    print(
      "[phenix.refine.mb] compute_f_calc using Mott-Bethe neutral hack "
      "(table=%s, voltage=%.1f kV, prefactor=%.6f)"
      % (_PATCH_CONFIG.xray_table, _PATCH_CONFIG.electron_voltage_kv, prefactor),
      file=sys.stderr,
    )

  # 4. Apply the Mott-Bethe relation:
  #    F_electron(h) = scale(h) * [F_nuclear(h) - F_xray(h)]
  return f_x.customized_copy(data=(f_n.data() - f_x.data()) * scales)


def _patched_gradients_wrt_atomic_parameters(
  self,
  selection=None,
  site=False,
  u_iso=False,
  u_aniso=False,
  occupancy=False,
  tan_b_iso_max=None,
  u_iso_refinable_params=None,
):
  """Replacement for ``target_result_mixin.gradients_wrt_atomic_parameters``.

  Computes analytical gradients of the refinement target with respect to
  atomic parameters (coordinates, B-factors, occupancies) under the
  Mott-Bethe electron model.

  The chain rule for the Mott-Bethe relation gives:

      d_target       d_target   d_F_e     C*gamma   d_target
      -------- = ------------ * ----- = --------- * ---------- * (d_F_n/d_p - d_F_x/d_p)
        d_p      d_F_calc_work  d_p       s^2      d_F_calc_work

  So the procedure is:
    1. Multiply d_target/d_F_calc by the per-reflection MB weights
    2. Feed the weighted derivatives into the CCTBX gradient engine twice:
       once for the nuclear structure, once for the X-ray structure
    3. Return the difference (nuclear - xray)

  Parameters
  ----------
  self : target_result_mixin
      The Phenix target result object (passed implicitly).
  selection : flex.bool, optional
      Atom selection mask.
  site, u_iso, u_aniso, occupancy : bool
      Which gradient type Phenix is requesting.
  tan_b_iso_max : float, optional
      Not supported by this wrapper (must be None or 0).
  u_iso_refinable_params : flex.double, optional
      Pre-filtered isotropic B parameters.

  Returns
  -------
  varies
      - If ``u_aniso``: flex array of d_target/d_u_cart
      - If ``u_iso``: flex array of d_target/d_u_iso
      - If ``occupancy``: flex array of d_target/d_occupancy
      - Otherwise: a ``_GradientDifference`` proxy with ``packed()`` etc.
  """
  if tan_b_iso_max is not None and tan_b_iso_max != 0:
    raise RuntimeError("tan_b_iso_max is not supported by phenix.refine.mb")

  timer = user_plus_sys_time()
  manager = self.manager
  xray_structure = manager.xray_structure

  # Optionally restrict to a subset of atoms.
  if selection is not None:
    xray_structure = xray_structure.select(selection)

  # Step 1: Weight d_target/d_F_calc by Mott-Bethe per-reflection factors.
  _, weighted = _mb_weighted_derivative_array(
    self.d_target_d_f_calc_work(),
    _PATCH_CONFIG.electron_voltage_kv,
  )

  # Build the two parallel structures for the gradient computation.
  xray_structure_mb = _build_xray_structure(xray_structure, _PATCH_CONFIG.xray_table)
  nuclear_structure_mb = _build_nuclear_structure(xray_structure)

  # Common keyword arguments for the CCTBX gradient engine.
  kwargs = dict(
    u_iso_refinable_params=u_iso_refinable_params,
    d_target_d_f_calc=weighted.data(),
    n_parameters=xray_structure.n_parameters(),
    miller_set=weighted,
    algorithm=manager.sfg_params.algorithm,
  )

  if u_aniso or u_iso or occupancy:
    # When Phenix asks for a single gradient type (B-factor, occupancy,
    # or aniso ADP), we don't need the packed() site-gradient vector.
    # Setting n_parameters=0 and u_iso_refinable_params=None tells the
    # CCTBX gradient engine to skip building the packed array, which is
    # safe because we only access the individual accessors
    # (d_target_d_u_iso, d_target_d_u_cart, d_target_d_occupancy) that
    # are computed independently of n_parameters.
    kwargs["u_iso_refinable_params"] = None
    kwargs["n_parameters"] = 0

  # Step 2: Compute gradients through both branches.
  nuclear_gradients = manager.structure_factor_gradients_w(
    xray_structure=nuclear_structure_mb,
    **kwargs
  )
  xray_gradients = manager.structure_factor_gradients_w(
    xray_structure=xray_structure_mb,
    **kwargs
  )

  # Accumulate elapsed time into Phenix's global timer.
  mmtbx.refinement.targets.time_gradients_wrt_atomic_parameters += timer.elapsed()

  # Step 3: Return the difference (nuclear - xray) for the requested type.
  if u_aniso:
    return nuclear_gradients.d_target_d_u_cart() - xray_gradients.d_target_d_u_cart()
  if u_iso:
    return nuclear_gradients.d_target_d_u_iso() - xray_gradients.d_target_d_u_iso()
  if occupancy:
    return nuclear_gradients.d_target_d_occupancy() - xray_gradients.d_target_d_occupancy()
  # For site refinement (the default), return a proxy that computes the
  # difference lazily for whichever accessor Phenix calls.
  return _GradientDifference(nuclear_gradients, xray_gradients)


# ---------------------------------------------------------------------------
#  Patch installation
# ---------------------------------------------------------------------------

def install_patch(config):
  """Install the Mott-Bethe patch on the Phenix/CCTBX internals.

  Replaces two methods:
    - ``f_model.manager.compute_f_calc``
    - ``target_result_mixin.gradients_wrt_atomic_parameters``

  The originals are saved so the patch is installed at most once.

  Parameters
  ----------
  config : _Config
      User-chosen settings (X-ray table, voltage, verbosity).
  """
  global _PATCH_CONFIG
  global _ORIGINAL_COMPUTE_F_CALC
  global _ORIGINAL_GRADIENTS

  _PATCH_CONFIG = config

  # Save the original methods exactly once (guard against double-patching).
  if _ORIGINAL_COMPUTE_F_CALC is None:
    _ORIGINAL_COMPUTE_F_CALC = f_model.manager.compute_f_calc
  if _ORIGINAL_GRADIENTS is None:
    _ORIGINAL_GRADIENTS = mmtbx.refinement.targets.target_result_mixin.gradients_wrt_atomic_parameters

  # Swap in the patched versions.
  f_model.manager.compute_f_calc = _patched_compute_f_calc
  mmtbx.refinement.targets.target_result_mixin.gradients_wrt_atomic_parameters = _patched_gradients_wrt_atomic_parameters

  if _PATCH_CONFIG.verbose:
    print(
      "[phenix.refine.mb] installed neutral-atom Mott-Bethe patch "
      "(xray_table=%s, voltage=%.1f kV)"
      % (_PATCH_CONFIG.xray_table, _PATCH_CONFIG.electron_voltage_kv),
      file=sys.stderr,
    )


# ---------------------------------------------------------------------------
#  CLI argument parsing
# ---------------------------------------------------------------------------

def _parse_wrapper_args(argv):
  """Separate Mott-Bethe wrapper flags from phenix.refine arguments.

  All ``--mb-*`` flags belong to this wrapper.  Everything else is passed
  through unchanged to ``phenix.refine``.

  Parameters
  ----------
  argv : list of str
      Full command-line argument list (including argv[0] = script name).

  Returns
  -------
  options : argparse.Namespace
      Parsed wrapper options (mb_xray_table, mb_electron_voltage_kv,
      mb_verbose, mb_help).
  remaining : list of str
      argv[0] plus all non-wrapper arguments, ready to pass to Phenix.
  """
  parser = argparse.ArgumentParser(add_help=False)
  parser.add_argument(
    "--mb-xray-table",
    choices=["wk1995", "it1992"],
    default="wk1995",
    help="X-ray form-factor table for the neutral-atom f_x(s) term.",
  )
  parser.add_argument(
    "--mb-electron-voltage-kv",
    type=float,
    default=200.0,
    help="Electron accelerating voltage in kV (for relativistic gamma).",
  )
  parser.add_argument(
    "--mb-verbose",
    action="store_true",
    help="Print diagnostic messages from the MB patch to stderr.",
  )
  parser.add_argument(
    "--mb-help",
    action="store_true",
    help="Show this help for the MB wrapper and exit.",
  )

  # parse_known_args lets phenix.refine arguments pass through untouched.
  options, remaining = parser.parse_known_args(argv[1:])

  if options.mb_help:
    parser.print_help(sys.stderr)
    sys.stderr.write(
      "\nAll other arguments are passed directly to phenix.refine after the patch is installed.\n"
    )
    raise SystemExit(0)

  # Reconstruct argv with the script name at position 0.
  return options, [argv[0]] + remaining


# ---------------------------------------------------------------------------
#  Entry point
# ---------------------------------------------------------------------------

def main(argv=None):
  """Main entry point: parse wrapper args, install the patch, run phenix.refine.

  This is called by the shell launcher ``bin/phenix.refine.mb``.  It:
    1. Separates ``--mb-*`` flags from ordinary phenix.refine arguments.
    2. Installs the Mott-Bethe patch on the Phenix internals.
    3. Overwrites sys.argv with the remaining arguments.
    4. Calls ``phenix.refine``'s own ``run_program`` entry point, which
       proceeds as a normal refinement — except that every F_calc and
       gradient evaluation now uses the Mott-Bethe electron model.

  Parameters
  ----------
  argv : list of str, optional
      Command-line arguments.  Defaults to ``sys.argv``.

  Returns
  -------
  int
      Exit code (0 on success).
  """
  if argv is None:
    argv = list(sys.argv)

  options, remaining_argv = _parse_wrapper_args(argv)

  # Install the patch before phenix.refine starts.
  install_patch(
    _Config(
      xray_table=options.mb_xray_table,
      electron_voltage_kv=options.mb_electron_voltage_kv,
      verbose=options.mb_verbose,
    )
  )

  # Replace sys.argv so phenix.refine sees only its own arguments.
  sys.argv[:] = remaining_argv

  # Hand off to the standard phenix.refine entry point.
  run_program(
    program_class=Program,
    parser_class=PhenixRefineParser,
    custom_process_arguments=custom_process_arguments,
    unused_phil_raises_sorry=False,
  )
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
