# src/systs.py

import scipy
import numpy

from .io import read_var_from_file, read_multisim_weights, read_spline_weights
from .config import SYSTS_CONFIG

def _get_weighted_histograms_from_cv(
  values: numpy.ndarray,
  weights_per_event: numpy.ndarray,
  bin_edges: numpy.ndarray,
) -> numpy.ndarray:
  '''
    Returns one histogram per universe.
  '''
  U = weights_per_event.shape[1]
  B = len(bin_edges) - 1
  out = numpy.zeros((U, B), dtype=numpy.float64)

  for u in range(U):
    out[u], _ = numpy.histogram(
      values, 
      bins = bin_edges,
      weights = weights_per_event[:, u]
    )
      
  return out

def _get_cov_from_universes(
  cv: numpy.ndarray, 
  hist_universes: numpy.ndarray
) -> numpy.ndarray:
  U     = hist_universes.shape[0]
  delta = hist_universes - cv[numpy.newaxis, :] 

  return (delta.T @ delta) / U

def _convert_spline_to_universes(
  wgt_array: numpy.ndarray,
  N_UNIVERSES: int,
  RNG: numpy.random.Generator,
) -> numpy.ndarray:
  '''
    Convert spline weights from (N, 7) to (N, N_UNIVERSES), by converting
    a cubic spline through the 7 sigma points for each event, drawing
    N_UNIVERSES Gaussian sigma values, and evaluating the spline at each drawn sigma.
  '''

  N = wgt_array.shape[0]
  drawn_sigmas = RNG.normal(0.0, 1.0, size=N_UNIVERSES)
  universes = numpy.ones((N, N_UNIVERSES), dtype=numpy.float64)

  for i in range(N):
    w = wgt_array[i]

    # skip all-1 universes
    if numpy.allclose(w, 1.0):
        continue
    
    spl = scipy.interpolate.interp1d(
        numpy.array([-3., -2., -1., 0., 1., 2., 3.]), 
        w,
        kind = "linear",
        bounds_error = False,
        fill_value = (w[0], w[-1]),
    )

    universes[i] = numpy.maximum(spl(drawn_sigmas), 0.0)

  return universes
 
def _get_cov_norm(
  cv: numpy.ndarray, 
  frac_unc: float
) -> numpy.ndarray:
  '''
    Covariance for a flat fractional normalization (diagonal term).
  '''
  return (frac_unc ** 2) * numpy.outer(cv, cv)
 
def build_covariance_matrix(
  file,
  BINNING: numpy.ndarray,
  var: str = "recoE",
  SYSTS_CONFIG: dict = SYSTS_CONFIG,
  N_UNIVERSES: int = 500,
  SEED: int = 42,
):
  # define binning, seed
  BINNING = numpy.asarray(BINNING)
  x       = 0.5 * (BINNING[:-1] + BINNING[1:])
  rng     = numpy.random.default_rng(SEED)

  # get unweighted CV histogram
  values  = read_var_from_file(file, var)
  cv, _   = numpy.histogram(values, bins=BINNING)
  cv      = cv.astype(numpy.float64)

  results = {}

  # multisim systematics
  for syst in SYSTS_CONFIG.get("multisim", []):
    name   = syst["name"]
    branch = syst.get("branch", name)
    tag    = syst.get("tag", "other")

    try:
      wgt = read_multisim_weights(file, branch, N_UNIVERSES)
    except Exception as e:
      print(f"[build_covariance_matrix] cannot read multisim branch '{branch}': {e}")
      continue
    hists = _get_weighted_histograms_from_cv(values, wgt, BINNING)

    results[name] = {"cov": _get_cov_from_universes(cv, hists), "tag": tag}

  # spline systematics
  for syst in SYSTS_CONFIG.get("spline", []):
    name   = syst["name"]
    branch = syst.get("branch", name)
    tag    = syst.get("tag", "other")

    try:
      wgt = read_spline_weights(file, branch)
    except Exception as e:
      print(f"[build_covariance_matrix] cannot read spline branch '{branch}': {e}")
      continue
    
    hists = _get_weighted_histograms_from_cv(
      values, 
      _convert_spline_to_universes(wgt, N_UNIVERSES, rng), 
      BINNING
    )

    results[name] = {"cov": _get_cov_from_universes(cv, hists), "tag": tag}

  # normalization systematics
  for syst in SYSTS_CONFIG.get("norm", []):
      name = syst["name"]
      results[name] = {"cov": _get_cov_norm(cv, syst["value"]), "tag": syst.get("tag", "other")}

  return results, cv, x

def get_error_band(
  file,
  BINNING: numpy.ndarray,
  variable: str = "recoE",
  SYSTS_CONFIG: dict = SYSTS_CONFIG,
  groups = None,
  N_UNIVERSES: int = 500,
  seed: int = 42,
):
  BINNING   = numpy.asarray(BINNING)
  bin_centers = 0.5 * (BINNING[:-1] + BINNING[1:])

  values    = read_var_from_file(file, variable)
  cv, _   = numpy.histogram(values, bins=BINNING)
  cv      = cv.astype(numpy.float64)
  cov_total = numpy.zeros((len(bin_centers),) * 2, dtype=numpy.float64)

  include_all  = (groups is None)
  include_stat = include_all or "stat" in (groups or [])
  syst_tags    = set(groups or []) - {"stat"}
  need_syst    = include_all or bool(syst_tags)

  if include_stat:
      cov_total += numpy.diag(cv)

  if need_syst:
      cov_results, _, _ = build_covariance_matrix(
          file, BINNING, SYSTS_CONFIG, variable, N_UNIVERSES, seed
      )
      for entry in cov_results.values():
          if include_all or entry["tag"] in syst_tags:
              cov_total += entry["cov"]

  err = numpy.sqrt(numpy.maximum(numpy.diag(cov_total), 0.0))

  return bin_centers, cv, err, err, cov_total

def get_fractional_uncertainty(
  cov_results: dict,
  cv: numpy.ndarray,
):
  safe_cv = numpy.where(cv > 0, cv, numpy.nan)

  def _frac(cov):
      return numpy.sqrt(numpy.maximum(numpy.diag(cov), 0.0)) / safe_cv

  B            = len(cv)
  tags         = sorted({v["tag"] for v in cov_results.values()})
  cov_by_group = {tag: numpy.zeros((B, B)) for tag in tags}
  for entry in cov_results.values():
      cov_by_group[entry["tag"]] += entry["cov"]
  cov_by_group["stat"] = numpy.diag(cv)

  frac_by_group = {tag: _frac(cov) for tag, cov in cov_by_group.items()}
  frac_total    = _frac(sum(cov_by_group.values()))

  return frac_by_group, frac_total