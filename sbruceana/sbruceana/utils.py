# src/utils.py

import numpy
import pandas
import uproot

def get_TH1(
  FILE_PATH: str,
  DIRECTORY: str
):
  '''
    Extract a ROOT `TH1` from the sbruce tree,
    in a pythonic way.
  '''
  file = uproot.open(f"{FILE_PATH}")
  hist = file[DIRECTORY] 

  values, edges = hist.to_numpy()

  return values, edges

def get_POT(
  FILE_PATH: str, 
):
  '''
    Get the POT associated with the sbruce tree.
    It is stored as a single-bin ROOT `TH1`.
  '''
  y, _ = get_TH1(
    FILE_PATH,
    "events/POT"
  )

  return y[0]

def get_livetime_offbeam(
  FILE_PATH: str, 
):
  '''
    Get the livetime associated with the sbruce tree from offbeam data.
    It is stored as a single-bin ROOT `TH1`.
  '''
  y, _ = get_TH1(
    FILE_PATH,
    "offbeam/Livetime"
  )

  return y[0]

def get_livetime_data(
  FILE_PATH: str, 
):
  '''
    Get the livetime associated with the sbruce tree from offbeam data.
    It is stored as a single-bin ROOT `TH1`.
  '''
  y, _ = get_TH1(
    FILE_PATH,
    "events/Livetime"
  )

  return y[0]

def get_offbeam_scale(
  livetime_offbeam,
  pot
):
  
  # POT-equivalent off-beam lifetime
  pot_offbeam = livetime_offbeam * 6.e13

  # scale factor from MC to offbeam in terms of POT
  MC_to_offbeam = pot_offbeam / pot

  return 1. / MC_to_offbeam

def get_ratio_of_vars(
  cv: pandas.DataFrame,
  pot_cv: float,
  var: pandas.DataFrame,
  pot_var: float,
  bins: numpy.array,
):
  '''
    Get a well-behaved `var`/`cv` ratio, in which the two are `pandas.Series`
    representing a single variable from the two `pandas.DataFrame`s. 
  '''
  y_cv, _      = numpy.histogram(cv, bins=bins)
  y_var, _     = numpy.histogram(var, bins=bins)
  y_var_scaled = y_var * (pot_cv / pot_var)

  with numpy.errstate(invalid = 'ignore'):
    ratio = numpy.where(
      y_cv > 0,
      y_var_scaled / y_cv,
      numpy.nan
    )
    ratio_err = numpy.where(
      y_cv > 0,
      (pot_cv / pot_var) * numpy.sqrt(y_var) / y_cv,
      numpy.nan
    )
  return ratio, ratio_err

def chi2_var(
    df_cv: pandas.DataFrame,
    df_var: pandas.DataFrame,
    bins: numpy.array,
    var: str,
    weight_cv: float = 1.0,
    weight_var: float = 1.0,
) -> tuple[float, float]:
    """
      Compute chi-squared between two weighted histograms using Poisson stat errors.
      Returns (chi2, ndof).
    """
    h_cv,  _ = numpy.histogram(df_cv[var],  bins=bins, weights=numpy.full(len(df_cv[var]),  weight_cv))
    h_var, _ = numpy.histogram(df_var[var], bins=bins, weights=numpy.full(len(df_var[var]), weight_var))

    # Poisson stat error on weighted histograms:
    # if each event carries weight w, sigma^2 = sum(w_i^2) = N * w^2
    n_cv,  _ = numpy.histogram(df_cv[var],  bins=bins)
    n_var, _ = numpy.histogram(df_var[var], bins=bins)
    sigma2_cv  = n_cv  * weight_cv**2
    sigma2_var = n_var * weight_var**2
    sigma2     = sigma2_cv + sigma2_var

    # only use bins where the combined variance is nonzero
    mask = sigma2 > 0
    chi2 = numpy.sum((h_cv[mask] - h_var[mask])**2 / sigma2[mask])
    ndof = mask.sum()

    return chi2, ndof