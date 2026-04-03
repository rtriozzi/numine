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

def chi2_histograms(
  df       : pandas.DataFrame,
  df_data  : pandas.DataFrame,
  var      : str,
  bins     : numpy.array,
  factor   : float = 1.
):
  '''
  Get the chi-sq between two histograms.
  If specified, the data df can be calibrated by a `factor`.
  '''
  mc_counts,   _ = numpy.histogram(df[var],               bins=bins)
  data_counts, _ = numpy.histogram(df_data[var] * factor, bins=bins)

  N_mc   = mc_counts.sum()
  N_data = data_counts.sum()

  mc_norm   = mc_counts / N_mc
  data_norm = data_counts / N_data

  # Only compare bins where both have entries
  mask = (mc_counts > 0) & (data_counts > 0)

  # Scale MC to data normalisation, then compare in count space
  mc_scaled = mc_norm[mask] * N_data
  observed  = data_counts[mask]

  chi2_val = numpy.sum((observed - mc_scaled) ** 2 / observed)
  ndof     = mask.sum()

  return chi2_val, ndof