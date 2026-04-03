# src/io.py

import numpy
import pandas

import awkward
import uproot

def get_root_file(
  FILE_PATH: str      
):
  '''
    Open a ROOT file via `uproot`.
  '''
  file = uproot.open(f"{FILE_PATH}")

  return file

def dump_sbruce_keys(
  FILE_PATH: str 
) -> None:
  '''
    Dump the keys of a ROOT file via `uproot`.
  '''
  file = uproot.open(f"{FILE_PATH}")
  print(file.keys())

  return

def convert_tree_to_df(
    FILE_PATH: str,
    DIRECTORY: str
) -> pandas.DataFrame:
  '''
    Convert a directory of a ROOT file to a `pandas.DataFrame`.
  '''
  file = uproot.open(f"{FILE_PATH}")
  arrays = file[f"{DIRECTORY}"].arrays(
      library = "np"
  )

  return pandas.DataFrame(arrays)

def convert_tree_to_df_from_file(
  file, 
  DIRECTORY: str
) -> pandas.DataFrame:
  '''
    Convert a directory of a ROOT file to a `pandas.DataFrame`,
    from an already-open file.
  '''
  arrays = file[DIRECTORY].arrays(
    library = "np"
  )
  return pandas.DataFrame(arrays)

def read_var_from_file(
  file,
  var: str,
):
  array = file['events/selectedNu'][var].array(
    library = "np"
  )

  return array

def read_spline_weights(
  file, 
  BRANCH: str
) -> numpy.ndarray:
  '''
    Read a multisigma systematic from a file.
    Each event stores a variable-length array of weights at sigma points
    [-3, -2, -1, 0, 1, 2, 3].
    Returns a `numpy.ndarray` of shape (N, 7).
  '''
  arr     = file['events/multisigmaTree'][BRANCH].array()
  arr_pad = awkward.pad_none(
    arr, 7, 
    clip = True
  )
  mat = awkward.to_numpy(
    arr_pad.to_numpy(allow_missing=True)
  ).astype(numpy.float64)

  mat = numpy.where(numpy.ma.getmaskarray(mat), 1.0, mat)

  return mat

def read_multisim_weights(
  file, 
  BRANCH: str, 
  N_UNIVERSES: int
) -> numpy.ndarray:
  '''
    Read a multisim systematic from a file.
    Pads with 1.0 when fewer universes are stored; truncates when more.
    Returns a `numpy.ndarray` of shape (N, n_universes), all non-negative values.
  '''
  arr     = file['events/multisimTree'][BRANCH].array()
  arr_pad = awkward.pad_none(
    arr, N_UNIVERSES, 
    clip = True
  )
  mat = awkward.to_numpy(
    arr_pad.to_numpy(allow_missing=True)
  ).astype(numpy.float64)

  mat = numpy.where(numpy.ma.getmaskarray(mat), 1.0, mat)

  return numpy.maximum(mat, 0.0)


def get_dfs_overlap(
    df_cv: pandas.DataFrame,
    df_var: pandas.DataFrame,
    key: str = 'trueE',
    common : bool = True,
):
  '''
    Get the (un)common instances between two `pandas.DataFrame`s,
    based on a specified `key`.
  '''
  common_events = pandas.merge(
      df_cv,
      df_var,
      on = key,
      how = 'inner'
  )

  common_in_cv = df_cv[key].isin(common_events[key])
  common_in_var1 = df_var[key].isin(common_events[key])

  if common:
    df_cv_out  = df_cv[common_in_cv]
    df_var_out = df_var[common_in_var1]
  else:
    df_cv_out  = df_cv[~common_in_cv]
    df_var_out = df_var[~common_in_var1]

  return df_cv_out, df_var_out