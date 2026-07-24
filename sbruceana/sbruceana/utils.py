# sbruceana/utils.py

import numpy
import pandas
import uproot
import scipy

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

def clopper_pearson(
    k, 
    n, 
    alpha=0.682
):
    """
        Return (efficiency, err_low, err_high) arrays using the Clopper-Pearson
        interval at coverage `alpha` (default 1-sigma, 68.2%).
        Bins with n==0 get eff=0 and symmetric zero errors.
    """
    eff     = numpy.where(n > 0, k / n, 0.)
    err_lo  = numpy.where(n > 0, eff - scipy.stats.beta.ppf((1 - alpha) / 2, k, n - k + 1), 0.)
    err_hi  = numpy.where(n > 0, scipy.stats.beta.ppf(1 - (1 - alpha) / 2, k + 1, n - k) - eff, 0.)

    # edge cases: k==0 or k==n
    err_lo  = numpy.where(k == 0, 0., err_lo)
    err_hi  = numpy.where(k == n, 0., err_hi)

    return eff, err_lo, err_hi

def augmented_bins(
    values, 
    bin_width   = 0.05, 
    offset      = 0.01, 
    range_min   = 0.0, 
    range_max   = 1.0
):

    values = numpy.asarray(values)
    starts = numpy.arange(range_min, range_max - bin_width + 1e-12, offset)

    bin_edges    = [(lo, lo + bin_width) for lo in starts]
    bin_contents = [values[(values >= lo) & (values < lo + bin_width)] for lo, _ in bin_edges]
    bin_counts   = numpy.array([len(b) for b in bin_contents])

    return bin_edges, bin_counts, bin_contents

def augmented_bins_with_breakpoints(
    values,
    bin_width  = 0.05,
    offset     = 0.01,
    range_min  = 0.0,
    range_max  = 1.0,
    bin_guides = None,
):
    """
    Collect values into overlapping sliding-window bins with optional
    variable bin widths and offsets above user-defined breakpoints.

    Parameters
    ----------
    values      : array-like of floats
    bin_width   : bin width used in the first segment (default 0.05)
    offset      : sliding-window step used in the first segment (default 0.01)
    range_min   : start of the first bin (default 0.0)
    range_max   : upper limit (default 1.0)
    bin_guides  : list of breakpoint specs, in any of these forms:
                    plain breakpoints  : [0.5, 0.75]
                      → width doubles, offset doubles at each step
                    (breakpoint, width): [(0.5, 0.10), (0.75, 0.25)]
                      → explicit width, offset doubles at each step
                    (breakpoint, width, offset): [(0.5, 0.10, 0.02), (0.75, 0.25, 0.05)]
                      → fully explicit per-segment width and offset

    Returns
    -------
    bin_edges   : list of (lo, hi) tuples
    bin_counts  : np.ndarray of counts per bin
    bin_contents: list of np.ndarray with the actual values per bin
    """
    values = numpy.asarray(values)

    segments = []  # (seg_lo, seg_hi, width, offset)

    if bin_guides is None:
        segments = [(range_min, range_max, bin_width, offset)]

    else:
        # normalise each guide to a (breakpoint, width, offset) triple
        guides = []
        for i, g in enumerate(bin_guides):
            if not hasattr(g, '__len__'):
                # plain scalar breakpoint: auto-double both width and offset
                guides.append((g,
                                bin_width * 2 ** (i + 1),
                                offset    * 2 ** (i + 1)))
            elif len(g) == 2:
                # (breakpoint, width): offset doubles automatically
                guides.append((g[0], g[1], offset * 2 ** (i + 1)))
            else:
                # (breakpoint, width, offset): fully specified
                guides.append((g[0], g[1], g[2]))

        guides = sorted(guides, key=lambda g: g[0])

        # first segment: base width and offset, up to the first breakpoint
        segments.append((range_min, guides[0][0], bin_width, offset))

        # intermediate segments
        for i, (bp, w, off) in enumerate(guides[:-1]):
            segments.append((bp, guides[i + 1][0], w, off))

        # last segment: to range_max
        bp, w, off = guides[-1]
        segments.append((bp, range_max, w, off))

    # build all bin start positions, using per-segment offset 
    starts_with_meta = []  # (start, width, offset) so we keep segment info
    for (seg_lo, seg_hi, w, off) in segments:
        s = numpy.arange(seg_lo, seg_hi - w + 1e-12, off)
        for lo in s:
            starts_with_meta.append((round(lo, 10), w, off))

    # deduplicate on start position, keeping first occurrence
    seen   = set()
    unique = []
    for item in sorted(starts_with_meta, key=lambda x: x[0]):
        if item[0] not in seen:
            seen.add(item[0])
            unique.append(item)

    # build edges, filter, collect 
    bin_edges    = [(lo, lo + w) for lo, w, _ in unique
                    if lo + w <= range_max + 1e-9]
    bin_contents = [values[(values >= lo) & (values < hi)]
                    for lo, hi in bin_edges]
    bin_counts   = numpy.array([len(b) for b in bin_contents])

    return bin_edges, bin_counts, bin_contents

def convert_ratio_TH1(
    name, 
    values, 
    edges, 
    filename, 
    mode = "recreate"
):
    """
      values: bin ratios (len = nbins)
      edges:  bin edges (len = nbins + 1)
    """
    with uproot.recreate(filename) if mode == "recreate" else uproot.update(filename) as f:
        f[name] = (values, edges)