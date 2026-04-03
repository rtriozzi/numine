# src/plotting.py

import numpy
import pandas

import matplotlib
import matplotlib.pyplot as plt

from .config import Category

def plot_var(
    ax,
    df: pandas.DataFrame,
    bins: numpy.array,
    var: str,
    label: str,
    weight: float = 1.0,
):
  ax.hist(
    df[var],
    bins = bins,
    label = label,
    histtype = 'step',
    linewidth = 2,
    weights = numpy.full(len(df[var]), weight),
  )

  return ax

def plot_var_by_category(
    ax,
    df: pandas.DataFrame,
    bins: numpy.array,
    var: str,
    band: bool = True,
    **kwargs,
):
  mask = (df.CC == 1) & (abs(df.truePDG) == 12)
  ax.hist(
    [df[mask][var], df[~mask][var]],
    stacked = True,
    bins    = bins,
    label   = ['$\\nu_e$CC', 'other']
  )

  counts, _ = numpy.histogram(df[var], bins=bins)
  errors    = numpy.sqrt(counts)

  x = 0.5 * (bins[:-1] + bins[1:])
  w = numpy.diff(bins)

  if band:
    ax.bar(
      x,
      2 * errors,
      width     = w,
      bottom    = counts - errors,
      fill      = True,
      linewidth = 0,
      color     = 'gray',
      alpha     = 0.5,
      **kwargs,
    )
  else:
    ax.errorbar(
      x,
      counts,
      yerr = errors,
      marker = '.',
      ls = '',
      c = 'black',
      capsize = 0
    )

  return ax

def plot_error_band(
  ax,
  cov_results: dict,
  cv: numpy.ndarray,
  bin_centers: numpy.ndarray,
  bin_edges: numpy.ndarray,
  groups = None,
  **kwargs,
):
  include_all = (groups is None)
  syst_tags   = set(groups or [])

  B = len(bin_centers)
  cov_total = numpy.zeros((B, B), dtype=numpy.float64)
  for entry in cov_results.values():
      if include_all or entry["tag"] in syst_tags:
          cov_total += entry["cov"]

  err   = numpy.sqrt(numpy.maximum(numpy.diag(cov_total), 0.0))
  width = numpy.diff(bin_edges) 

  bar_kwargs = dict(color="gray", alpha=0.4, linewidth=0)
  bar_kwargs.update(kwargs)

  ax.bar(
    bin_centers,
    2 * err,
    width  = width,
    bottom = cv - err,
    align  = "center",
    **bar_kwargs,
  )

  return ax

def plot_covariance_matrix(
  ax,
  cov_results: dict,
  cv: numpy.ndarray,
  bin_centers: numpy.ndarray,
  groups = None,
  label: str = '',
):
  include_all = (groups is None)
  syst_tags   = set(groups or [])

  B = len(bin_centers)
  cov_total = numpy.zeros((B, B), dtype=numpy.float64)
  for entry in cov_results.values():
    if include_all or entry["tag"] in syst_tags:
      cov_total += entry["cov"]

  # symmetric color scale centred on zero
  vmax = numpy.abs(cov_total).max()

  im = ax.imshow(
    cov_total,
    origin    = "lower",
    aspect    = "auto",
    cmap      = "coolwarm",
    vmin      = -vmax,
    vmax      =  vmax,
    extent    = [bin_centers[0],  bin_centers[-1], bin_centers[0],  bin_centers[-1]],
  )

  return ax, im

def plot_by_category(
    ax,
    df: pandas.DataFrame,
    categories: list[Category],
    bins: numpy.array,
    var: str,
    calib_factor: float = 1,
    yscale: float = 1,
    area_normalized: bool = True,
    band: bool = False,
    clip: bool = False,
):
  
  # get data and category information
  data = [df[cat.mask(df)][var] * calib_factor for cat in categories]
  if clip:
      data = [numpy.clip(df[cat.mask(df)][var] * calib_factor, bins[0], bins[-1]) for cat in categories]
  labels     = [cat.label     for cat in categories]
  colors     = [cat.color     for cat in categories]
  hatches    = [cat.hatch     for cat in categories]
  edgecolors = [cat.edgecolor for cat in categories]

  x = 0.5 * (bins[:-1] + bins[1:])
  w = numpy.diff(bins)

  # weighted total for normalization and error band
  total_counts = numpy.zeros(len(bins) - 1)
  for cat, d in zip(categories, data):
      c, _ = numpy.histogram(d, bins=bins)
      total_counts += c * cat.scale
  Ntot = total_counts.sum()

  errors = numpy.sqrt(total_counts)
  if area_normalized:
      scale        = numpy.where(Ntot * w > 0, Ntot * w, 1.0)
      y            = total_counts / scale
      yerr         = errors / scale
  else:
      y    = total_counts
      yerr = errors

  # per-event weights
  weights = []
  for cat, d in zip(categories, data):

    # for area-normalized: cat.scale * yscale / (Ntot * bin_width)
    if area_normalized:
      bin_idx = numpy.clip(
        numpy.searchsorted(bins[:-1], numpy.clip(d, bins[0], bins[-1] - 1e-10), side='right') - 1,
        0, len(bins) - 2
      )
      norm = numpy.where(
        Ntot > 0,
        cat.scale * yscale / (Ntot * w[bin_idx]),
        0
      )
      weights.append(norm)

    # for POT-normalized
    else:
      weights.append(numpy.full(len(d), cat.scale * yscale))

  # plot MC with stacked categories
  ax.hist(
    data,
    stacked   = True,
    histtype  = 'stepfilled',
    bins      = bins,
    label     = labels,
    color     = colors,
    weights   = weights,
    hatch     = hatches,
    edgecolor = edgecolors,
  )

  # error band, if you'd like
  if band:
    ax.bar(
      x,
      2 * yerr * yscale,
      width     = w,
      bottom    = (y - yerr) * yscale,
      fill      = True,
      linewidth = 0,
      color     = 'gray',
      alpha     = 0.5,
    )
  return ax

def plot_by_category_with_offbeam(
    ax,
    df: pandas.DataFrame,
    categories: list[Category],
    bins: numpy.array,
    var: str,
    df_offbeam: pandas.DataFrame,
    offbeam_scale: float = 1,
    calib_factor: float = 1,
    yscale: float = 1,
    area_normalized: bool = True,
    band: bool = False,
    clip: bool = False,
    **kwargs,
):
  
  # define the off-beam category
  offbeam_cat = Category(
    mask      = lambda df: pandas.Series([True] * len(df), index=df.index),
    label     = 'off-beam',
    color     = 'white',
    hatch     = '//',
    edgecolor = 'dodgerblue',
    scale     = offbeam_scale,
  )

  # append the off-beam category, and plot
  plot_by_category(
    ax,
    pandas.concat([df, df_offbeam], ignore_index=True),
    categories + [offbeam_cat],
    bins,
    var,
    calib_factor    = calib_factor,
    yscale          = yscale,
    area_normalized = area_normalized,
    band            = band,
    clip            = clip,
    **kwargs,
  )
  return ax

def plot_data(
  ax,
  df: pandas.DataFrame,
  bins: numpy.array,
  var: str,
  area_normalized: bool = True,
  clip: bool = False,
  **kwargs,
):
  counts, _ = numpy.histogram(df[var], bins=bins)
  if clip:
    counts, _ = numpy.histogram(
      numpy.clip(df[var], bins[0], bins[1]), 
      bins = bins
    )
  errors = numpy.sqrt(counts)

  x = 0.5 * (bins[:-1] + bins[1:])
  widths = numpy.diff(bins)

  if area_normalized:
    Ntot = counts.sum()

    scale = Ntot * widths
    scale[scale == 0] = 1.0

    y = counts / scale
    yerr = errors / scale
  else:
    y = counts
    yerr = errors

  ax.errorbar(
    x,
    y,
    yerr = yerr,
    marker = '.',
    ls = '',
    color = 'black',
    label = 'data',
    **kwargs,
  )

  return ax

def place_cut(
  ax,
  CUT,
  lower = True,
):
  ax.axvline(
    CUT, 
    lw=0.75, c='gray', 
    zorder=-2
  )
  if lower:
    ax.axvspan(
      CUT, ax.get_xlim()[1],  
      fc='None', hatch='//', ec='lightgray',
      zorder=-3
    )
  else: 
    ax.axvspan(
      ax.get_xlim()[0], CUT,
      fc='None', hatch='//', ec='lightgray',
      zorder=-3
    )

  return ax