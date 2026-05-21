# src/plotting.py

import numpy
import pandas

import matplotlib
import matplotlib.pyplot as plt

from .config import Category

import dataclasses

def plot_var(
    ax,
    df: pandas.DataFrame,
    bins: numpy.array,
    var: str,
    weight: float = 1.0,
    **kwargs
):
  ax.hist(
    df[var],
    bins = bins,
    weights = numpy.full(len(df[var]), weight),
    **kwargs
  )

  return ax

def plot_var_by_category(
    ax,
    df: pandas.DataFrame,
    bins: numpy.array,
    var: str,
    flavor: str = 'nue',
    band: bool = True,
    **kwargs,
):
  match flavor:
    case 'nue':
      mask = (df.CC == 1) & (abs(df.truePDG) == 12)
      label = ['$\\nu_e$CC', 'other']
    case 'numu':
      mask = (df.CC == 1) & (abs(df.truePDG) == 14)
      label = ['$\\nu_{\\mu}$CC', 'other'] 
    case _:
      print('Flavor tag not supported, yet.')

  ax.hist(
    [df[mask][var], df[~mask][var]],
    stacked = True,
    histtype = 'stepfilled',
    bins    = bins,
    label   = label
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
        facecolor = 'None',
        hatch     = '\\\\\\',
        edgecolor = 'gray',
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
  raw_counts   = numpy.zeros(len(bins) - 1)
  total_counts = numpy.zeros(len(bins) - 1)
  for cat, d in zip(categories, data):
    c, _ = numpy.histogram(d, bins=bins)
    raw_counts   += c
    total_counts += c * cat.scale
  Ntot = total_counts.sum()

  if area_normalized:
    scale = numpy.where(Ntot * w > 0, Ntot * w, 1.0)
    y    = total_counts / scale
    yerr = numpy.sqrt(raw_counts) / scale
  else:
    y    = total_counts
    safe = numpy.where(raw_counts > 0, raw_counts, 1.0)
    yerr = numpy.sqrt(raw_counts) * total_counts / safe

  # per-event weights
  weights = []
  for cat, d in zip(categories, data):
    if area_normalized:
      bin_idx = numpy.clip(
        numpy.searchsorted(bins[:-1], numpy.clip(d, bins[0], bins[-1] - 1e-10), side='right') - 1,
        0, len(bins) - 2
      )
      weights.append(numpy.where(Ntot > 0, cat.scale * yscale / (Ntot * w[bin_idx]), 0))
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
        2 * yerr,
        width     = w,
        bottom    = y - yerr,
        fill      = True,
        linewidth = 0,
        facecolor = 'None',
        hatch     = '\\\\\\\\',
        edgecolor = 'gray',
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
):
    # get data and category information
    data = [df[cat.mask(df)][var].to_numpy() for cat in categories]
    if clip:
        data = [numpy.clip(d, bins[0], bins[-1]) for d in data]
    labels     = [cat.label     for cat in categories]
    colors     = [cat.color     for cat in categories]
    hatches    = [cat.hatch     for cat in categories]
    edgecolors = [cat.edgecolor for cat in categories]
    scales     = [cat.scale * yscale for cat in categories]

    # offbeam
    off_data = df_offbeam[var].to_numpy() * calib_factor
    if clip:
        off_data = numpy.clip(off_data, bins[0], bins[-1])
    data.append(off_data)
    labels.append('off-beam')
    colors.append('white')
    hatches.append('//')
    edgecolors.append('dodgerblue')
    scales.append(offbeam_scale)

    # print scaled counts per category and total
    max_label_len = max(len(l) for l in labels)
    print(f"{'Category':<{max_label_len}}  {'Raw':>10}  {'Scaled':>12}  {'%':>7}")
    print("-" * (max_label_len + 36))
    scaled_counts = [len(d) * s for d, s in zip(data, scales)]
    total_scaled  = sum(scaled_counts)
    pcts = [100.0 * sc / total_scaled if total_scaled > 0 else 0.0 for sc in scaled_counts]
    for label, d, scaled, pct in zip(labels, data, scaled_counts, pcts):
        print(f"{label:<{max_label_len}}  {len(d):>10d}  {scaled:>12.2f}  {pct:>6.2f}%")
    print("-" * (max_label_len + 36))
    print(f"{'Total':<{max_label_len}}  {'':>10}  {total_scaled:>12.2f}  {'100.00%':>7}")

    # append percentage to legend labels
    labels = [f"{l} ({p:.1f}%)" for l, p in zip(labels, pcts)]

    x = 0.5 * (bins[:-1] + bins[1:])
    w = numpy.diff(bins)

    # compute total counts and sum of weights^2
    total_counts = numpy.zeros(len(bins) - 1)
    sumw2        = numpy.zeros(len(bins) - 1)
    for s, d in zip(scales, data):
        total_counts += numpy.histogram(d, bins=bins, weights=numpy.full(len(d), s))[0]
        sumw2        += numpy.histogram(d, bins=bins, weights=numpy.full(len(d), s**2))[0]

    # normalization
    if area_normalized:
        Ntot  = total_counts.sum()
        scale = numpy.where(Ntot * w > 0, Ntot * w, 1.0)
        y    = total_counts / scale
        yerr = numpy.sqrt(sumw2) / scale
    else:
        y    = total_counts
        yerr = numpy.sqrt(sumw2)

    # per-event weights for plotting
    weights = []
    for s, d in zip(scales, data):
        if area_normalized:
            bin_idx = numpy.clip(
                numpy.searchsorted(bins[:-1], numpy.clip(d, bins[0], bins[-1] - 1e-10), side='right') - 1,
                0, len(bins) - 2
            )
            weights.append(numpy.where(
                total_counts.sum() > 0,
                s / (total_counts.sum() * w[bin_idx]),
                0
            ))
        else:
            weights.append(numpy.full(len(d), s))

    # stacked histogram, with offbeam
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

    # error band
    if band:
        ax.bar(
            x,
            2 * yerr,
            width     = w,
            bottom    = y - yerr,
            fill      = True,
            linewidth = 0,
            facecolor = 'None',
            hatch     = '\\\\\\\\',
            edgecolor = 'gray',
        )

    return ax

def plot_data(
  ax,
  df: pandas.DataFrame,
  bins: numpy.array,
  var: str,
  area_normalized: bool = True,
  clip: bool = False,
  cutoff: float = None,
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

  if cutoff is not None:
    yerr = yerr[x <= cutoff]
    y    = y[x <= cutoff]
    x    = x[x <= cutoff]

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

def plot_by_category(
    ax,
    df: pandas.DataFrame,
    categories: list[Category],
    bins: numpy.array,
    var: str,
    yscale: float = 1,
    band: bool = False,
):
    # get data and category information
    data       = [df[cat.mask(df)][var].to_numpy() for cat in categories]
    labels     = [cat.label     for cat in categories]
    colors     = [cat.color     for cat in categories]
    hatches    = [cat.hatch     for cat in categories]
    edgecolors = [cat.edgecolor for cat in categories]
    scales     = [cat.scale * yscale for cat in categories]

    # per-event weights
    weights = [numpy.full(len(d), s) for s, d in zip(scales, data)]

    # print scaled counts per category and total
    max_label_len = max(len(l) for l in labels)
    print(f"{'Category':<{max_label_len}}  {'Raw':>10}  {'Scaled':>12}  {'%':>7}")
    print("-" * (max_label_len + 36))
    scaled_counts = [len(d) * s for d, s in zip(data, scales)]
    total_scaled  = sum(scaled_counts)
    pcts = [100.0 * sc / total_scaled if total_scaled > 0 else 0.0 for sc in scaled_counts]
    for label, d, scaled, pct in zip(labels, data, scaled_counts, pcts):
        print(f"{label:<{max_label_len}}  {len(d):>10d}  {scaled:>12.2f}  {pct:>6.2f}%")
    print("-" * (max_label_len + 36))
    print(f"{'Total':<{max_label_len}}  {'':>10}  {total_scaled:>12.2f}  {'100.00%':>7}")

    # append percentage to legend labels
    labels = [f"{l} ({p:.1f}%)" for l, p in zip(labels, pcts)]
    
    # stacked histogram
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

    # optional error band
    if band:
        x = 0.5 * (bins[:-1] + bins[1:])
        w = numpy.diff(bins)

        total_counts = numpy.zeros(len(bins) - 1)
        sumw2        = numpy.zeros(len(bins) - 1)
        for s, d in zip(scales, data):
            total_counts += numpy.histogram(d, bins=bins, weights=numpy.full(len(d), s))[0]
            sumw2        += numpy.histogram(d, bins=bins, weights=numpy.full(len(d), s**2))[0]

        y    = total_counts
        yerr = numpy.sqrt(sumw2)

        ax.bar(
            x,
            2 * yerr,
            width     = w,
            bottom    = y - yerr,
            fill      = True,
            linewidth = 0,
            facecolor = 'None',
            hatch     = '\\\\\\',
            edgecolor = 'gray',
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