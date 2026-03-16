# src/evaluation.py

import numpy
import pandas
import xgboost

import matplotlib.pyplot as plt

from src.plotting import ICARUSPlotter
from src.plotting import Category

def make_test_df(
    df     : pandas.DataFrame,
    X_test : pandas.DataFrame,
    model  : xgboost.XGBClassifier,
) -> pandas.DataFrame:

    # go back to the whole dataset
    # from the test indices
    df_test = df.loc[X_test.index].copy()

    # compute predictions and probabilities
    df_test['pred'] = model.predict(X_test)
    df_test['prob'] = model.predict_proba(X_test)[:,1]

    return df_test
    
def plot_precision_recall_threshold(
    thrs      : numpy.ndarray,
    precs     : numpy.ndarray,
    recs      : numpy.ndarray,
    threshold : float = None,
    title     = ICARUSPlotter.DEFAULT_TITLE,
    figsize   = ICARUSPlotter.DEFAULT_FIGSIZE,
) -> tuple[plt.Figure, plt.Axes]:

    # compute F1 score as the harmonic mean
    # of recall and precisions
    f1 = 2 * recs[:-1] * precs[:-1] / (recs[:-1] + precs[:-1])

    fig, ax = ICARUSPlotter.make_fig(figsize)

    ax.plot(thrs, recs[:-1], lw=2, c='C0', label='Recall')
    ax.plot(thrs, precs[:-1], lw=2, c='C1', label='Precision')
    ax.plot(thrs, f1,         lw=2, c='C2', label='F1')

    if threshold is not None:
        ax.axvline(threshold, lw=2, c='red', zorder=-3)
        rec_at_thr  = numpy.interp(threshold, thrs, recs[:-1])
        prec_at_thr = numpy.interp(threshold, thrs, precs[:-1])
        print(f"Recall at {threshold}:    {rec_at_thr:.4f}")
        print(f"Precision at {threshold}: {prec_at_thr:.4f}")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    ICARUSPlotter.apply_style(ax, title=title, xlabel="Threshold", ylabel="Score")
    ICARUSPlotter.apply_legend(ax, loc='lower right', title='Nom. flux')

    return fig, ax

def plot_score_distribution(
    df         : pandas.DataFrame,
    categories : list[Category],
    var        : str = 'prob',
    binning    : list = None,
    xlabel     : str = "BDT score",
    title      : str = ICARUSPlotter.DEFAULT_TITLE,
    figsize    : tuple = ICARUSPlotter.DEFAULT_FIGSIZE,
) -> tuple[plt.Figure, plt.Axes]:

    # default binning for default score variable...
    if binning is None:
        binning = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,
                   0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    data   = [df[cat.mask(df)][var] for cat in categories]
    labels = [cat.label for cat in categories]
    colors = [cat.color for cat in categories]

    fig, ax = ICARUSPlotter.make_fig(figsize)

    ax.hist(data, stacked=True, bins=binning, label=labels, color=colors)
    ax.hist(df[var], bins=binning, histtype='step', ec='black', linewidth=0.75)

    ax.set_xlim(0, 1)
    ax.set_yscale('log')

    ICARUSPlotter.apply_style(ax, title=title, xlabel=xlabel, ylabel="Slices [#]")
    leg = ICARUSPlotter.apply_legend(ax, fontsize=9, loc='upper right', title='Nom. flux')
    leg.get_title().set_fontsize(11)

    return fig, ax