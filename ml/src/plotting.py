# src/plotting.py

import matplotlib.pyplot as plt

from dataclasses import dataclass

"""
    This class provides basic utilities and styling 
    for matplotlib plots, with a ICARUS flavor.
"""
class ICARUSPlotter:

    DEFAULT_TITLE   = "ICARUS NuMI v9 MC"
    DEFAULT_FIGSIZE = (4, 3)

    @staticmethod
    def make_fig(figsize=DEFAULT_FIGSIZE):
        return plt.subplots(figsize=figsize, ncols=1, layout='constrained')

    @staticmethod
    def apply_style(ax, title=DEFAULT_TITLE, xlabel="", ylabel="", sci_y=False):
        ax.set_title(title, fontsize=16, loc='right', color='gray')
        ax.set_xlabel(xlabel, fontsize=16, loc='right')
        ax.set_ylabel(ylabel, fontsize=16, loc='top')

        ax.tick_params(labelsize=14, length=6, direction='in', right=True, top=True)
        ax.tick_params(which='minor', length=3, direction='in', right=True, top=True)
        ax.minorticks_on()

        if sci_y:
            ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0), useMathText=True)
            ax.yaxis.offsetText.set_fontsize(13)

    @staticmethod
    def apply_legend(ax, **kwargs):
        defaults = dict(frameon=True, fancybox=False, fontsize=10,
                        ncol=1, loc='lower right', handlelength=1, columnspacing=1)
        return ax.legend(**{**defaults, **kwargs})
        

"""
    This dataclass can be used to defined true 
    categories in plots broken by truth definitions.
"""
@dataclass
class Category:
    label : str
    color : str
    mask  : callable