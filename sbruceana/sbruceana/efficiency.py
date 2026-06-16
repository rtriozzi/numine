# sbruceana/efficiency.py

import numpy
import matplotlib

from sbruceana.utils import clopper_pearson

def plot_efficiency(
    ax,
    df,
    var         = 'trueE',
    binning     = numpy.arange(0, 4+0.2, 0.2),
    cuts_idx    = None,    # list of cut indices to overlay, e.g. [0,2,5]; `None` is all
    cut_labels  = None,    # list of cut labels; `None` is defaulted
    alpha_cp    = 0.682,   # Clopper-Pearson coverage
):
    # handle cuts gracefully
    if cuts_idx is None:
        cuts_idx = list(range(9))

    if cut_labels is None:
        cut_labels = [
            'Presel.', 'FM', 'Electron ID', 'dE/dx',
            'Angle', 'Gap', 'Proton ID', '0π', 'LE-μ veto',
        ]
    
    assert len(cuts_idx) == len(cut_labels), \
        'You messed up either the cuts, or the cut labels.'
 
    bin_centres = 0.5 * (binning[:-1] + binning[1:])
 
    # denominator of the efficiency, scale it to [0, 1]
    # so that it's a nice background
    n_total, _ = numpy.histogram(df[var], bins=binning)
    ax_twin = ax.twinx()
    ax_twin.hist(df[var], bins=binning, histtype='stepfilled',
                 color='lightgray', alpha=0.5, zorder=-3)
    ax_twin.set_ylabel('counts [#]', color='gray')
    ax_twin.tick_params(axis='y', labelcolor='gray')
    ax_twin.set_ylim(bottom=0)
    fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
    fmt.set_scientific(True)
    fmt.set_powerlimits((0, 0))
    ax_twin.yaxis.set_major_formatter(fmt)

    # efficiency curves
    counter = 0
    for i in cuts_idx:
        cut_name = f'cut{i}'
        n_pass, _ = numpy.histogram(df.loc[df[cut_name] == 1, var], bins=binning)
        eff, err_lo, err_hi = clopper_pearson(n_pass, n_total, alpha=alpha_cp)
 
        overall_eff = (df[cut_name] == 1).sum() / max(len(df), 1)
        label = f'{cut_labels[counter]} ({100*overall_eff:.1f}%)'
        
        ax.errorbar(
            bin_centres, eff,
            xerr = 0.5 * numpy.diff(binning),
            yerr        = [err_lo, err_hi],
            fmt         = '',
            ls          = '',
            # markersize  = 3,
            linewidth   = 1.25,
            label       = label,
            zorder      = 1,
        )

        counter += 1
 
    # gfx
    ax.axhline(1, color='black', linewidth=0.75, linestyle='-', zorder=0)
    ax.set(
        xlim = (binning[0], binning[-1]),
        ylim = (0, 1.05),
        xlabel = var,
        ylabel = 'efficiency',

    )

    # make the count histogram visible
    ax.set_zorder(ax_twin.get_zorder() + 1)
    ax.patch.set_visible(False)

    # legend
    # ax.legend(fontsize=7, ncol=2, loc='upper right',
    #           framealpha=0.9, edgecolor='none')
 
    return ax