# src/config.py

from dataclasses import dataclass

'''
  Configuration to handle plotting nicely.
'''

@dataclass
class Category:
  label     : str           # label in plots
  color     : str           # color of the histograms
  mask      : callable      # mask to define the category, e.g., in truth
  scale     : float = 1.0   # how much to scale the category (e.g., for data-driven cosmics from off-beam)
  hatch     : str   = None  # whether to hatch the histogram
  edgecolor : str   = None  # edgecolor of the step-like histograms

CC1E0PI_CATEGORIES = [
  Category(
    label = '1eNp0π',
    color = 'orangered',
    mask  = lambda df: (df.signal == 1),
  ),
  Category(
    label = '$\\nu_e$CC',
    color = 'darkorange',
    mask  = lambda df: (df.signal != 1) & (df.nue == 1) & (df.CC == 1) & (df.trueFV == 1),
  ),
  Category(
    label = '$\\nu_{\\mu}$CC$\\pi^0$',
    color = 'violet',
    mask  = lambda df: (df.signal != 1) & (df.numu == 1) & (df.CC == 1) & (df.trueFV == 1) & (df.ispi0 == 1),
  ),
  Category(
    label = '$\\nu_{\\mu}$CC',
    color = 'darkorchid',
    mask  = lambda df: (df.signal != 1) & (df.numu == 1) & (df.CC == 1) & (df.trueFV == 1) & (df.ispi0 == 0),
  ),
  Category(
    label = '$\\nu$NC',
    color = 'pink',
    mask  = lambda df: (df.signal != 1) & (df.CC == 0) & (df.trueFV == 1),
  ),
  Category(
    label = 'OoFV',
    color = 'cyan',
    mask  = lambda df: (df.signal != 1) & (df.trueOOFV == 1),
  ),
  Category(
    label = 'cosmics',
    color = 'dodgerblue',
    mask  = lambda df: (df.cosmic == 1),
  ),
]

NCPI0_CATEGORIES = [
  Category(
    label = 'NC$\\pi^0$',
    color = 'orangered',
    mask  = lambda df: (df.signal == 1),
  ),
  Category(
    label = '$\\nu_e$CC',
    color = 'darkorange',
    mask  = lambda df: (df.signal != 1) & (df.nue == 1) & (df.CC == 1) & (df.trueFV == 1),
  ),
  Category(
    label = '$\\nu_{\\mu}$CC$\\pi^0$',
    color = 'violet',
    mask  = lambda df: (df.signal != 1) & (df.numu == 1) & (df.CC == 1) & (df.trueFV == 1) & (df.ispi0 == 1),
  ),
  Category(
    label = '$\\nu_{\\mu}$CC',
    color = 'darkorchid',
    mask  = lambda df: (df.signal != 1) & (df.numu == 1) & (df.CC == 1) & (df.trueFV == 1) & (df.ispi0 == 0),
  ),
  Category(
    label = '$\\nu$NC',
    color = 'pink',
    mask  = lambda df: (df.signal != 1) & (df.CC == 0) & (df.trueFV == 1),
  ),
  Category(
    label = 'OoFV',
    color = 'cyan',
    mask  = lambda df: (df.signal != 1) & (df.trueOOFV == 1),
  ),
  Category(
    label = 'cosmics',
    color = 'dodgerblue',
    mask  = lambda df: (df.cosmic == 1),
  ),
]

'''
  Configuration for handling systematics.
'''
SYSTS_CONFIG = {
  'spline': [
    # x-sec
    {'name': 'MaCCRES', 'tag': 'OtherXSec', 'branch': 'GENIEReWeight_SBN_v1_multisigma_MaCCRES'},
    {'name': 'MvCCRES', 'tag': 'OtherXSec', 'branch': 'GENIEReWeight_SBN_v1_multisigma_MvCCRES'},
    {'name': 'RPA_CCQE', 'tag': 'QE-MEC', 'branch': 'GENIEReWeight_SBN_v1_multisigma_RPA_CCQE'},
    {'name': 'CoulombCCQE', 'tag': 'QE-MEC', 'branch': 'GENIEReWeight_SBN_v1_multisigma_CoulombCCQE'},
    {'name': 'NormCCMEC', 'tag': 'QE-MEC', 'branch': 'GENIEReWeight_SBN_v1_multisigma_NormCCMEC'},
    {'name': 'FrCEx_N', 'tag': 'FSI', 'branch': 'GENIEReWeight_SBN_v1_multisigma_FrCEx_N'},
    {'name': 'FrCEx_pi', 'tag': 'FSI', 'branch': 'GENIEReWeight_SBN_v1_multisigma_FrCEx_pi'},
    {'name': 'DecayAngMEC', 'tag': 'QE-MEC', 'branch': 'GENIEReWeight_SBN_v1_multisigma_DecayAngMEC'},
    # flux
    *[{'name': f'numi_pc_{i}', 'tag': 'fluxsyst'} for i in range(8)],
    {'name': 'numi_beam_div', 'tag': 'fluxsyst'},
    {'name': 'numi_beam_power', 'tag': 'fluxsyst'},
    {'name': 'numi_beam_shift_y_minus', 'tag': 'fluxsyst'},
    {'name': 'numi_beam_shift_y_plus', 'tag': 'fluxsyst'},
    {'name': 'numi_beam_shift_x', 'tag': 'fluxsyst'},
    {'name': 'numi_beam_spot', 'tag': 'fluxsyst'},
    {'name': 'numi_horn1_x', 'tag': 'fluxsyst'},
    {'name': 'numi_horn1_y', 'tag': 'fluxsyst'},
    {'name': 'numi_horn_current_plus', 'tag': 'fluxsyst'},
    {'name': 'numi_water_layer', 'tag': 'fluxsyst'},
  ],
  'multisim': [
    {'name': 'NCELVariationResponse', 'tag': 'OtherXSec', 'branch': 'GENIEReWeight_SBN_v1_multisim_NCELVariationResponse'},
    {'name': 'NCRESVariationResponse', 'tag': 'OtherXSec', 'branch': 'GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse'},
    {'name': 'COHVariationResponse',  'tag': 'OtherXSec', 'branch': 'GENIEReWeight_SBN_v1_multisim_COHVariationResponse'},
    {'name': 'DISBYVariationResponse', 'tag': 'OtherXSec', 'branch': 'GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse'},
    {'name': 'FSI_pi', 'tag': 'FSI', 'branch': 'GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse'},
    {'name': 'FSI_N', 'tag': 'FSI', 'branch': 'GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse'},
    {'name': 'NormNCMEC', 'tag': 'QE-MEC', 'branch': 'GENIEReWeight_SBN_v1_multisim_NormNCMEC'},
    *[{'name': f'NonRESBGvpCC{i}pi', 'tag': 'NonRESBG', 'branch': f'GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC{i}pi'} for i in [1,2]],
    *[{'name': f'NonRESBGvpNC{i}pi', 'tag': 'NonRESBG', 'branch': f'GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC{i}pi'} for i in [1,2]],
    *[{'name': f'NonRESBGvnCC{i}pi', 'tag': 'NonRESBG', 'branch': f'GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC{i}pi'} for i in [1,2]],
    *[{'name': f'NonRESBGvnNC{i}pi', 'tag': 'NonRESBG', 'branch': f'GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC{i}pi'} for i in [1,2]],
    *[{'name': f'NonRESBGvbarp{cc}{i}pi', 'tag': 'NonRESBG', 'branch': f'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarp{cc}{i}pi'} for cc in ['CC','NC'] for i in [1,2]],
    *[{'name': f'NonRESBGvbarn{cc}{i}pi', 'tag': 'NonRESBG', 'branch': f'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarn{cc}{i}pi'} for cc in ['CC','NC'] for i in [1,2]],
  ],
  'norm': [
    {'name': 'FluxNorm', 'tag': 'other', 'value': 0.02},
    {'name': 'FiducialVol', 'tag': 'other', 'value': 0.01},
    # {'name': 'DummyDetSyst','tag': 'other', 'value': 0.10},
  ],
}