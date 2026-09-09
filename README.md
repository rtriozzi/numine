# `numine`

This repository contains algorithmic and ML code to select (the `mine` part) neutrino interactions (the `nu` part) in the ICARUS liquid argon time projection chamber (LAr-TPC) imaging detector at Fermilab, using neutrinos from the BNB and NuMI beams (the `numi` part).

The repository covers multiple selection channels (1eNp0π, 1μNp0π, NC, Michels) and provides two complementary approaches: a CAFAna-based cut selection framework (see `cafana`) and machine learning pipelines (see `ml`).
Currently, the CAFAna-based analysis is supported from the simulation of the events, through the selection and creation of trees with systematic uncertainties, and oscillation studies.

A separate Python package, `sbruceana`, is provided to process analysis-level SBN trees with more tools for selections and evaluations, and to achieve an accurate data and simulation normalization.

### Repository structure

```
numine/
├── cafana/       # CAFAna-based selections and tree makers
    └── cc1e0pi/, cc1mu0pi/, ncpi0/, ...    # per-channel selections and tree makers
└── fcl/          # all icaruscode configurations to produce simulated NuMI events and detector model variations
└── fitting/      # all PROfit configurations to carry out oscillation studies
└── ml/           # ML-based selections (e.g., XGBoost BDTs)
    ├── src/      # shared Python modules (preprocessing, training, evaluation, plotting)
    ├── models/   # base model definitions, if needed
    └── cc1e0pi/, ...    # per-channel notebooks and evaluations
└── sbruceana/    # Python analyzers for analysis-level sbruce trees (systematics, selections, normalizations)
    ├── sbruceana/       # actual code...
    └── cc1e0pi/, ncpi0/, ...    # per-channel notebooks
```

## Setup

This analysis heavily relies on the NuGraph2 graph neural network for selections. This is now largely integrated throughout the ICARUS simulation and reconstruction packages. 
Here, setup instructions are provided based on code development branches: note that once the integration with released software is finished, one just would need to setup a recent-enough tag.

Still, here you go.

### ICARUS/SBN software

Our work is based off of icaruscode `v10_06_00_01p01` (and the corresponding LArSoft base), so grab that:

```
# you'll need the SL7 container - use the one that 
# supports LArBatch for submitting jobs to the grid
sh /exp/$(id -ng)/data/users/vito/podman/start_SL7dev_jsl.sh

# setup the icaruscode version
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup icaruscode v10_06_00_01p01 -q e26:prof
```

Go to your app area, setup a development directory, and:
```
mrb newDev
source localProducts*/setup
```

Grab the code you need in srcs:
* icaruscode (where the magic happens):
```
mrb g -t v10_06_00_01p01 icaruscode
git checkout feature/rtriozzi_cerati_NuGraph2_Filter
```
* sbnanaobj (definitions for new CAF variables):
```
mrb g -t $SBNANAOBJ_VERSION sbnanaobj
git checkout feature/rtriozzi_NuGraph2_MoveVars
```
* sbncode (code for actually filling the CAF variables):
```
mrb g -t $SBNCODE_VERSION sbncode
git checkout feature/llena_CAFMaker_PandoraAfterNuGraph
```
* lardataobj (some class definition updates, make sure to grab this directly from Leonardo's branch - I didn't port this to a branch because there were no modifications needed + it was eventually merged):
```
mrb g -t $LARDATAOBJ_VERSION lardataobj
# make sure to point to Leonardo's repo here...
# git@github.com:leonardo-lena/lardataobj.git
git checkout feature-nugraph-multislice
```
* larrecodnn (some base NuGraph2 code - same situation as for lardataobj):
```
mrb g -t $LARRECODNN_VERSION larrecodnn
# make sure to point to Leonardo's repo here...
# git@github.com:leonardo-lena/larrecodnn.git
git checkout feature/nugraph-multislice
```
* larpandora (this is the interface between LArSoft, i.e. our icaruscode magic, and Pandora):
```
mrb g -t $LARPANDORA_VERSION larpandora
git checkout feature/rtriozzi_NuGraphInterface
```
* larpandoracontent (this is purely Pandora code):
```
mrb g -t $LARPANDORACONTENT_VERSION larpandoracontent
git checkout feature/rtriozzi_NuGraphInterface
```

With this, you can produce any event you want with all the NuGraph2 information and some of the recent Pandora developments.
Note that [Pandora](https://github.com/PandoraPFAOrg) code is not strictly needed, unless you want to do reconstruction work.
You could use also my local code (here provided on CNAF machines):
```
/storage/gpfs_data/icarus/local/users/rtriozzi/dev/dev_v10_06_00_01p01_rtriozzi_cerati_NuGraph
```

### CAFAna/SBNAna

The CAFAna-based code requires a SL7 environment, available at FNAL or CNAF.
Set up `sbnana`:
```
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana v10_01_00 -q e26:prof
```

Note that if you want to rely on the novel NuGraph2 capabilities for selections, you might need to rely on the new development branches, so that CAFAna is aware of the new variables in the CAF files. 
From the [ICARUS/SBN software](#icarus-SBN-software), just grab sbnana and sbnanaobj.
Otherwise, you could use also my local code (here provided on CNAF machines):
```
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana v10_01_00 -q e26:prof
setup fife_utils
source /storage/gpfs_data/icarus/local/users/rtriozzi/sbnana/sbnana_NG2/localProducts_larsoft_v10_06_00_e26_prof/setup
```

The CAFAna modules provide tools to select neutrinos, plot distributions for relevant variables (organized by true interaction type or final-state topology), and produce selection efficiency and purity curves at different cut stages, accounting for pile-up. 
CAFAna tree makers are also provided to flatten the CAF structure into ROOT `TTree`s, correctly handling systematic uncertainties.

Run any CAFAna macro with:
```
cafe -bq <CAFAna_Macro.C>
```

making sure input files (defined via `TargetFile`) are accessible via path, wildcard, or `sam`.

## Fitting

Everything about fitting is provided in the [PROfit repository](https://github.com/markrosslonergan/PROfit).
The `fitting/config` folder provides every configuration you could need via PROfit to perform the electron neutrino disappearance (and more) analysis, both with the simple disappearance-only model, and with 3+1.
The `fitting/datamc` provides an XML to perform data/MC comparisons with all the systematic uncertainties.
Everything was done with PROfit 2.4.0, but do expect compatibility with newer tags (modulo some minor XML tweaks).

Note that all the tree makers in `cafana` provide trees with all the systematic uncertainties you need.
A couple of knobs (z-expansion, CC-MEC x-sec shape) are not available in standard CAFs: you can ignore those for a first-pass analysis.
If you want to add those knobs, prepare for pain: you'll need a local GENIE distribution and [sbnnusyst](https://github.com/jedori0228/sbnnusyst). 

## ML

The ML pipeline is Python-based and self-contained. Set up the environment and
install the `src` package:
```
pip install -e ml/
```

Per-channel BDT trainings and evaluations live in the channel subdirectories as
Jupyter notebooks. Shared utilities (feature definitions, preprocessing, training,
evaluation, plotting) are in `ml/src/`.
