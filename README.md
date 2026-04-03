### `numine`

This repository contains algorithmic and ML code to select (the `mine` part) neutrino interactions (the `nu` part)
in the ICARUS liquid argon time projection chamber (LAr-TPC) imaging detector at Fermilab,
using neutrinos from the BNB and NuMI beams (the `numi` part).

The repository covers multiple selection channels (1eNp0π, 1μNp0π, NC, Michels)
and provides two complementary approaches: a CAFAna-based cut selection framework
(see `cafana`) and machine learning pipelines (see `ml`).

A separate Python package, `sbruceana`, is provided to process analysis-level SBN sbruce trees,
including a first-pass treatment of systematics uncertainties based on covariance matrices, 
more tools for selections and evaluations, and to achieve an accurate data and simulation normalization.

#### Repository structure

```
numine/
├── cafana/       # CAFAna-based selections and tree makers
    └── cc1e0pi/, cc1mu0pi/, ncpi0/, ...    # per-channel selections and tree makers
└── ml/           # ML-based selections (e.g., XGBoost BDTs)
    ├── src/      # shared Python modules (preprocessing, training, evaluation, plotting)
    ├── models/   # base model definitions, if needed
    └── cc1e0pi/, ...    # per-channel notebooks and evaluations
└── sbruceana/    # Python analyzers for analysis-level sbruce trees (systematics, selections, normalizations)
    ├── sbruceana/       # actual code...
    └── cc1e0pi/, ncpi0/, ...    # per-channel notebooks
```

#### CAFAna setup

The CAFAna-based code requires a SL7 environment, available at FNAL or CNAF.
Set up `sbnana`:
```
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana v10_01_00 -q e26:prof
```

The CAFAna modules provide tools to select neutrinos, plot distributions for
relevant variables (organized by true interaction type or final-state topology),
and produce selection efficiency and purity curves at different cut stages,
accounting for pile-up. CAFAna tree makers are also provided to flatten the CAF
structure into ROOT `TTree`s, correctly handling systematic uncertainties.

Run any CAFAna macro with:
```
cafe -bq <CAFAna_Macro.C>
```

making sure input files (defined via `TargetFile`) are accessible via path, wildcard, or `sam`.

#### ML setup

The ML pipeline is Python-based and self-contained. Set up the environment and
install the `src` package:
```
pip install -e ml/
```

Per-channel BDT trainings and evaluations live in the channel subdirectories as
Jupyter notebooks. Shared utilities (feature definitions, preprocessing, training,
evaluation, plotting) are in `ml/src/`.
