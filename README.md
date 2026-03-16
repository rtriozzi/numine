This repository contains code written within the CAFAna framework for selecting neutrinos in the ICARUS liquid argon time projection chamber detector.

You should be able to activate a SL7 environment, either at FNAL or CNAF.
Then, set up `sbnana`:
```
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana v10_01_00 -q e26:prof
```

The repository provides tools to select neutrinos, plot distributions for relevant variables (organized by true interaction type, or true final-state topology), produce well-defined selection efficiency and purity curves at different selection steps, accounting for pile-up.
It also provides some CAFAna tree makers, which can be used to flatten the CAF structure into ROOT `TTree`s, also correctly handling systematic uncertainties.
Just run the tool you need, making sure that input files (defined through `TargetFile`) are there (you can use a path, wildcards, or `sam`):
```
cafe -bq <CAFAna_Macro.C>
```
