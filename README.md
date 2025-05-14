This repository contains code written within the CAFAna framework for selecting neutrinos in the ICARUS liquid argon time projection chamber detector.

You should be able to activate a SL7 environment, either at FNAL or CNAF.
Then, set up `sbnana`:
```
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana v09_93_01 -q e26:prof
setup larsoft_data v1_02_02
```

The repository provides tools to select neutrinos, plot distributions for relevant variables (organized by true interaction tyoe), produce well-defined selection efficiency and purity curves.
Just run the tool you need, making sure that input files (defined through `TargetFile`) are there:
```
cafe -bq <CAFAna_Macro.C>
```
