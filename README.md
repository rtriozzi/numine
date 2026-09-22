# `numine`

This repository contains algorithmic and ML code to select (the `mine` part) neutrino interactions (the `nu` part) in the ICARUS liquid argon time projection chamber (LAr-TPC) imaging detector at Fermilab, using neutrinos from the BNB and NuMI beams (the `numi` part).

The repository covers multiple selection channels (1eNp0π, 1μNp0π, NC, Michels) and provides two complementary approaches: a CAFAna-based cut selection framework (see `cafana`) and machine learning pipelines (see `ml`).
Currently, the CAFAna-based analysis is supported from the simulation of the events, through the selection and creation of trees with systematic uncertainties, and oscillation studies.

A separate Python package, `sbruceana`, is provided to process analysis-level SBN trees with more tools for selections and evaluations, and to achieve an accurate data and simulation normalization.

### Table of contents

- [Repository structure](#repository-structure)
- [Setup](#setup)
  - [ICARUS/SBN software](#icarussbn-software)
    - [Full setup](#full-setup)
    - [Pre-installed codebases](#pre-installed-codebases)
  - [CAFAna/SBNAna](#cafanasbnana)
  - [TITUS](#titus)
- [Systs & Fits](#systs--fits)
  - [Systematics](#systematics)
    - [Using a distributed GENIE version](#using-a-distributed-genie-version)
    - [Using a local GENIE version](#using-a-local-genie-version)
  - [PROfit](#profit)
- [ML](#ml)

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

You'll need to setup the code, or use my installed codebase and/or the corresponding tarball, if you want to simulate events or process data with the latest reconstruction updates.

#### Full setup

Our work is based off of icaruscode `v10_06_00_01p01` (and the corresponding LArSoft base), so grab that for starters:

```
# you'll need the SL7 container - use the one that 
# supports LArBatch for submitting jobs to the grid, e.g.:
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
Note that [Pandora](https://github.com/PandoraPFAOrg) code is not strictly needed for Monte Carlo, unless you want to do reconstruction work.
Some Pandora updates are indeed needed for calibrating for the finite electron lifetime and other non-uniformities in electromagnetic shower energy.
You can use those branches, in that case:
```
# icaruscode (basically fcl updates and database settings)
git checkout feature/rtriozzi_cerati_NuGraph2_Filter_WithShowerCalibration

# larpandora (norm tools interface)
git checkout rtriozzi_NuGraphInterface_ShowerNormTools
```
Please remember that the Pandora branches come from my own forks, set the remote correspondingly.


#### Pre-installed codebases:

You could use rely on pre-installed code in my areas:
* CNAF:
```
/storage/gpfs_data/icarus/local/users/rtriozzi/dev/dev_v10_06_00_01p01_rtriozzi_cerati_NuGraph
```
* FNAL (here including data calibration code):
```
/exp/icarus/app/users/rtriozzi/dev/dev_v10_06_00_01p01_llena_rtriozzi_cerati_NuGraph
```


### CAFAna/SBNAna

The CAFAna-based code requires a SL7 environment, available at FNAL or CNAF.
At FNAL, for example:
```
# this *had* LArBatch support - ask Vito if this is not the case now
sh /exp/$(id -ng)/data/users/vito/podman/start_SL7dev_jsl.sh

# while we're at it, get your token
htgettoken -a htvaultprod.fnal.gov -i icarus
```

Set up `sbnana`:
```
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana v10_01_00 -q e26:prof
```

Note that if you want to rely on the novel NuGraph2 capabilities for selections, you might need to rely on the new development branches, so that CAFAna is aware of the new variables in the already-produced CAF files. 
From the #icarussbn-software setup guide, just grab the sbnana and sbnanaobj branches and you're set: you might need to bump the dependencies in this case (just do so in the `ups/product_deps).
Otherwise, use my pre-installed code:
* FNAL:
```
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana v10_01_00 -q e26:prof
source sbnana_NG2/localProducts_larsoft_v10_06_00_e26_prof/setup

# properly set SAM...
setup fife_utils
export SAM_EXPERIMENT=sbn
export SAM_GROUP=sbn
export SAM_STATION=sbn
export IFDH_BASE_URI=http://samsbn.fnal.gov:8480/sam/sbn/api/
```
* CNAF:
```
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana v10_01_00 -q e26:prof
source /storage/gpfs_data/icarus/local/users/rtriozzi/sbnana/sbnana_NG2/localProducts_larsoft_v10_06_00_e26_prof/setup
```

The CAFAna modules provide tools to select neutrinos, plot distributions for relevant variables (organized by true interaction type or final-state topology), and produce selection efficiency and purity curves at different cut stages, accounting for pile-up. 
CAFAna tree makers are also provided to flatten the CAF structure into ROOT `TTree`s, correctly handling systematic uncertainties.

Run any CAFAna macro with:
```
cafe -bq <CAFAna_Macro.C>
```

making sure input files (defined via `TargetFile`) are accessible via path, wildcard, or `sam`.

### TITUS

The best thing about electron neutrino interactions in LAr-TPCs is that they're beautiful to see, and also pretty easy to spot (by eye, at least).
Please, validate your event selections by looking at event displays, in both Monte Carlo and data.

General information and source code can be found in the [TITUS GitHub page](https://github.com/TITUS-EVD/gallery-framework).
Some instructions for CNAF can be found [here](https://wiki.infn.it/progetti/icarus/display).
Some instructions for the GPVMs can be found [here](https://sbnsoftware.github.io/sbndcode_wiki/TITUS_Event_Display.html).
The following works at FNAL.

Set up the VNC on the GPVM:
```
# check for a free port, e.g. 20
vncserver :20 -localhost -bs

# get the list of running VNCs
vncserver -list

# push output of VNC
export DISPLAY=localhost:20
```
Then, locally:
```
# create the VNC tunnel
ssh -L 5920:localhost:5920 -N -f -l rtriozzi icarusgpvm02.fnal.gov

# if the port is busy on you end, just free it 
lsof -i :5920
kill -9 <PID-of-process-you-want-to-kill>

# open the VNC
open vnc://localhost:5920

# open the terminal, SL7 container, etc. within the VNC session!
```

On the GPVM, you can eventually terminate the VNC:
```
vncserver -kill :20
```
In general, the VNC should not be started from inside the SL7 container, bur from the GPVMs. 
Once a VNC session is open, a terminal can be started from there and a SL7 container opened.
Once there, you can setup `icaruscode` and a python virtual environment.

To use TITUS, Stage0 and Stage1 files will work out of the box, by plotting ChannelROIs (after wire deconvolution).
For Stage1, RawDigits can also be retained from the Stage0 (or the DAQ files) to have a proper EVD of the raw event.
The deconvolution is more than enough in most cases.

On the VNC:
```
### move to the SL7 container!!!
source /exp/$(id -ng)/data/users/vito/podman/start_SL7dev.sh

### move to a Python virtual environment!!!
source env/bin/activate

# just in case, update pip to latest version
pip install --upgrade pip

# install pyqt and pyqtgraph in Python3
pip install PyQt5==5.15.6
pip install pyqtgraph==0.12.4

# updated source file
source /exp/sbnd/app/users/sbnd/static_evd/setup.sh
```

You're done:
```
# run it!
# -i stands for ICARUS (-s for SBND)
# use stage1 files with RawDigits (or Wires)
# stage0 files can be used to display ChannelROIs (after deconvolution!)
evd.py -i <stage(0|1)_file.root>
```

A lifesaver for me was `xclip`, to get high-resolution displays saved:
```
# take a screenshot via `File`

# then in the terminal:
xclip -selection clipboard -t image/png -o > screenshot.png
```

## Systs & Fits

### Systematics

Note that all the tree makers in `cafana` provide trees with all the systematic uncertainties you need.
A couple of knobs (z-expansion, CC-MEC x-sec shape) are not available in standard CAFs: you can ignore those for a first-pass analysis.

If you want to add those knobs yourself, prepare for pain: you'll need a local GENIE distribution and [sbnnusyst](https://github.com/jedori0228/sbnnusyst). An exhaustive (exhausting, perhaps?) guide that worked for me follows.

#### Using a distributed GENIE version

You can use a standard GENIE distribution. At this stage, it will have the z-expansion systematics that are pretty relevant for QE-related systematics (note that if you want to do a QE-like analysis, it's pretty important to get QE systematics right). 
Some MEC knobs are still missing from the main GENIE distributions, though. 

The only thing to pay attention to here is to get `genie` and `genie_xsec` tags that are compatible with the `sbnanaobj` version you need for dealing with your CAFs. 
In general, you can rely on `ups` starting from the `icaruscode` release you need. 
You'll need to setup some more things, but they're less important: if there are conflicts, just unsetup and setup the products you get told about in the logs.

Initial setup:
```
# make sure to be in the SL7 container, and setup the standard ICARUS products

# setup sbnanaobj for dealing with your CAFs
setup sbnanaobj v10_00_04 -q e26:prof

# setup an updated CMake
# we'll need an updated version for NuSystematics
setup cmake v3_27_4

# setup GENIE distributions compatible with this sbnanaobj
setup genie v3_04_02 -qe26:prof
setup genie_xsec v3_04_00 -q AR2320i00000:e1000:k250

# you'll need a recent boost for NuSystematics
setup boost v1_82_0 -qe26:prof
setup eigen v23_08_01_66e8f
setup fhiclcpp v4_18_04 -qe26:prof
```

Now, we'll install NuSystematics:
```
# install NuSystematics
### this is needed only for the first setup!!!
### once you go through this purgatory, you'll just need to setup the already-installed packages
cd ${mywd} # go to your working directory
mkdir nusystematics; cd nusystematics
git clone git@github.com:NuSystematics/nusystematics.git nusystematics-src
mkdir build; cd build
cmake ../nusystematics-src/
make install

# setup NuSystematics
source /exp/icarus/app/users/rtriozzi/nusystematics/build/Linux/bin/setup.fhicl_cpp_standalone.sh
source /exp/icarus/app/users/rtriozzi/nusystematics/build/Linux/bin/setup.systematicstools.sh
source /exp/icarus/app/users/rtriozzi/nusystematics/build/Linux/bin/setup.nusystematics.sh
```

And now, we'll install SBNNuSyst:
```
# install SBNNuSyst
### this is needed only for the first setup!!!
### once you go through this purgatory, you'll just need to setup the already-installed packages
# #{mywd} is your working area
cd ${mywd} # go to your working directory
mkdir sbnnusyst; cd sbnnusyst;
git clone git@github.com:jedori0228/sbnnusyst.git sbnnusyst-src
mkdir build; cd build
cmake ../sbnnusyst-src/
make install

# setup SBNNuSyst
source /exp/icarus/app/users/rtriozzi/sbnnusyst/build/Linux/bin/setup.sbnnusyst.sh
```

We'll need a fcl specifying what systematics we'd like to add. 
You can find those in `numine/fcl/systs`, or, at FNAL, in /exp/icarus/app/users/rtriozzi/All.ParamHeader_NoMEC.fcl.

We'll need to add some template files manually, as they weigh a lot and cannot be distributed through GitHub. 
Either add them to you area or change the corresponding paths in the fcl:
```
# FSI
/exp/icarus/app/users/jskim/sbnnusyst/nusystematics/build/Linux/data/FSI_reweight_template.root

# RPA
/exp/icarus/app/users/PRoy/packages/GENIE/GENIE_3_04_00/nusystematics/nusystematics-src/data/output_RPAReweight.root
```

And then run the weight updater on a list of non-flat CAFs (here called `input_cafs.txt`):`
```
UpdateReweight -c All.ParamHeader_NoMEC.fcl -i input_cafs.txt -o output_flat.caf.root
```

Note that the default GENIE distribution misses some MEC-related knobs.

#### Using a local GENIE version

If you got here, you're already doing great. But, there's more pain now.

Using a local GENIE distribution requires building locally a compatible tag; building with the wanted GENIE Reweight branch, and then re-building NuSystematics with a specific branch:
```
# of course, setup the SL7 container and the ICARUS products

# setup sbnanaobj for dealing with CAFs
setup sbnanaobj v10_00_04 -q e26:prof

# setup an updated CMake 
# we'll need it for NuSystematics
setup cmake v3_27_4

# setup GENIE distributions compatible with your sbnanaobj
setup genie v3_04_02 -qe26:prof
setup genie_xsec   v3_04_00 -q AR2320i00000:e1000:k250

# you'll need a recent boost for NuSystematics
setup boost v1_82_0 -qe26:prof
setup eigen v23_08_01_66e8f
setup fhiclcpp v4_18_04 -qe26:prof

# prepare GENIE directory
cd /exp/icarus/app/users/rtriozzi/
mkdir GENIE; cd GENIE
mkdir GENIE_3_04_02; cd GENIE_3_04_02

##########
# STUFF BELOW CAN GO IN A SETUP SCRIPT AFTER THE 
# FIRST TIME YOU'VE DEALT WITH THIS
# IN GENERAL, THIS STUFF SHOULD BE RAN FROM /exp/icarus/app/users/rtriozzi/GENIE/GENIE_3_04_02

GENIEVERSION=3_04_02

export GENIE_VERSION=v${GENIEVERSION}

## PUT YOUR OWN DIRECTORY HERE
export GENIE_FQ_DIR=`pwd` # /icarus/app/users/jskim/GENIE/GENIE_3_04_00/

## GENIE:
export GENIE=${GENIE_FQ_DIR}/Generator/

## GENIE_DIR: /cvmfs/larsoft.opensciencegrid.org/products/genie/v3_00_06p/
unset GENIE_DIR

## GENIE_INC: /cvmfs/larsoft.opensciencegrid.org/products/genie/v3_00_06p/Linux64bit+3.10-2.17-e20-prof/include
export GENIE_INC=${GENIE_FQ_DIR}/include/

## GENIE_LIB: /cvmfs/larsoft.opensciencegrid.org/products/genie/v3_00_06p/Linux64bit+3.10-2.17-e20-prof/lib
export GENIE_LIB=${GENIE_FQ_DIR}/lib/

## GENIE_REWEIGHT: /cvmfs/larsoft.opensciencegrid.org/products/genie/v3_00_06p/Linux64bit+3.10-2.17-e20-prof/GENIE-Reweight
export GENIE_REWEIGHT=${GENIE_FQ_DIR}/Reweight/

export LD_LIBRARY_PATH=${GENIE_FQ_DIR}/lib/:${LD_LIBRARY_PATH}
export PATH=${GENIE_FQ_DIR}/bin/:${PATH}
export ROOT_INCLUDE_PATH=${GENIE_FQ_DIR}/include/GENIE/:${ROOT_INCLUDE_PATH}
export CMAKE_PREFIX_PATH=${GENIE_FQ_DIR}:${CMAKE_PREFIX_PATH}

# STUFF ABOVE CAN GO IN A SETUP SCRIPT AFTER THE 
# FIRST TIME YOU'VE DEALT WITH THIS
##########
```

Build GENIE:
```
# Get GENIE Generator from github
git clone git@github.com:GENIE-MC/Generator.git -b R-3_04_02

## configure
cd Generator/

./configure --prefix=${GENIE_FQ_DIR} \
--disable-profiler \
--disable-validation-tools \
--disable-cernlib \
--disable-lhapdf5 \
--enable-lhapdf6 \
--enable-gfortran \
--enable-flux-drivers \
--enable-geom-drivers \
--disable-doxygen \
--enable-test \
--enable-mueloss \
--enable-dylibversion \
--enable-t2k \
--enable-fnal \
--enable-atmo \
--enable-nucleon-decay \
--disable-masterclass \
--disable-debug \
--with-optimiz-level=O2

## make
make

## install
make install

## At this point, we have to remove "Generator/bin":
## "genie-config" firstly detects Generator/lib, and use that as library path.
## However, we want to use the ones that we "install"-ed,
rm -rf lib
```

To build `Reweight`, pay attention to getting the right branch `larsbp_feature_2p2h` from `larsb-p`:
```
cd ${GENIE_FQ_DIR}
git clone git@github.com:larsb-p/Reweight.git 
cd Reweight

# checkout the right tag for GENIE Reweight
# and stay in develop
git checkout -b larsbp_feature_2p2h_on_R-1_02_04 R-1_02_04

# merge lar's branch into this tag
git merge origin/larsbp_feature_2p2h

# make
make

# install
make install
```

Finally, go to your nusystematics installation, remove the build directory if present at all, and install again.
Note that for the MEC dials, you need the following NuSystematics branch: [nusystematics `feature/jskim_MECDecAngDep`](https://github.com/NuSystematics/nusystematics/tree/feature/jskim_MECDecAngDep)
```
# install NuSystematics
mkdir nusystematics; cd nusystematics
git clone git@github.com:NuSystematics/nusystematics.git nusystematics-src
cd nusystematics-src
git pull feature/jskim_MECDecAngDep

# go back and build
cd ..
mkdir build; cd build
cmake ../nusystematics-src/
make install
```

Remember to set the GENIE tune, even if I've observed it's not really needed anymore:
```
export GENIE_XSEC_TUNE=AR23_20i_00_000
```
and use the full configuration with MEC weights via `/exp/icarus/app/users/rtriozzi/All.ParamHeader.fcl`(or look at `numine/fcl/systs`):
```
UpdateReweight -c All.ParamHeader.fcl -i input_cafs.txt -o output_flat.caf.root
```

### PROfit

Everything about fitting is provided in the [PROfit repository](https://github.com/markrosslonergan/PROfit).
The `fitting/config` folder provides every configuration you could need via PROfit to perform the electron neutrino disappearance (and more) analysis, both with the simple disappearance-only model, and with 3+1.
The `fitting/datamc` provides an XML to perform data/MC comparisons with all the systematic uncertainties.
Note that everything was done with PROfit 2.4.0, but _do_ expect compatibility with newer tags (modulo some minor XML tweaks).

An exhaustive list of commands to proceed with the analysis via PROfit follows.

* I/O:
```
# create binary I/O from CAFAna trees
../Elephant_Vanishes/build/bin/PROfit -x /exp/icarus/app/users/rtriozzi/profit/Elephant_Vanishes/xml/PROfit_ICARUS-NuMI_PandoraNuGraph2_NuEDis.xml -t nues --log process.log process
```
* plotting, with systematics and covariance matrices, giving a ROOT file for pretty-plotting:
```
# plot variabiles w/ channels and systematics
../Elephant_Vanishes/build/bin/PROfit -x /exp/icarus/app/users/rtriozzi/profit/Elephant_Vanishes/xml/PROfit_ICARUS-NuMI_PandoraNuGraph2_NuEDis.xml -t nues --log plot.log plot
# --with-splines at the end to dump *all* the splines...

# inject CV in plot, e.g. for appearance studies (note different XML)
../Elephant_Vanishes/build/bin/PROfit -x /storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/profit/nueapp/PROfit_NuMI_nueapp.xml -t nue_app --inject-cv dmsq 2 sinsq2thee 0.2 sinsqth24 0.05 --log plot.log plot
```
* actually fitting, with profiling:
```
# global fit with an injected signal + Poisson throws
../Elephant_Vanishes/build/bin/PROfit -x /exp/icarus/app/users/rtriozzi/profit/Elephant_Vanishes/xml/PROfit_ICARUS-NuMI_PandoraNuGraph2_NuEDis.xml -t nues --poisson-throw --inject dmsq 1.5 sinsq2thee 0.3 --progress global

# profile over physics and syst parameters with an injected signal + Poisson throws
../Elephant_Vanishes/build/bin/PROfit -x /exp/icarus/app/users/rtriozzi/profit/Elephant_Vanishes/xml/PROfit_ICARUS-NuMI_PandoraNuGraph2_NuEDis.xml -t nues --poisson-throw --inject dmsq 1.5 sinsq2thee 0.3 --progress -n 8 profile
```
* asimov sensitivity for exclusion contours:
```
# surface
../Elephant_Vanishes/build/bin/PROfit -x /exp/icarus/app/users/rtriozzi/profit/Elephant_Vanishes/xml/PROfit_ICARUS-NuMI_PandoraNuGraph2_NuEDis.xml -t nues -o surfaces --log log.surf0 -v 2 -w 3 -n 8 surface -g 50 --xlo 0.01 --xhi 1 --ylo 0.1 --yhi 100

# surface, excluding, e.g., detector systematics
../Elephant_Vanishes/build/bin/PROfit -x /exp/icarus/app/users/rtriozzi/profit/Elephant_Vanishes/xml/PROfit_ICARUS-NuMI_PandoraNuGraph2_NuEDis.xml -t nues -o surfaces --log log.surf0 -v 2 -w 3 -n 8 --exclude-systs DetVar surface -g 50 --xlo 0.01 --xhi 1 --ylo 0.1 --yhi 100

# surface, stat-only errors
../Elephant_Vanishes/build/bin/PROfit -x /exp/icarus/app/users/rtriozzi/profit/Elephant_Vanishes/xml/PROfit_ICARUS-NuMI_PandoraNuGraph2_NuEDis.xml -t nues -o surfaces --log log.surf0 -v 2 -w 3 -n 8 --statonly surface -g 50 --xlo 0.01 --xhi 1 --ylo 0.1 --yhi 100
```
* ...and for discovery contours:
```
../../Elephant_Vanishes/build/bin/PROfit -x /storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/profit/nuedis/PROfit_NuMI_nuedis_DetSyst_Fixes.xml -t nue_detsysts -o surfaces_inject_gallium --inject dmsq 1.254237 sinsq2thee 0.3188 --log log_inject_gallium.surf0 -v 2 -w 3 -n 8 surface -g 50 --xlo 0.01 --xhi 1 --ylo 0.1 --yhi 100
```
* data/MC comparisons:
```
# process
../../Elephant_Vanishes/build/bin/PROfit -x /storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/profit/nuedis/PROfit_NuMI_MC_CompareWithData.xml --data /storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/profit/nuedis/PROfit_NuMI_data.xml

# plot
../../Elephant_Vanishes/build/bin/PROfit -x /storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/profit/nuedis/PROfit_NuMI_MC_CompareWithData.xml --data /storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/profit/nuedis/PROfit_NuMI_data.xml --log plot.log plot

# area-normalized plot
../../Elephant_Vanishes/build/bin/PROfit -x /storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/profit/nuedis/PROfit_NuMI_MC_CompareWithData.xml --data /storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/profit/nuedis/PROfit_NuMI_data.xml --area-norm --log plot.log plot
```

## ML

The ML pipeline is Python-based and self-contained. Set up the environment and
install the `src` package:
```
pip install -e ml/
```

Per-channel BDT trainings and evaluations live in the channel subdirectories as
Jupyter notebooks. Shared utilities (feature definitions, preprocessing, training,
evaluation, plotting) are in `ml/src/`.
