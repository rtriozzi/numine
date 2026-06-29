#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Systs/SBNOnOffSysts.h"
#include "sbnana/CAFAna/Systs/UniverseOracle.h"
// #include "sbnana/CAFAna/Systs/IcarusRun2DetectorSysts.h"

#include "sbnana/CAFAna/Systs/NuMIFluxSysts.h"

#include "sbnana/SBNAna/Vars/Vars.h"
#include "CC1e0piSelection_Cuts.h"
#include "CC1e0piSelection_TruthCuts.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TDirectory.h"
#include "TFile.h"

#include <string>
#include <utility>
#include <vector>
#include <glob.h>

using namespace ana;

double data_livetime = 0;
const SpillVar kDataLivetime([](const caf::SRSpillProxy *sr) {
  if(kGoodRunCut(sr))
    data_livetime += sr->hdr.numiinfo.size();
  return 1;
});

void make_tree_data(std::string outname = "CNAF_Data_1eNp0pi_NuMI.root")
{

  // data
  SpectrumLoader data("/pnfs/icarus/persistent/users/rtriozzi/nuedis/data/202607/caf/caf/*/*Unblind.DONOTLOOK.flat.caf.root");

  // some simple truth variables on the fly
  const Var kTrueE = SIMPLEVAR(truth.E);
  const Var kTrueL = SIMPLEVAR(truth.baseline);
  const Var kTruePDG = SIMPLEVAR(truth.pdg);
  const Var kTrueCC = SIMPLEVAR(truth.iscc);
  const Var kIndex = SIMPLEVAR(truth.index);
  const Var kSlcVX = SIMPLEVAR(vertex.x);
  const Var kSlcVY = SIMPLEVAR(vertex.y);
  const Var kSlcVZ = SIMPLEVAR(vertex.z);

  const Cut kTrueNu = SIMPLEVAR(truth.index) >= 0;

  // event selection
  const SpillCut kSpillSelection = kNoSpillCut;
  const Cut kSliceSelection = kAutomaticSelection;

  // variables
  std::vector<std::string> branch_names = {
    "index",
    // neutrino
    "recoE", "recopT", "recopT_NuMI",
    "recoepT", "recoepT_NuMI", "recoppT", "recoppT_NuMI",
    "direp3d", "direpT", "direpT_NuMI",
    // event
    "vtxx", "vtxy", "vtxz",
    // electron
    "elrecoE", "gap", "angle", "colldEdx", "elength",
    "eldirx", "eldiry", "eldirz", "eldirnumi",
    "elendx", "elendy", "elendz",
    // proton
    "np", "pmom", "slpmom",
    "pendx", "pendy", "pendz",
    "pdirx", "pdiry", "pdirz", "pdirnumi"
  };

  std::vector<Var> vars = {
    kIndex,
    // neutrino
    kRecoNeutrino_CC0piEnergy, kRecoNeutrino_CC0piTransverseMomentum, kRecoNeutrino_CC0piTransverseMomentum_NuMI,
    kRecoNeutrino_ElectronTransverseMomentum, kRecoNeutrino_ElectronTransverseMomentum_NuMI, kRecoNeutrino_ProtonTransverseMomentum, kRecoNeutrino_ProtonTransverseMomentum_NuMI,
    kRecoNeutrino_epCosAngle_3D, kRecoNeutrino_epCosAngle_Transverse, kRecoNeutrino_epCosAngle_Transverse_NuMI,
    // event
    kSlcVX, kSlcVY, kSlcVZ, 
    // electron
    kLargestRecoShower_CollEnergy, kLargestRecoShower_ConvGap, kLargestRecoShower_OpenAngle, kLargestRecoShower_ColldEdx, kLargestRecoShower_Length,
    kLargestRecoShower_DirX, kLargestRecoShower_DirY, kLargestRecoShower_DirZ, kLargestRecoShower_DirNuMI,
    kLargestRecoShower_EndX, kLargestRecoShower_EndY, kLargestRecoShower_EndZ,
    // proton
    kNSelectedProtons_N, kLeadingProtonMomentum, kSubLeadingProtonMomentum,
    kLeadingProton_EndX, kLeadingProton_EndY, kLeadingProton_EndZ,
    kLeadingProton_DirX, kLeadingProton_DirY, kLeadingProton_DirZ, kLeadingProton_DirNuMI,
  };

  Tree nutree("selectedData", branch_names, data, vars, kSpillSelection && kNuMISpillQualityCut && kGoodRunCut, kSliceSelection, kNoShift, true, true);
  Spectrum dummy_spec_data("", Binning::Simple(2,0,2), data, kDataLivetime, kNoSpillCut);

  data.Go();

  nutree.OverrideLivetime(data_livetime);

  TFile fout(outname.c_str(), "RECREATE");
  TDirectory* dir = fout.mkdir("events");
  nutree.SaveTo(dir); 
}
