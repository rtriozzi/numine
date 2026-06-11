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

// #include "sbnana/SBNAna/Cuts/NumuCutsIcarus202401.h"
// #include "sbnana/SBNAna/Cuts/ICARUSDataQualityCuts.h"
// #include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "NCPi0_Cuts.h"
#include "NCPi0_TruthCuts.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TDirectory.h"
#include "TFile.h"

#include <string>
#include <utility>
#include <vector>
#include <glob.h>

using namespace ana;

double offbeam_livetime = 0;
const SpillVar kOffbeamLivetime([](const caf::SRSpillProxy *sr) {
  if(kGoodRunCut(sr))
    offbeam_livetime += sr->hdr.noffbeamnumi;
  return 1;
});

void make_tree_offbeam(std::string outname = "CNAF_OffBeam_NCpi0_NuMI_NoSysts_ShowerCaloTools_EnergyFix.root")
{

  SpectrumLoader offbeam("/pnfs/icarus/scratch/users/rtriozzi/NuGraph_NuMIOffBeam_v10_06_00_01p01_1D_NuGraphReco_NueDis_ShowerCalo/caf/*/*Unblind.DONOTLOOK.flat.caf.root");
  
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
    "vtxx", "vtxy", "vtxz",
    "pi0gap", "pi0open", "pi0mass",
    "leadcollE", "leaddEdx",
    "sublcollE", "subldEdx"
  };

  std::vector<Var> vars = {
    kSlcVX, kSlcVY, kSlcVZ,
    kPi0_LargestConvGap, kPi0_CosPhotonOpenAngle, kPi0_InvariantMass,
    kLargestRecoShower_CollEnergy, kLargestRecoShower_ColldEdx,
    kSubleadRecoShower_CollEnergy, kSubleadRecoShower_ColldEdx
  };

  Tree offbeamtree("selectedOffbeam", branch_names, offbeam, vars, kSpillSelection && kGoodRunCut, kSliceSelection, kNoShift, true, true);
  Spectrum dummy_spec("", Binning::Simple(2,0,2), offbeam, kOffbeamLivetime, kNoSpillCut); 

  offbeam.Go();

  offbeamtree.OverrideLivetime(offbeam_livetime);

  TFile fout(outname.c_str(), "RECREATE");
  TDirectory *offdir = fout.mkdir("offbeam");
  offbeamtree.SaveTo(offdir);
}
