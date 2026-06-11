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
  // if(icarus::kGoodRunsRun2(sr))
    offbeam_livetime += sr->hdr.noffbeamnumi;
  return 1;
});

void make_tree_NoSysts(std::string outname = "CNAF_CV_NCpi0_NuMI_NoSysts.root")
{

  SpectrumLoader mc("/pnfs/icarus/scratch/users/rtriozzi/NuGraph2/NueDis_CAFs_NuSystematics/caf_wMEC/*/*caf.root");
  
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

  // neutrino variables, including truth
  std::vector<std::string> nu_branch_names = {
    "trueE", "trueL", "truePDG", "CC", 
    "signal", "nue", "numu", "ispi0", "trueFV", "trueOOFV",
    "trueleadE", "truesublE",
    "index", 
    "pi0gap", "pi0open", "pi0mass",
    "leadcollE", "leaddEdx",
    "sublcollE", "subldEdx",
    "leadconvgap", "sublconvgap",
    "leadopenangle", "sublopenangle",
    "leadcosnuph", "sublcosnuph",
    "nshwpfps", "nshwpfpthr", "nmipfps",
    "np", "notherpfps", "leadpmom", "sublpmom"
  };

  std::vector<Var> nu_vars = {
    kTrueE, kTrueL, kTruePDG, kTrueCC, 
    static_cast<const Var>(kTrueNCPi0), static_cast<const Var>(kIsNue), static_cast<const Var>(kIsNuMu), static_cast<const Var>(kIsTherePi0), static_cast<const Var>(kTrueVertexInFV), static_cast<const Var>(kIsNuOOFV),
    kLargestRecoShower_TrueEnergy, kSubleadRecoShower_TrueEnergy,
    kIndex, 
    kPi0_LargestConvGap, kPi0_CosPhotonOpenAngle, kPi0_InvariantMass,
    kLargestRecoShower_CollEnergy, kLargestRecoShower_ColldEdx,
    kSubleadRecoShower_CollEnergy, kSubleadRecoShower_ColldEdx,
    kLargestRecoShower_ConvGap, kSubleadRecoShower_ConvGap,
    kLargestRecoShower_OpenAngle, kSubleadRecoShower_OpenAngle,
    kLargestRecoShower_CosVertexAngle, kSubleadRecoShower_CosVertexAngle,
    kNuGraph_NShowerPFPs, kNuGraph_NShowerPFPs_AboveThreshold, kNuGraph_NMIPPFPs,
    kNSelectedProtons_N, kNSelectedProtonsIdx_NOtherParticles, kLeadingProtonMomentum, kSubLeadingProtonMomentum
  };

  // cosmics (MC and off-beam)
  std::vector<std::string> branch_names = {
    "index", 
    "pi0gap", "pi0open", "pi0mass",
    "leadcollE", "leaddEdx",
    "sublcollE", "subldEdx",
    "leadconvgap", "sublconvgap",
    "leadopenangle", "sublopenangle",
    "leadcosnuph", "sublcosnuph",
    "nshwpfps", "nshwpfpthr", "nmipfps", 
    "np", "notherpfps", "leadpmom", "sublpmom"
  };

  std::vector<Var> vars = {
    kIndex, 
    kPi0_LargestConvGap, kPi0_CosPhotonOpenAngle, kPi0_InvariantMass,
    kLargestRecoShower_CollEnergy, kLargestRecoShower_ColldEdx,
    kSubleadRecoShower_CollEnergy, kSubleadRecoShower_ColldEdx,
    kLargestRecoShower_ConvGap, kSubleadRecoShower_ConvGap,
    kLargestRecoShower_OpenAngle, kSubleadRecoShower_OpenAngle,
    kLargestRecoShower_CosVertexAngle, kSubleadRecoShower_CosVertexAngle,
    kNuGraph_NShowerPFPs, kNuGraph_NShowerPFPs_AboveThreshold, kNuGraph_NMIPPFPs,
    kNSelectedProtons_N, kNSelectedProtonsIdx_NOtherParticles, kLeadingProtonMomentum, kSubLeadingProtonMomentum
  };

  Tree nutree("selectedNu", nu_branch_names, mc, nu_vars, kSpillSelection, kSliceSelection && kTrueNu, kNoShift, true, true);
  Tree costree("selectedCos", branch_names, mc, vars, kSpillSelection, kSliceSelection && !kTrueNu, kNoShift, true, true);

  mc.Go();

  TFile fout(outname.c_str(), "RECREATE");
  TDirectory* dir = fout.mkdir("events");
  nutree.SaveTo(dir); 
  costree.SaveTo(dir);
}
