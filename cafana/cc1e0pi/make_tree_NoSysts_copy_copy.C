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

double offbeam_livetime = 0;
const SpillVar kOffbeamLivetime([](const caf::SRSpillProxy *sr) {
  // if(icarus::kGoodRunsRun2(sr))
    offbeam_livetime += sr->hdr.noffbeambnb;
  return 1;
});

std::vector<std::string> expand_glob(const std::string& pattern) {
    glob_t result;
    std::vector<std::string> files;
    if (glob(pattern.c_str(), GLOB_TILDE, nullptr, &result) == 0) {
        for (size_t i = 0; i < result.gl_pathc; ++i)
            files.push_back(result.gl_pathv[i]);
    }
    globfree(&result);
    return files;
}

void make_tree_NoSysts_copy_copy(std::string outname = "CNAF_NuE_1eNp0pi_NuMI_var14_NoCut.root")
{
  // CNAF nuedis - nominal flux
  // SpectrumLoader mc("/storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/concats/cv/cv_run*.flat.caf.root");
  // SpectrumLoader mc("/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap-cv/run*/nuedis_cafmakerjob*/*.flat.caf.root");
  // SpectrumLoader mc("/storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/concats/var1_hitcohnoise/var1_run*.flat.caf.root");
  // SpectrumLoader mc("/storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/concats/var2_hiintnoise/var2_run*.flat.caf.root");
  
  // CNAF nuedis - nue-only flux
  // SpectrumLoader mc("/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap-cv-nueonly/run*/nuedis_cafmakerjob*/*.flat.caf.root"); // CV
  SpectrumLoader mc("/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap_variations/nue_var14_opaque/run*/nuedis_cafmakerjob*/*.flat.caf.root"); // CV

  // std::vector<std::string> files;
  // for (int run = 1; run <= 2100; run++) {
  //   // auto expanded = expand_glob("/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap-cv-nueonly/run"+std::to_string(run)+"/nuedis_cafmakerjob*/*.flat.caf.root");
  //   // auto expanded = expand_glob("/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap_variations/nue_var2_hiintnoise/run"+std::to_string(run)+"/nuedis_cafmakerjob*/*.flat.caf.root");
  //   // auto expanded = expand_glob("/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap_variations/nue_var3_recomb/run"+std::to_string(run)+"/nuedis_cafmakerjob*/*.flat.caf.root");
  //   auto expanded = expand_glob("/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap_variations/nue_var4_diff/run"+std::to_string(run)+"/nuedis_cafmakerjob*/*.flat.caf.root");
  //   // auto expanded = expand_glob("/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap_variations/nue_var5_null/run"+std::to_string(run)+"/nuedis_cafmakerjob*/*.flat.caf.root");
  //   files.insert(files.end(), expanded.begin(), expanded.end());
  // }
  // SpectrumLoader mc(files);
    
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
  const Cut kSliceSelection = kNoCut;

  // neutrino variables, including truth
  std::vector<std::string> nu_branch_names = {
    "signal", "nue", "numu", "ispi0", "trueFV", "trueOOFV",
    "trueE", "truevisE", "trueL", "truePDG", "CC", "index", "selected",
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

  std::vector<Var> nu_vars = {
    static_cast<const Var>(kTrueCC1e0pi), static_cast<const Var>(kIsNue), static_cast<const Var>(kIsNuMu), static_cast<const Var>(kIsTherePi0), static_cast<const Var>(kTrueVertexInFV), static_cast<const Var>(kIsNuOOFV),
    kTrueE, kTrueCC1e0piVisibleEnergy, kTrueL, kTruePDG, kTrueCC, kIndex, static_cast<Var>(kAutomaticSelection),
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

  // cosmics (MC and off-beam)
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

  Tree nutree("selectedNu", nu_branch_names, mc, nu_vars, kSpillSelection, kSliceSelection && kTrueNu, kNoShift, true, true);
  Tree costree("selectedCos", branch_names, mc, vars, kSpillSelection, kSliceSelection && !kTrueNu, kNoShift, true, true);

  mc.Go();

  TFile fout(outname.c_str(), "RECREATE");
  TDirectory* dir = fout.mkdir("events");
  nutree.SaveTo(dir); 
  costree.SaveTo(dir);
}
