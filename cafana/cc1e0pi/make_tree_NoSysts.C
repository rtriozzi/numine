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

void make_tree_NoSysts(std::string outname = "CNAF_CV_1eNp0pi_NuMI_NoSysts_PreselectionElectron.root")
{
  // CNAF nuedis - nominal flux
  SpectrumLoader mc("/storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/concats/cv/cv_run*.flat.caf.root");
  // SpectrumLoader mc("/storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/concats/var1_hitcohnoise/var1_run*.flat.caf.root");
  // SpectrumLoader mc("/storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/concats/var2_hiintnoise/var2_run*.flat.caf.root");
  
  // CNAF nuedis - nue-only flux
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
  
  // SpectrumLoader offbeam("/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco/numi_offbeam.unblind.flat.caf.root");
  
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
  // const Cut kSliceSelection = kAutomaticSelection;
  const Cut kSliceSelection = kPreSelectionElectron;

  // neutrino variables, including truth
  std::vector<std::string> nu_branch_names = {
    "trueE", "trueL", "truePDG", "CC", 
    "signal", "nue", "numu", "ispi0", "trueFV", "trueOOFV",
    "index", "recoE", 
    "deltaZ", "deltaZ_Trigger", "flashTime",
    "collE", "colldEdx", "openangle", "convgap"
  };

  std::vector<Var> nu_vars = {
    kTrueE, kTrueL, kTruePDG, kTrueCC, 
    static_cast<const Var>(kTrueCC1e0pi), static_cast<const Var>(kIsNue), static_cast<const Var>(kIsNuMu), static_cast<const Var>(kIsTherePi0), static_cast<const Var>(kTrueVertexInFV), static_cast<const Var>(kIsNuOOFV),
    kIndex, kRecoNeutrino_CC0piEnergy, 
    kBarycenterFM_DeltaZ, kBarycenterFM_DeltaZ_Trigger, kBarycenterFM_FlashTime,
    kLargestRecoShower_CollEnergy, kLargestRecoShower_ColldEdx, kLargestRecoShower_OpenAngle, kLargestRecoShower_ConvGap
  };

  // cosmics (MC and off-beam)
  std::vector<std::string> branch_names = {
    "index", "recoE",
    "deltaZ", "deltaZ_Trigger", "flashTime",
    "collE", "colldEdx", "openangle", "convgap"
  };

  std::vector<Var> vars = {
    kIndex, kRecoNeutrino_CC0piEnergy,
    kBarycenterFM_DeltaZ, kBarycenterFM_DeltaZ_Trigger, kBarycenterFM_FlashTime,
    kLargestRecoShower_CollEnergy, kLargestRecoShower_ColldEdx, kLargestRecoShower_OpenAngle, kLargestRecoShower_ConvGap
  };

  Tree nutree("selectedNu", nu_branch_names, mc, nu_vars, kSpillSelection, kSliceSelection && kTrueNu, kNoShift, true, true);
  Tree costree("selectedCos", branch_names, mc, vars, kSpillSelection, kSliceSelection && !kTrueNu, kNoShift, true, true);
  // Tree offbeamtree("selectedOffbeam", branch_names, offbeam, vars, kSpillSelection, kSliceSelection, kNoShift, true, true);
  // Spectrum dummy_spec("", Binning::Simple(2,0,2), offbeam, kOffbeamLivetime, kNoSpillCut); 

  mc.Go();
  // offbeam.Go();

  // offbeamtree.OverrideLivetime(offbeam_livetime);

  TFile fout(outname.c_str(), "RECREATE");
  TDirectory* dir = fout.mkdir("events");
  nutree.SaveTo(dir); 
  costree.SaveTo(dir);
  // TDirectory *offdir = fout.mkdir("offbeam");
  // offbeamtree.SaveTo(offdir);
}
