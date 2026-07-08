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
#include "CC1mu0pi_Cuts.h"
#include "CC1mu0pi_TruthCuts.h"

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

void make_tree_NoSysts_copy(std::string outname = "CNAF_CV_1muNp0pi_NuMI_var9_NoCut.root")
{
  // CNAF nuedis - nominal flux
  // SpectrumLoader mc("/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap-cv/run*/nuedis_cafmakerjob*/*.flat.caf.root");
  SpectrumLoader mc("/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap_variations/var9_hilifetime_fixed/run*/nuedis_cafmakerjob*/*.flat.caf.root");
  
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
    "trueE", "trueL", "truePDG", "CC", 
    "signal", "nue", "numu", "ispi0", "ischpi", "trueFV", "trueOOFV",
    "index",
    "selected",
    "recoE", "recopT_NuMI",
    "muonl", "muonp", "muonke",
    "muchi2mu", "muchi2pr",
    "leadpmom", "sleadpmom",
    "pchi2mu", "pchi2pr",
  };

  std::vector<Var> nu_vars = {
    kTrueE, kTrueL, kTruePDG, kTrueCC,
    static_cast<const Var>(kTrueCC1mu0pi), static_cast<const Var>(kIsNue), static_cast<const Var>(kIsNuMu), static_cast<const Var>(kIsTherePi0), static_cast<const Var>(kIsThereChPi), static_cast<const Var>(kTrueVertexInFV), static_cast<const Var>(kIsNuOOFV),
    kIndex,
    static_cast<const Var>(kAutomaticNuMuSelection),
    kRecoNeutrino_NuMuCC0piEnergy, kRecoNeutrino_NuMuCC0piTransverseMomentum_NuMI,
    kMuon_Length, kMuon_Momentum, kMuon_KE,
    kMuon_Chi2Muon, kMuon_Chi2Proton,
    kLeadingProtonMomentum, kSubLeadingProtonMomentum,
    kLeadingProton_Chi2Muon, kLeadingProton_Chi2Proton
  };

  // cosmics (MC and off-beam)
  std::vector<std::string> branch_names = {
    "index",
    "selected",
    "recoE", "recopT_NuMI",
    "muonl", "muonp", "muonke",
    "chi2mu", "chi2pr",
    "leadpmom", "sleadpmom",
    "pchi2mu", "pchi2pr",
  };

  std::vector<Var> vars = {
    kIndex,
    static_cast<const Var>(kAutomaticNuMuSelection),
    kRecoNeutrino_NuMuCC0piEnergy, kRecoNeutrino_NuMuCC0piTransverseMomentum_NuMI,
    kMuon_Length, kMuon_Momentum, kMuon_KE,
    kMuon_Chi2Muon, kMuon_Chi2Proton,
    kLeadingProtonMomentum, kSubLeadingProtonMomentum,
    kLeadingProton_Chi2Muon, kLeadingProton_Chi2Proton
  };

  Tree nutree("selectedNu", nu_branch_names, mc, nu_vars, kSpillSelection, kSliceSelection && kTrueNu, kNoShift, true, true);
  Tree costree("selectedCos", branch_names, mc, vars, kSpillSelection, kSliceSelection && !kTrueNu, kNoShift, true, true);

  mc.Go();

  TFile fout(outname.c_str(), "RECREATE");
  TDirectory* dir = fout.mkdir("events");
  nutree.SaveTo(dir); 
  costree.SaveTo(dir);
}
