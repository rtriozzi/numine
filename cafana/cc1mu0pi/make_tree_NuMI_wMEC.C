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

std::vector<std::string> bound_at_neg2{
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NormCCMEC",
        "GENIEReWeight_SBN_v1_multisigma_NormNCMEC",
};

double offbeam_livetime = 0;
const SpillVar kOffbeamLivetime([](const caf::SRSpillProxy *sr) {
  // if(icarus::kGoodRunsRun2(sr))
    offbeam_livetime += sr->hdr.noffbeambnb;
  return 1;
});

const Var kDebugWgts([](const caf::SRSliceProxy* slc) -> double {
  bool bad = slc->truth.wgt.size() > 147 && slc->truth.wgt[147].univ.size() < 6;
  if(bad) {
    std::cout << "BAD EVENT: nwgts=" << slc->truth.wgt.size() << "\n";
    for(size_t i = 0; i < slc->truth.wgt.size(); ++i)
      std::cout << "  wgt[" << i << "].univ.size() = " 
                << slc->truth.wgt[i].univ.size() << "\n";
  }
  return 1.0;
});

void make_tree_NuMI_wMEC(std::string outname = "CNAF_CV_1muNp0pi_NuMI_wMEC.root")
{
  // FNAL
  // SpectrumLoader mc("/pnfs/icarus/persistent/users/rtriozzi/nuedis/MC_wMEC/*flat.caf.root");

  // CNAF
  SpectrumLoader mc("/storage/gpfs_data/icarus/local/users/cfarnese/Produzioni_Riccardo_NUMInue_2026/caf_wMEC/*/*flat.caf.root");

  // some simple truth variables on the fly
  const Var kTrueE = SIMPLEVAR(truth.E);
  const Var kTrueL = SIMPLEVAR(truth.baseline);
  const Var kTruePDG = SIMPLEVAR(truth.pdg);
  const Var kTrueCC = SIMPLEVAR(truth.iscc);
  const Var kIndex = SIMPLEVAR(truth.index);
  const Var kMode = SIMPLEVAR(truth.genie_mode);
  const Var kSlcVX = SIMPLEVAR(vertex.x);
  const Var kSlcVY = SIMPLEVAR(vertex.y);
  const Var kSlcVZ = SIMPLEVAR(vertex.z);

  const Cut kTrueNu = SIMPLEVAR(truth.index) >= 0;

  // event selection
  const SpillCut kSpillSelection = kNoSpillCut;
  const Cut kSliceSelection = kAutomaticNuMuSelection;

  // neutrino variables, including truth
  std::vector<std::string> nu_branch_names = {
    "trueE", "trueL", "truePDG", "CC", 
    "signal", "nue", "numu", "ispi0", "ischpi", "trueFV", "trueOOFV",
    "index",
    "selected",
    "recoE", "recopT_NuMI",
    "muonl", "muonp", "muonke",
    "truemuonl", "truemuonke",
    "muchi2mu", "muchi2pr",
    "mutrkscore", "muendx", "muendy", "muendz",
    "nprotons",
    "leadpmom", "sleadpmom",
    "pchi2mu", "pchi2pr",
    "ptrkscore",
  };

  std::vector<Var> nu_vars = {
    kTrueE, kTrueL, kTruePDG, kTrueCC,
    static_cast<const Var>(kTrueCC1mu0pi), static_cast<const Var>(kIsNue), static_cast<const Var>(kIsNuMu), static_cast<const Var>(kIsTherePi0), static_cast<const Var>(kIsThereChPi), static_cast<const Var>(kTrueVertexInFV), static_cast<const Var>(kIsNuOOFV),
    kIndex,
    static_cast<const Var>(kAutomaticNuMuSelection),
    kRecoNeutrino_NuMuCC0piEnergy, kRecoNeutrino_NuMuCC0piTransverseMomentum_NuMI,
    kMuon_Length, kMuon_Momentum, kMuon_KE,
    kMuon_TrueLength, kMuon_TrueKE,
    kMuon_Chi2Muon, kMuon_Chi2Proton,
    kMuon_TrackScore, kMuon_EndX, kMuon_EndY, kMuon_EndZ,
    kNSelectedProtons_N,
    kLeadingProtonMomentum, kSubLeadingProtonMomentum,
    kLeadingProton_Chi2Muon, kLeadingProton_Chi2Proton,
    kLeadingProton_TrackScore
  };

  // cosmics (MC and off-beam)
  std::vector<std::string> branch_names = {
    "index",
    "selected",
    "recoE", "recopT_NuMI",
    "muonl", "muonp", "muonke",
    "truemuonl", "truemuonke",
    "muchi2mu", "muchi2pr",
    "mutrkscore", "muendx", "muendy", "muendz",
    "nprotons",
    "leadpmom", "sleadpmom",
    "pchi2mu", "pchi2pr",
    "ptrkscore",
  };

  std::vector<Var> vars = {
    kIndex,
    static_cast<const Var>(kAutomaticNuMuSelection),
    kRecoNeutrino_NuMuCC0piEnergy, kRecoNeutrino_NuMuCC0piTransverseMomentum_NuMI,
    kMuon_Length, kMuon_Momentum, kMuon_KE,
    kMuon_TrueLength, kMuon_TrueKE,
    kMuon_Chi2Muon, kMuon_Chi2Proton,
    kMuon_TrackScore, kMuon_EndX, kMuon_EndY, kMuon_EndZ,
    kNSelectedProtons_N,
    kLeadingProtonMomentum, kSubLeadingProtonMomentum,
    kLeadingProton_Chi2Muon, kLeadingProton_Chi2Proton,
    kLeadingProton_TrackScore
  };

  Tree nutree("selectedNu", nu_branch_names, mc, nu_vars, kSpillSelection, kSliceSelection && kTrueNu, kNoShift, true, true);
  Tree costree("selectedCos", branch_names, mc, vars, kSpillSelection, kSliceSelection && !kTrueNu, kNoShift, true, true);

  std::vector<std::string> genie_names = GetSBNGenieWeightNames();
  genie_names.push_back("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b1");
  genie_names.push_back("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b2");
  genie_names.push_back("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b3");
  genie_names.push_back("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b4");
  genie_names.push_back("GENIEReWeight_SBNNuSyst_GENIE_multisigma_FracPN_CCMEC");
  SBNWeightSyst b1("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b1");
  SBNWeightSyst b2("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b2");
  SBNWeightSyst b3("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b3");
  SBNWeightSyst b4("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b4");
  SBNWeightSyst fracPN("GENIEReWeight_SBNNuSyst_GENIE_multisigma_FracPN_CCMEC");
  std::vector<const ISyst*> genie_systs = GetSBNGenieWeightSysts();
  genie_systs.push_back(&b1);
  genie_systs.push_back(&b2);
  genie_systs.push_back(&b3);
  genie_systs.push_back(&b4);
  genie_systs.push_back(&fracPN);
  std::vector<std::string> onoff_names = GetSBNOnOffNames();
  std::vector<const ISyst*> onoff_systs = GetSBNOnOffSysts();
  onoff_names.push_back("GENIEReWeight_SBNNuSyst_GENIE_multisigma_XSecShape_CCMEC");
  SBNOnOffSyst mecshape("GENIEReWeight_SBNNuSyst_GENIE_multisigma_XSecShape_CCMEC");
  onoff_systs.push_back(&mecshape);
  std::vector<std::string> nsigma_names = genie_names;
  std::vector<const ISyst*> nsigma_systs = genie_systs;
  for(const std::string &name: onoff_names) nsigma_names.push_back(name);
  for(const ISyst *s: onoff_systs) nsigma_systs.push_back(s);
  std::vector<std::pair<int, int>> min_max;
  for(size_t i = 0; i < genie_names.size(); ++i) min_max.emplace_back(-3, 3);
  for(size_t i = 0; i < onoff_names.size(); ++i) min_max.emplace_back(0, 1);
  for(size_t i = 0; i < genie_names.size(); ++i) {
    if(std::find(bound_at_neg2.begin(), bound_at_neg2.end(), genie_names[i]) != bound_at_neg2.end())
      min_max[i] = std::make_pair(-2, 3);
  }

  // NuMI beamline flux systs
  std::vector<const ISyst*> numi_beamline_systs = GetNuMIBeamlineFluxSysts();
  for(const ISyst* s : numi_beamline_systs) {
    nsigma_names.push_back(s->ShortName());
    nsigma_systs.push_back(s);
    min_max.emplace_back(-3, 3);
  }

  // NuMI PCA flux systs
  std::vector<const ISyst*> numi_pca_systs = GetNuMIPCAFluxSysts(15);
  for(const ISyst* s : numi_pca_systs) {
    nsigma_names.push_back(s->ShortName());
    nsigma_systs.push_back(s);
    min_max.emplace_back(-3, 3);
  }

  NSigmasTree nsigtree("multisigmaTree", nsigma_names, mc, nsigma_systs, min_max, kSpillSelection, kSliceSelection && kTrueNu, kNoShift, true, true);
  
  const std::vector<std::string> xsec_multisim_names{
    "GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_RPA_CCQE",
    "GENIEReWeight_SBN_v1_multisim_CoulombCCQE",
    "GENIEReWeight_SBN_v1_multisim_NormCCMEC",
    "GENIEReWeight_SBN_v1_multisim_NormNCMEC",
    "GENIEReWeight_SBN_v1_multisim_NCELVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi",
    "GENIEReWeight_SBN_v1_multisim_COHVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse",
    "GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse",
  };
  std::vector<std::string> multisim_names;
  std::vector<std::vector<Var>> univsKnobs;
  std::vector<unsigned int> nuniverses;

  for(const auto& name: xsec_multisim_names) {
    multisim_names.push_back(name);
    size_t nuniv = 100;
    nuniverses.push_back(nuniv);
    univsKnobs.emplace_back();
    for(size_t i = 0; i < nuniv; ++i) 
      univsKnobs.back().push_back(GetUniverseWeight(name, i));
  }

  NUniversesTree nunivtree("multisimTree", multisim_names, mc, univsKnobs, nuniverses, kSpillSelection, kSliceSelection && kTrueNu, kNoShift, true, true);

  mc.Go();

  // nsigtree.MergeTree(nutree);
  // nunivtree.MergeTree(nutree);

  TFile fout(outname.c_str(), "RECREATE");
  TDirectory* dir = fout.mkdir("events");
  nutree.SaveTo(dir); 
  costree.SaveTo(dir);
  nsigtree.SaveTo(dir);
  nunivtree.SaveTo(dir);
}
