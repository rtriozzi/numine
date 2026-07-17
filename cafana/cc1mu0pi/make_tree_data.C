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

double data_livetime = 0;
const SpillVar kDataLivetime([](const caf::SRSpillProxy *sr) {
  if(kGoodRunCut(sr))
    data_livetime += sr->hdr.numiinfo.size();
  return 1;
});

void make_tree_data(std::string outname = "CNAF_Data_1muNp0pi_NuMI.root")
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
  const Cut kSliceSelection = kAutomaticNuMuSelection;

  // variables
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

  Tree nutree("selectedData", branch_names, data, vars, kSpillSelection && kNuMISpillQualityCut && kGoodRunCut, kSliceSelection, kNoShift, true, true);
  Spectrum dummy_spec_data("", Binning::Simple(2,0,2), data, kDataLivetime, kNoSpillCut);

  data.Go();

  nutree.OverrideLivetime(data_livetime);

  TFile fout(outname.c_str(), "RECREATE");
  TDirectory* dir = fout.mkdir("events");
  nutree.SaveTo(dir); 
}
