#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

// helpers
#include "CC1mu0pi_Cuts.h"
#include "CC1mu0pi_TruthCuts.h"
#include "CC1mu0pi_Efficiency.h"

// ROOT stuff
#include "TCanvas.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "THStack.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"

using namespace ana;

void CC1mu0pi_MakeEfficiencyOutput() {

    // CNAF nuedis - nom flux
    SpectrumLoader NuLoader("/storage/gpfs_data/icarus/local/users/cfarnese/Produzioni_Riccardo_NUMInue_2026/caf_wMEC/*/*flat.caf.root");

    Spectrum sDummy("dummy", Binning::Simple(1, 0, 1), NuLoader, kCC1mu0piSignal_MakeSelectionEfficiencyOutput, kNoSpillCut);

    NuLoader.Go();

    return;
}
