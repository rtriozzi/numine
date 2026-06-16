#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

// helpers
#include "Pi0X_Cuts.h"
#include "Pi0X_TruthCuts.h"
#include "Debuggers.h"

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

void Pi0X_MakeEfficiencyOutput() {

    // CNAF nuedis - nue-only flux
    SpectrumLoader mc("/pnfs/icarus/scratch/users/rtriozzi/NuGraph2/NueDis_CAFs_NuSystematics/caf_wMEC/*/*caf.root");

    Spectrum sDummy("dummy", Binning::Simple(1, 0, 1), mc, kPi0X_MakeSelectionEfficiencyOutput, kNoSpillCut);

    mc.Go();

    return;
}
