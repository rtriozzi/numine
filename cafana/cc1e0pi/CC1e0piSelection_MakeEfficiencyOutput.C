#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

// helpers
#include "CC1e0piSelection_Cuts.h"
#include "CC1e0piSelection_TruthCuts.h"
#include "CC1e0piSelection_Efficiency.h"

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

void CC1e0piSelection_MakeEfficiencyOutput() {

    // CNAF nuedis - nue-only flux
    SpectrumLoader NuLoader("/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap-cv-nueonly/run*/nuedis_cafmakerjob*/*.flat.caf.root"); // CV

    Spectrum sDummy("dummy", Binning::Simple(1, 0, 1), NuLoader, kCC1e0p1Signal_MakeSelectionEfficiencyOutput, kNoSpillCut);

    NuLoader.Go();

    return;
}
