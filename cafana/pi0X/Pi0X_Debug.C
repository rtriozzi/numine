#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "Pi0X_Cuts.h"
#include "Pi0X_TruthCuts.h"
#include "Debuggers.h"

// root stuff
#include "TCanvas.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "THStack.h"

using namespace ana;

void Pi0X_Debug() {

    SpectrumLoader data("/pnfs/icarus/scratch/users/rtriozzi/NuGraph_NuMIPrescaled_v10_06_00_01p01_1D_NuGraphReco_NueDis_ShowerCalo/caf/*/*Unblind.DONOTLOOK.flat.caf.root");

    Spectrum *sDataEventDump = new Spectrum("", Binning::Simple(3, 0, 3), data, kEventDump, kNoSpillCut);

    data.Go();

    return;
}
