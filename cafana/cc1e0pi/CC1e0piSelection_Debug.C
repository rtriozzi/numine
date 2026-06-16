#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "CC1e0piSelection_Cuts.h"
#include "CC1e0piSelection_TruthCuts.h"
#include "CC1e0piSelection_Efficiency.h"
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

void CC1e0piSelection_Debug() {

    SpectrumLoader data("/pnfs/icarus/persistent/users/rtriozzi/nuedis/data/FNAL_NuMI_Run2_PrescaledData_20260403.flat.caf.root");

    Spectrum *sDataEventDump = new Spectrum("", Binning::Simple(3, 0, 3), data, kNuMIBeamQuality, kNoSpillCut);

    data.Go();

    return;
}
