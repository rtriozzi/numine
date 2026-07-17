#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "CC1mu0pi_Cuts.h"
#include "CC1mu0pi_TruthCuts.h"
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

void CC1mu0pi_Debug() {

    const std::string TargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nuedis/data/202607/caf/caf/*/*Unblind.DONOTLOOK.flat.caf.root";

    SpectrumLoader NuLoader(TargetFile);

    // Spectrum *sMCEventDump = new Spectrum("", Binning::Simple(3, 0, 3), NuLoader, kEventTruthDump, kNoSpillCut);
    Spectrum *sMCEventDump = new Spectrum("", Binning::Simple(3, 0, 3), NuLoader, kDataDump, kNoSpillCut);

    NuLoader.Go();

    return;
}
