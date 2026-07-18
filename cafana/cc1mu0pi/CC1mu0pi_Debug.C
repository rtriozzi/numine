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

    // FNAL data
    // const std::string TargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nuedis/data/202607/caf/caf/*/*Unblind.DONOTLOOK.flat.caf.root";

    // CNAF MC
    const std::string TargetFile = "/storage/gpfs_data/icarus/local/users/cfarnese/Produzioni_Riccardo_NUMInue_2026/caf_wMEC/*/*flat.caf.root";

    SpectrumLoader NuLoader(TargetFile);

    // truth dumping
    // Spectrum *sMCEventDump = new Spectrum("", Binning::Simple(3, 0, 3), NuLoader, kEventTruthDump, kNoSpillCut);

    // data dumping
    // Spectrum *sMCEventDump = new Spectrum("", Binning::Simple(3, 0, 3), NuLoader, kDataDump, kNoSpillCut);

    // muons crossing boundaries
    Spectrum *sMCEventDump = new Spectrum("", Binning::Simple(4, 0, 4), NuLoader, kMuonCrossingDump, kNoSpillCut);

    NuLoader.Go();

    // muons crossing boundaries
    std::cout << "No-crossing "         << sMCEventDump->ToTH1(sMCEventDump->POT())->GetBinContent(1) << std::endl;
    std::cout << "Cathode-crossing: "   << sMCEventDump->ToTH1(sMCEventDump->POT())->GetBinContent(2) << std::endl;
    std::cout << "z=0-crossing: "       << sMCEventDump->ToTH1(sMCEventDump->POT())->GetBinContent(3) << std::endl;
    std::cout << "Both crossings: "     << sMCEventDump->ToTH1(sMCEventDump->POT())->GetBinContent(4) << std::endl;

    return;
}
