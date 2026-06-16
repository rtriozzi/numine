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

    const std::string TargetFile = "/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap-cv/run*/nuedis_cafmakerjob*/*.flat.caf.root";

    SpectrumLoader NuLoader(TargetFile);

    Spectrum *sMCEventDump = new Spectrum("", Binning::Simple(3, 0, 3), NuLoader, kEventTruthDump, kNoSpillCut);

    NuLoader.Go();

    return;
}
