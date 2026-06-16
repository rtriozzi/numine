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

    // // FNAL
    // const std::string TargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nuedis/data/FNAL_NuMI_Run2_PrescaledData_20260403.flat.caf.root";

    // CV nominal
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap-cv/run*/nuedis_cafmakerjob*/*.flat.caf.root";

    // CV nue-only 
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap-cv-nueonly/run*/nuedis_cafmakerjob*/*.flat.caf.root"; 
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap_variations/nue_var9_hilifetime/run*/nuedis_cafmakerjob*/*.flat.caf.root";
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap_variations/nue_var9_hilifetime_test/run*/nuedis_cafmakerjob*/*.flat.caf.root";
    const std::string TargetFile = "/storage/gpfs_data/icarus/plain/data/mc/mc-v10_06_00_01p01-202603-cnaf-numi-nue-disap_variations/nue_var9_hilifetime_fixed/run*/nuedis_cafmakerjob*/*.flat.caf.root";

    SpectrumLoader NuLoader(TargetFile);

    // Spectrum *sMCEventDump = new Spectrum("", Binning::Simple(3, 0, 3), NuLoader, kEventTruthDump, kNoSpillCut);
    Spectrum *sMCEventDump = new Spectrum("", Binning::Simple(3, 0, 3), NuLoader, kElectronStudy, kNoSpillCut);

    NuLoader.Go();

    return;
}
