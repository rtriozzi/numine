#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

// #include "QELikeNuESelection_Vars.h"
#include "CC1e0piSelection_Cuts.h"
#include "CC1e0piSelection_TruthCuts.h"
#include "CC1e0piSelection_Efficiency.h"

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

void CC1e0piSelection() {

    // CNAF NuMI MC
    const std::string TargetFile = "/storage/gpfs_data/icarus/local/users/cfarnese/NUMI/NUMI_MC/0.root";

    SpectrumLoader NuLoader(TargetFile);

    Spectrum *sCounting = new Spectrum("", 
                                      Binning::Simple(3, 0, 3),
                                      NuLoader,
                                      kCounting, 
                                      kCRTPMTNeutrino,
                                      kNotClearCosmic); 

    NuLoader.Go();

    TFile FOut("CC1e0piSelection.root", "recreate");

    TCanvas *c = new TCanvas("counting", "counting");
    double TargetPOT = sCounting->POT();
    TH1 *h = sCounting->ToTH1(TargetPOT);
    h->Draw("HISTE0");
    c->Write();

    FOut.Close();

    return;
}