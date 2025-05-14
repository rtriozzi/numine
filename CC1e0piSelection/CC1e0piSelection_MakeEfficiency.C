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

using namespace ana;

void CC1e0piSelection_MakeEfficiency() {

    // CNAF NuMI MC
    const std::string TargetFile = "/storage/gpfs_data/icarus/local/users/cfarnese/NUMI/NUMI_MC/*.root";

    SpectrumLoader NuLoader(TargetFile);

    // selection efficiency
    std::vector<double> TrueEnergyBinning = {0., 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2., 2.4, 2.8, 3.2, 4};

    Spectrum* sTrueNeutrinoEnergy = new Spectrum("E_{#nu} [GeV]", 
                                                 Binning::Custom(TrueEnergyBinning), //Binning::Simple(20, 0, 4),
                                                 NuLoader,
                                                 kCC1e0p1Signal_TrueNeutrinoEnergy, 
                                                 kNoSpillCut);

    Spectrum* sTrueNeutrinoEnergy_CRTPMT = new Spectrum("E_{#nu} [GeV]", 
                                                 Binning::Custom(TrueEnergyBinning),
                                                 NuLoader,
                                                 kCC1e0p1Signal_TrueNeutrinoEnergy, 
                                                 kCRTPMTNeutrino);

    const unsigned int kNSelectionSteps = SelectionSteps.size();
    Spectrum *sTrueNeutrinoEnergy_SelectionSteps[kNSelectionSteps];

    for(unsigned int iSel = 0; iSel < kNSelectionSteps; ++iSel) {
        sTrueNeutrinoEnergy_SelectionSteps[iSel] = new Spectrum("E_{#nu} [GeV]", 
                                                                Binning::Custom(TrueEnergyBinning),
                                                                NuLoader, 
                                                                kCC1e0p1Signal_TrueNeutrinoEnergy_MakeSelectionStep(SelectionSteps[iSel].cut), 
                                                                kCRTPMTNeutrino); 
    }

    NuLoader.Go();

    TFile FOut("CC1e0piSelection_Efficiency.root", "recreate");
    double TargetPOT;

    // selection efficiency
    TCanvas* cEffProg = new TCanvas("efficiencyCC1e0piselection", "efficiencyCC1e0piselection", 500, 500);
    TLegend* lEffProg = new TLegend(0.6, 0.5, 0.8, 0.85, "NuMI CV");

    TargetPOT = sTrueNeutrinoEnergy->POT();                   
    TH1* hTrue = sTrueNeutrinoEnergy->ToTH1(TargetPOT);
    hTrue->SetFillColorAlpha(kGray, 0.5);
    hTrue->SetLineColor(0); 
    hTrue->SetTitle(";E_{#nu} [GeV];Selection efficiency");
    TH1* hTrue_Scaled = (TH1*) hTrue->Clone();
    hTrue_Scaled->Scale(1. / hTrue->GetMaximum());
    hTrue_Scaled->Draw("HIST SAME");

    for(unsigned int iSel = 0; iSel < kNSelectionSteps; ++iSel) {
        TargetPOT = sTrueNeutrinoEnergy_SelectionSteps[iSel]->POT();
        TH1* h = sTrueNeutrinoEnergy_SelectionSteps[iSel]->ToTH1(TargetPOT);

        std::cout << h->GetEntries() << "\t" << hTrue->GetEntries() << "\t" << h->GetEntries() / hTrue->GetEntries() << std::endl;
        TString labelEffProg = SelectionSteps[iSel].label + Form(" (%.0f%%)", 100.* h->GetEntries() / hTrue->GetEntries());

        TEfficiency* eff = new TEfficiency(*h, *hTrue);
        eff->SetLineWidth(2);
        eff->SetLineColor(SelectionSteps[iSel].color);
        eff->SetMarkerColor(SelectionSteps[iSel].color);

        h->SetLineColor(SelectionSteps[iSel].color);
        lEffProg->AddEntry(h, labelEffProg, "f");

        eff->Draw("AP SAME");
    }

    lEffProg->SetTextSize(0.03);
    lEffProg->SetFillStyle(0);
    lEffProg->Draw();
    cEffProg->Write();

    FOut.Close();

    return;
}