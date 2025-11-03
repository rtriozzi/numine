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

void CC1e0piSelection_MakeEfficiency_MultiSample() {

    // FNAL development NuMI MC / standard reconstruction
    // const std::string TargetFile_NuE = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/standard/numinue.flat.caf.root"; ///< NuE
    // const std::string TargetFile_Nom = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/standard/numinom.flat.caf.root"; ///< nominal flux, mostly NuMu 

    // FNAL development NuMI MC / NG2 filter + NG2 PID
    const std::string TargetFile_NuE = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco/numinue_NuGraphReco.flat.caf.root"; ///< NuE
    const std::string TargetFile_Nom = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco/numinom_NuGraphReco.flat.caf.root"; ///< nominal flux, mostly NuMu 

    SpectrumLoader NuLoader_NuE(TargetFile_NuE);
    SpectrumLoader NuLoader_Nom(TargetFile_Nom);

    // selection efficiency
    std::vector<double> TrueEnergyBinning = {0., 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2., 2.4, 2.8, 3.2, 4};

    Spectrum* sTrueNeutrinoEnergy = new Spectrum("E_{#nu} [GeV]", 
                                                 Binning::Custom(TrueEnergyBinning),
                                                 NuLoader_NuE,
                                                 kCC1e0p1Signal_TrueNeutrinoEnergy, 
                                                 kNoSpillCut);
    Spectrum* sTrueNeutrinoEnergy_Nom = new Spectrum("E_{#nu} [GeV]", 
                                                 Binning::Custom(TrueEnergyBinning),
                                                 NuLoader_Nom,
                                                 kCC1e0p1Signal_TrueNeutrinoEnergy, 
                                                 kNoSpillCut);

    const unsigned int kNSelectionSteps = SelectionSteps_NoTrigger.size();
    Spectrum *sTrueNeutrinoEnergy_SelectionSteps[kNSelectionSteps];
    Spectrum *sTrueNeutrinoEnergy_SelectionSteps_Nom[kNSelectionSteps];

    for(unsigned int iSel = 0; iSel < kNSelectionSteps; ++iSel) {
        sTrueNeutrinoEnergy_SelectionSteps[iSel] = new Spectrum("E_{#nu} [GeV]", 
                                                                Binning::Custom(TrueEnergyBinning),
                                                                NuLoader_NuE, 
                                                                kCC1e0p1Signal_TrueNeutrinoEnergy_MakeSelectionStep(SelectionSteps_NoTrigger[iSel].cut),
                                                                kNoSpillCut);//kCRTPMTNeutrino); 
        sTrueNeutrinoEnergy_SelectionSteps_Nom[iSel] = new Spectrum("E_{#nu} [GeV]", 
                                                                Binning::Custom(TrueEnergyBinning),
                                                                NuLoader_Nom, 
                                                                kCC1e0p1Signal_TrueNeutrinoEnergy_MakeSelectionStep(SelectionSteps_NoTrigger[iSel].cut),
                                                                kNoSpillCut);//kCRTPMTNeutrino); 
    }

    NuLoader_NuE.Go();
    NuLoader_Nom.Go();

    TFile FOut("CC1e0piSelection_Efficiency.root", "recreate");
    // double TargetPOT = std::min(sTrueNeutrinoEnergy->POT(), sTrueNeutrinoEnergy_Nom->POT());

    // selection efficiency
    TCanvas* cEffProg = new TCanvas("efficiencyCC1e0piselection", "efficiencyCC1e0piselection", 300, 300);
    TLegend* lEffProg = new TLegend(0.125, 0.7, 0.9, 0.875, "NuMI CV");
    lEffProg->SetNColumns(3);

    TH1* hTrue = sTrueNeutrinoEnergy->ToTH1(sTrueNeutrinoEnergy->POT());
    TH1* hTrue_Nom = sTrueNeutrinoEnergy_Nom->ToTH1(sTrueNeutrinoEnergy_Nom->POT());
    hTrue->Add(hTrue_Nom);
    hTrue->SetFillColorAlpha(kGray, 0.5);
    hTrue->SetLineColor(0); 
    hTrue->SetTitle(";E_{#nu} [GeV];Selection efficiency");
    TH1* hTrue_Scaled = (TH1*) hTrue->Clone();
    hTrue_Scaled->Scale(1. / hTrue->GetMaximum());
    hTrue_Scaled->GetYaxis()->SetRangeUser(0, 1.4);
    hTrue_Scaled->Draw("HIST SAME");

    TLine *line = new TLine(0, 1, 4, 1);
    line->Draw("SAME");

    for(unsigned int iSel = 0; iSel < kNSelectionSteps; ++iSel) {
        TH1* h = sTrueNeutrinoEnergy_SelectionSteps[iSel]->ToTH1(sTrueNeutrinoEnergy_SelectionSteps[iSel]->POT());
        TH1* h_Nom = sTrueNeutrinoEnergy_SelectionSteps_Nom[iSel]->ToTH1(sTrueNeutrinoEnergy_SelectionSteps_Nom[iSel]->POT());
        h->Add(h_Nom);
        h->SetLineColor(SelectionSteps_NoTrigger[iSel].color);

        TEfficiency* eff = new TEfficiency(*h, *hTrue);
        eff->SetLineWidth(2);
        eff->SetLineColor(SelectionSteps_NoTrigger[iSel].color);
        eff->SetMarkerColor(SelectionSteps_NoTrigger[iSel].color);

        std::cout << h->GetEntries() << "\t" << hTrue->GetEntries() << "\t" << h->GetEntries() / hTrue->GetEntries() << std::endl;
        TString labelEffProg = SelectionSteps_NoTrigger[iSel].label + Form(" (%.0f%%)", 100.* h->GetEntries() / hTrue->GetEntries());
        lEffProg->AddEntry(h, labelEffProg, "f");

        eff->Draw("P SAME");
        gPad->Update();
    }

    lEffProg->SetTextSize(0.0275);
    lEffProg->SetFillStyle(0);
    lEffProg->Draw();
    cEffProg->Write();

    gStyle->SetLineScalePS(5);
    cEffProg->SaveAs("plots/CC1e0piSelection_Efficiency.pdf");

    FOut.Close();

    return;
}
