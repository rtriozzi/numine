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
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

using namespace ana;

void CC1e0piSelection_MakeEfficiency() {

    // CNAF NuMI MC
    const std::string TargetFile = "/storage/gpfs_data/icarus/local/users/cfarnese/NUMI/NUMI_MC/*.root";

    SpectrumLoader NuLoader(TargetFile);

    // selection efficiency
    Spectrum* sTrueNeutrinoEnergy = new Spectrum("E_{#nu} [GeV]", 
                                                 Binning::Simple(20, 0, 4),
                                                 NuLoader,
                                                 kCC1e0p1Signal_TrueNeutrinoEnergy, 
                                                 kNoSpillCut);

    Spectrum* sTrueNeutrinoEnergy_CRTPMT = new Spectrum("E_{#nu} [GeV]", 
                                                 Binning::Simple(20, 0, 4),
                                                 NuLoader,
                                                 kCC1e0p1Signal_TrueNeutrinoEnergy, 
                                                 kCRTPMTNeutrino);

    const unsigned int kNSelectionSteps = SelectionSteps.size();
    Spectrum *sTrueNeutrinoEnergy_SelectionSteps[kNSelectionSteps];

    for(unsigned int iSel = 0; iSel < kNSelectionSteps; ++iSel) {

        sTrueNeutrinoEnergy_SelectionSteps[iSel] = new Spectrum("E_{#nu} [GeV]", 
                                                                Binning::Simple(20, 0, 4), 
                                                                NuLoader, 
                                                                kCC1e0p1Signal_TrueNeutrinoEnergy_MakeSelectionStep(SelectionSteps[iSel].cut), 
                                                                kCRTPMTNeutrino); 
    }

    NuLoader.Go();

    TFile FOut("CC1e0piSelection_Efficiency.root", "recreate");
    double TargetPOT;

    // true energy of generated neutrinos
    TCanvas* cTrue = new TCanvas("energyCC1e0pisignal", "energyCC1e0pisignal");
    TargetPOT = sTrueNeutrinoEnergy->POT();                   
    TH1* hTrue = sTrueNeutrinoEnergy->ToTH1(TargetPOT);

    // selection efficiency
    TCanvas* cEffProg = new TCanvas("efficiencyCC1e0piselection", "efficiencyCC1e0piselection");
    TLegend *lEffProg = new TLegend(0.55, 0.55, 0.85, 0.85);

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
        gPad->Update();
        eff->GetPaintedGraph()->GetXaxis()->SetTitle("E_{#nu} [GeV]");
        eff->GetPaintedGraph()->GetYaxis()->SetTitle("Selection efficiency");
        gPad->Update();
    }

    lEffProg->Draw();
    cEffProg->Write();

    FOut.Close();

    return;
}