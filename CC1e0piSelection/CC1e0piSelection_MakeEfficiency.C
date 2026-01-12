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

void CC1e0piSelection_MakeEfficiency() {

    // CNAF NuMI MC
    // const std::string TargetFile = "/storage/gpfs_data/icarus/local/users/cfarnese/NUMI/NUMI_MC/*.root"; ///< CV
    
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/user/cfarnese/RT_production_NuMIreference_20May25/mc*/caf_here/*.flat.caf.root"; ///< reference
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/user/cfarnese/RT_NUMI_nuonly_May18/run*/cafmakerjob_here_2d_updated/*.flat.caf.root"; ///< new BDT
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/user/cfarnese/RT_production_NuMIcheating_20May25/mc*/caf_here/*.flat.caf.root"; ///< vertex cheated
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/user/cfarnese/RT_production_NuMIcheatingnew_21May25/mc*/caf_here/*.flat.caf.root"; ///< BDT vertex closest to truth
    // const std::string TargetFile = "/storage/gpfs_data/icarus/local/users/rtriozzi/concats/NuMI_CV_MopUp_NewBDT/*.root";
    
    // FNAL NuMI MC
    // const std::string TargetFile = "/exp/icarus/data/users/rtriozzi/mc/numi_FRFIX/concat_NuMI_MC_FRFIX_*.root";

    // FNAL development MC
    // const std::string TargetFile = "/pnfs/icarus/scratch/users/rtriozzi/NuGraphCAFs_NuE/StandardReco/*.root"; ///< Standard Pandora
    // const std::string TargetFile = "/pnfs/icarus/scratch/users/rtriozzi/NuGraphCAFs_NuE/NuGraphReco/*.root"; ///< Pandora + NuGraph

    // FNAL development NuMI MC
    // const std::string TargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/standard/numinue.flat.caf.root"; ///< NuE
    // const std::string TargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/standard/numinom.flat.caf.root"; ///< nominal flux, mostly NuMu

    // FNAL development NuMI MC / standard reconstruction
    // const std::string TargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/standard/numinue.flat.caf.root"; ///< NuE
    
    // FNAL development NuMI MC / NG2 filter + NG2 PID
    const std::string TargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco/numinue_NuGraphReco.flat.caf.root"; ///< NuE

    SpectrumLoader NuLoader(TargetFile);

    // selection efficiency
    std::vector<double> TrueEnergyBinning = {0., 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2., 2.4, 2.8, 3.2, 4};

    Spectrum* sTrueNeutrinoEnergy = new Spectrum("E_{#nu} [GeV]", 
                                                 Binning::Custom(TrueEnergyBinning),
                                                 NuLoader,
                                                 kCC1e0p1Signal_TrueNeutrinoEnergy, 
                                                 kNoSpillCut);

    const unsigned int kNSelectionSteps = SelectionSteps_NoTrigger.size();
    Spectrum *sTrueNeutrinoEnergy_SelectionSteps[kNSelectionSteps];

    for(unsigned int iSel = 0; iSel < kNSelectionSteps; ++iSel) {
        sTrueNeutrinoEnergy_SelectionSteps[iSel] = new Spectrum("E_{#nu} [GeV]", 
                                                                Binning::Custom(TrueEnergyBinning),
                                                                NuLoader, 
                                                                kCC1e0p1Signal_TrueNeutrinoEnergy_MakeSelectionStep(SelectionSteps_NoTrigger[iSel].cut),
                                                                kNoSpillCut);//kCRTPMTNeutrino); 
    }

    NuLoader.Go();

    TFile FOut("CC1e0piSelection_Efficiency.root", "recreate");
    double TargetPOT;

    // selection efficiency
    TCanvas* cEffProg = new TCanvas("efficiencyCC1e0piselection", "efficiencyCC1e0piselection", 500, 500);
    TLegend* lEffProg = new TLegend(0.125, 0.7, 0.9, 0.875, "NuMI CV");
    lEffProg->SetNColumns(3);

    TargetPOT = sTrueNeutrinoEnergy->POT();                   
    TH1* hTrue = sTrueNeutrinoEnergy->ToTH1(TargetPOT);
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
        TargetPOT = sTrueNeutrinoEnergy_SelectionSteps[iSel]->POT();
        TH1* h = sTrueNeutrinoEnergy_SelectionSteps[iSel]->ToTH1(TargetPOT);
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
