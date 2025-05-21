#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

// helpers
#include "CC1e0piSelection_Cuts.h"
#include "CC1e0piSelection_TruthCuts.h"

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

void CC1e0piSelection_MakePurity() {
    
    // CNAF NuMI MC
    // const std::string TargetFile = "/storage/gpfs_data/icarus/local/users/cfarnese/NUMI/NUMI_MC/*.root"; ///< CV
    const std::string TargetFile = "/storage/gpfs_data/icarus/plain/user/cfarnese/RT_production_NuMIreference_20May25/mc*/caf_here/*.flat.caf.root"; ///< reference
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/user/cfarnese/RT_NUMI_nuonly_May18/run*/cafmakerjob_here_2d_updated/*.flat.caf.root"; ///< new BDT

    SpectrumLoader NuLoader(TargetFile);

    // selection purity at each selection step
    const unsigned int kNVar = SelectionPlots.size();
    const unsigned int kNSel = SelectionSteps.size();
    Spectrum *spectra_Signal[kNVar];
    Spectrum *spectra_SelectedAll[kNVar][kNSel];
    Spectrum *spectra_SelectedSignal[kNVar][kNSel];

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        spectra_Signal[iVar] = new Spectrum(SelectionPlots[iVar].label, 
                                            SelectionPlots[iVar].bins, 
                                            NuLoader, 
                                            SelectionPlots[iVar].var, 
                                            kNoSpillCut,
                                            kTrueCC1e0pi);  

        for(unsigned int jSel = 0; jSel < kNSel; ++jSel) {
            spectra_SelectedAll[iVar][jSel] = new Spectrum(SelectionPlots[iVar].label, 
                                                           SelectionPlots[iVar].bins, 
                                                           NuLoader, 
                                                           SelectionPlots[iVar].var, 
                                                           kCRTPMTNeutrino,
                                                           SelectionSteps[jSel].cut);     

            spectra_SelectedSignal[iVar][jSel] = new Spectrum(SelectionPlots[iVar].label, 
                                                              SelectionPlots[iVar].bins, 
                                                              NuLoader, 
                                                              SelectionPlots[iVar].var, 
                                                              kCRTPMTNeutrino,
                                                              SelectionSteps[jSel].cut && kTrueCC1e0pi);    
        }
    }

    NuLoader.Go();
 
    TFile FOut("CC1e0piSelection_Purity.root", "recreate");

    TCanvas *c[kNVar];
    TLegend *l[kNVar];
    double TargetPOT;
    std::string title;

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {

        c[iVar] = new TCanvas(SelectionPlots[iVar].suffix.c_str(), SelectionPlots[iVar].suffix.c_str(), 400, 400);
        l[iVar] = new TLegend(0.6, 0.5, 0.8, 0.85, "NuMI CV");

        TargetPOT = spectra_Signal[iVar]->POT();
        TH1* hSignal = spectra_Signal[iVar]->ToTH1(TargetPOT);
        hSignal->SetFillColorAlpha(kGray, 0.5);
        hSignal->SetLineColor(0); 
        title = std::string(";") + SelectionPlots[iVar].label + std::string(";Selection purity");
        hSignal->SetTitle(title.c_str());
        TH1* hSignal_Scaled = (TH1*) hSignal->Clone();
        hSignal_Scaled->Scale(1. / hSignal->GetMaximum());
        hSignal_Scaled->Draw("HIST SAME");

        for(unsigned int jSel = 0; jSel < kNSel; ++jSel) {
            // all the selected reconstructed interactions
            TargetPOT = spectra_SelectedAll[iVar][jSel]->POT();
            TH1* hAll = spectra_SelectedAll[iVar][jSel]->ToTH1(TargetPOT);

            // selected reconstructed interactions matching signal
            TargetPOT = spectra_SelectedSignal[iVar][jSel]->POT();
            TH1* hSelected = spectra_SelectedSignal[iVar][jSel]->ToTH1(TargetPOT);

            TString labelEffProg = SelectionSteps[jSel].label + Form(" (%.0f%%)", 100.* hSelected->GetEntries() / hAll->GetEntries());

            title = std::string(";") + SelectionPlots[iVar].label + std::string(";Purity");
            TEfficiency* eff = new TEfficiency("eff", title.c_str(), SelectionPlots[iVar].bins.NBins(), &SelectionPlots[iVar].bins.Edges()[0]);

            for (int i = 1; i <= SelectionPlots[iVar].bins.NBins(); i++) {
                if (hAll->GetBinContent(i) == 0) continue;
                eff->SetTotalEvents(i, hAll->GetBinContent(i));  
                eff->SetPassedEvents(i, hSelected->GetBinContent(i));
            }

            eff->SetLineWidth(2);
            eff->SetLineColor(SelectionSteps[jSel].color);
            eff->SetMarkerColor(SelectionSteps[jSel].color);

            l[iVar]->AddEntry(eff, labelEffProg, "f");
            eff->Draw("AP SAME");
        }        

        l[iVar]->SetFillStyle(0);
        l[iVar]->SetTextSize(0.04);
        l[iVar]->Draw();
        c[iVar]->Write();
        //title = std::string("plots/") + SelectionPlots[iVar].suffix + std::string(".pdf");
        //c[iVar]->SaveAs(title.c_str());
    }

    FOut.Close();

    return;
}