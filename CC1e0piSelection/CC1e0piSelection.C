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

// c++ stuff
#include <vector>
#include <algorithm>
#include <numeric>

using namespace ana;

void CC1e0piSelection() {

    // CNAF NuMI MC
    // const std::string TargetFile = "/storage/gpfs_data/icarus/local/users/cfarnese/NUMI/NUMI_MC/*.root"; ///< CV
    
    const std::string TargetFile = "/storage/gpfs_data/icarus/plain/user/cfarnese/RT_production_NuMIreference_20May25/mc*/caf_here/*.flat.caf.root"; ///< reference
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/user/cfarnese/RT_NUMI_nuonly_May18/run*/cafmakerjob_here_2d_updated/*.flat.caf.root"; ///< new BDT
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/user/cfarnese/RT_production_NuMIcheating_20May25/mc*/caf_here/*.flat.caf.root"; ///< vertex cheated
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/user/cfarnese/RT_production_NuMIcheatingnew_21May25/mc*/caf_here/*.flat.caf.root"; ///< BDT vertex closest to truth

    SpectrumLoader NuLoader(TargetFile);

    const unsigned int kNVar = SelectionPlots.size();
    const unsigned int kNSel = InteractionTypes.size();
    Spectrum *spectra[kNVar][kNSel];

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        for(unsigned int jSel = 0; jSel < kNSel; ++jSel){
            spectra[iVar][jSel] = new Spectrum(SelectionPlots[iVar].label, 
                                               SelectionPlots[iVar].bins, 
                                               NuLoader, 
                                               SelectionPlots[iVar].var, 
                                               kCRTPMTNeutrino, //kCRTPMTNeutrino && kNoPileUp
                                               kAutomaticSelection_NoTrigger && InteractionTypes[jSel].cut);  ///< change selection here if needed    
        }
    }

    NuLoader.Go();

    TFile FOut("CC1e0piSelection.root", "recreate");

    TCanvas *c[kNVar];
    TLegend *l[kNVar];
    THStack *hs[kNVar];
    double TargetPOT = 0.;
    std::string title;

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        c[iVar] = new TCanvas(SelectionPlots[iVar].suffix.c_str(), SelectionPlots[iVar].suffix.c_str(), 500, 500);
        hs[iVar] = new THStack(SelectionPlots[iVar].suffix.c_str(), SelectionPlots[iVar].label.c_str());
        l[iVar] = new TLegend(0.65, 0.5, 0.85, 0.85, "NuMI CV");

        // all slices with margins
        TargetPOT = spectra[iVar][0]->POT();
        TH1* hAll = spectra[iVar][0]->ToTH1(TargetPOT);

        int yMax = 0;
        for (int i = 1; i <= hAll->GetNbinsX(); ++i) {
            double y = hAll->GetBinContent(i);
            double err = hAll->GetBinError(i);
            if (y + err > yMax) yMax = y + err;
        }
        
        // stack by interaction type
        for(unsigned int jSel = 1; jSel < kNSel; ++jSel) {
            TargetPOT = spectra[iVar][jSel]->POT();
            TH1* h = spectra[iVar][jSel]->ToTH1(TargetPOT);

            h->SetFillColor(InteractionTypes[jSel].color);
            h->SetLineWidth(0);
            h->SetLineColor(kBlack);
            l[iVar]->AddEntry(h, InteractionTypes[jSel].label.c_str(), "f");

            hs[iVar]->Add(h);
        }

        hs[iVar]->SetMaximum(yMax + 0.1*yMax);
        hs[iVar]->Draw("HIST");

        title = std::string(";") + 
                SelectionPlots[iVar].label + std::string(";") + 
                Form("Slices / %.1e POT", TargetPOT);
        hs[iVar]->SetTitle(title.c_str());
        gPad->Modified();
        gPad->Update();

        // errors
        gStyle->SetHatchesLineWidth(3);
        gStyle->SetHatchesSpacing(1.5);
        for (int i = 1; i <= hAll->GetNbinsX(); ++i) {
            double xlow = hAll->GetBinLowEdge(i);
            double xup = xlow + hAll->GetBinWidth(i);
            double y = hAll->GetBinContent(i);
            double err = hAll->GetBinError(i);

            TBox* box = new TBox(xlow, y - err, xup, y + err);
            box->SetFillStyle(3003); 
            box->SetFillColor(InteractionTypes[0].color);
            box->SetLineColor(InteractionTypes[0].color);
            box->Draw("SAME");
        }
        hAll->SetLineColor(InteractionTypes[0].color);
        hAll->SetLineWidth(2);
        hAll->Draw("HIST SAME");

        l[iVar]->SetTextSize(0.04);
        l[iVar]->Draw();
        c[iVar]->Write();

        gStyle->SetLineScalePS(5);
        title = std::string("plots/") + SelectionPlots[iVar].suffix + std::string(".pdf");
        c[iVar]->SaveAs(title.c_str());
    }

    FOut.Close();

    return;
}