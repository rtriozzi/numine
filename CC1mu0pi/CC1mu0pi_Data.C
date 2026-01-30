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

void CC1mu0pi_Data() {

    // FNAL development NuMI prescaled data / NG2 filter + NG2 PID
    const std::string DataTargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_HIPTagger/numi_prescaled_NuGraphReco_HIPTagger.unblind.flat.caf.root";

    SpectrumLoader dataNuLoader(DataTargetFile);

    const unsigned int kNVar = NuMuSelectionPlots.size();
    Spectrum *dataSpectra[kNVar];

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
            dataSpectra[iVar] = new Spectrum(NuMuSelectionPlots[iVar].label, 
                                             NuMuSelectionPlots[iVar].bins, 
                                             dataNuLoader, 
                                             NuMuSelectionPlots[iVar].var, 
                                             kNoSpillCut, //kCRTPMTNeutrino,
                                             kAutomaticNuMuSelection);          
    }

    dataNuLoader.Go();

    // FNAL development NuMI MC / NG2 filter + NG2 PID
    const std::string MCTargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_MoreVars/numinom*.flat.caf.root"; ///< nominal flux, mostly NuMu 

    SpectrumLoader NuLoader(MCTargetFile);

    const unsigned int kNSel = InteractionTypes.size();
    Spectrum *spectra[kNVar][kNSel];

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        for(unsigned int jSel = 0; jSel < kNSel; ++jSel){
            spectra[iVar][jSel] = new Spectrum(NuMuSelectionPlots[iVar].label, 
                                               NuMuSelectionPlots[iVar].bins, 
                                               NuLoader, 
                                               NuMuSelectionPlots[iVar].var, 
                                               kNoSpillCut, // kCRTPMTNeutrino,
                                               kAutomaticNuMuSelection && InteractionTypes[jSel].cut);          
        }
    }

    // Spectrum *sEventDump = new Spectrum("", Binning::Simple(3, 0, 3), NuLoader, kDebugger, kNoSpillCut); 

    NuLoader.Go();

    TFile FOut("CC1mu0pi_Data.root", "recreate");

    TCanvas *c[kNVar];
    TLegend *l[kNVar];
    THStack *hs[kNVar];
    std::string title;
    double TargetMCPOT = spectra[0][0]->POT();

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        gStyle->SetCanvasDefW(250); gStyle->SetCanvasDefH(250); 
        c[iVar] = new TCanvas(NuMuSelectionPlots[iVar].suffix.c_str(), NuMuSelectionPlots[iVar].suffix.c_str(), 300, 300);
        c[iVar]->SetTopMargin(0.025); c[iVar]->SetRightMargin(0.025); c[iVar]->SetBottomMargin(0.225); c[iVar]->SetLeftMargin(0.225);
        hs[iVar] = new THStack(NuMuSelectionPlots[iVar].suffix.c_str(), NuMuSelectionPlots[iVar].label.c_str());
        l[iVar] = new TLegend(0.625, 0.425, 0.85, 0.925, "NuMI CV");

        // set up data plot
        TH1* hData = dataSpectra[iVar]->ToTH1(1e18);
        hData->Scale(1.0 / hData->Integral());
        float yMax = 0;
        for (int i = 1; i <= hData->GetNbinsX(); ++i) {
            double y = hData->GetBinContent(i);
            double err = hData->GetBinError(i);
            if (y + err > yMax) yMax = y + err;
        }

        // all slices with margins
        TH1* hAll = spectra[iVar][0]->ToTH1(TargetMCPOT);
        float MCIntegral = hAll->Integral();
        hAll->Scale(1.0 / MCIntegral);

        // stack by interaction type
        for(unsigned int jSel = 1; jSel < kNSel; ++jSel) {
            TH1* h = spectra[iVar][jSel]->ToTH1(TargetMCPOT);

            h->SetFillColor(InteractionTypes[jSel].color);
            h->SetFillStyle(1001);
            h->SetLineColor(h->GetFillColor());
            h->SetLineWidth(0);
            l[iVar]->AddEntry(h, InteractionTypes[jSel].label.c_str(), "f");

            h->Scale(1.0 / MCIntegral);
            hs[iVar]->Add(h);
        }

        if (yMax == 0) yMax = 0.05;
        hs[iVar]->SetMaximum(yMax + 0.1*yMax + 0.025);
        hs[iVar]->Draw("HIST");

        title = std::string(";") + 
                NuMuSelectionPlots[iVar].label + std::string(";") + 
                Form("Slices [a.n.]");
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
            box->SetFillStyle(3004); 
            box->SetFillColor(InteractionTypes[0].color);
            box->SetLineColor(InteractionTypes[0].color);
            box->Draw("SAME");
        }
        hAll->SetLineColor(InteractionTypes[0].color);
        hAll->SetLineWidth(2);

        // plot data
        hData->SetMarkerStyle(20); 
        hData->SetLineWidth(2);
        hData->SetMarkerSize(0.5);
        hData->SetMarkerColor(kBlack);
        hData->Draw("EX0 SAME");
        l[iVar]->AddEntry(hData, "Data", "pe");

        l[iVar]->SetTextSize(0.06);
        l[iVar]->Draw();
        c[iVar]->Write();

        gStyle->SetLabelSize(0.07, "XY"); gStyle->SetTitleSize(0.08, "XY");
        gStyle->SetTitleOffset(1.25, "Y"); gStyle->SetTitleOffset(1.25, "X");

        gStyle->SetLineScalePS(5);
        title = std::string("plots/") + NuMuSelectionPlots[iVar].suffix + std::string(".pdf");
        c[iVar]->SaveAs(title.c_str());
    }

    FOut.Close();

    return;
}
