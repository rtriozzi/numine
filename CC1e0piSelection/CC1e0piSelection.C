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

    const unsigned int kNVar = SelectionPlots.size();
    const unsigned int kNSel = InteractionTypes.size();
    Spectrum *spectra[kNVar][kNSel];

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        for(unsigned int jSel = 0; jSel < kNSel; ++jSel){
            spectra[iVar][jSel] = new Spectrum(SelectionPlots[iVar].label, 
                                               SelectionPlots[iVar].bins, 
                                               NuLoader, 
                                               SelectionPlots[iVar].var, 
                                               kCRTPMTNeutrino,
                                               InteractionTypes[jSel].cut);          
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
        c[iVar] = new TCanvas(SelectionPlots[iVar].suffix.c_str(), SelectionPlots[iVar].suffix.c_str(), 400, 400);
        hs[iVar] = new THStack(SelectionPlots[iVar].suffix.c_str(), SelectionPlots[iVar].label.c_str());
        l[iVar] = new TLegend(0.65, 0.5, 0.85, 0.85, "NuMI CV");

        // all slices with margins
        TargetPOT = spectra[iVar][0]->POT();
        TH1* hAll = spectra[iVar][0]->ToTH1(TargetPOT);

        // // handle underflow and overflow
        // float xMin = hAll->GetXaxis()->GetXmin(), xMax = hAll->GetXaxis()->GetXmax();
        // float binWidth = (xMax - xMin) / hAll->GetNbinsX();
        // TH1F* hAll_WithUF = new TH1F("", "", hAll->GetNbinsX()+2, xMin-binWidth, xMax+binWidth);
        // for (int i = 1; i <= hAll->GetNbinsX(); ++i)
        //     hAll_WithUF->SetBinContent(i+1, hAll->GetBinContent(i));
        // hAll_WithUF->SetBinContent(1, hAll->GetBinContent(0));

        // int yMax = 0;
        // for (int i = 1; i <= hAll_WithUF->GetNbinsX(); ++i) {
        //     double y = hAll_WithUF->GetBinContent(i);
        //     double err = hAll_WithUF->GetBinError(i);
        //     if (y + err > yMax) yMax = y + err;
        // }

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

            // // handle underflow and overflow
            // float xMin = h->GetXaxis()->GetXmin(), xMax = h->GetXaxis()->GetXmax();
            // float binWidth = (xMax - xMin) / h->GetNbinsX();
            // TH1F* hWithUF = new TH1F("", "", h->GetNbinsX()+2, xMin-binWidth, xMax+binWidth);
            // for (int i = 1; i <= h->GetNbinsX(); ++i)
            //     hWithUF->SetBinContent(i+1, h->GetBinContent(i));
            // hWithUF->SetBinContent(1, h->GetBinContent(0));
            // hWithUF->SetFillColor(InteractionTypes[jSel].color);
            // hWithUF->SetLineWidth(0);
            // hWithUF->SetLineColor(kBlack);
            // l[iVar]->AddEntry(hWithUF, InteractionTypes[jSel].label.c_str(), "f");

            // hs[iVar]->Add(hWithUF);

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

        // // all slices with errors
        // gStyle->SetHatchesLineWidth(2);
        // gStyle->SetHatchesSpacing(0.5);
        // for (int i = 1; i <= hAll_WithUF->GetNbinsX(); ++i) {
        //     double xlow = hAll_WithUF->GetBinLowEdge(i);
        //     double xup = xlow + hAll_WithUF->GetBinWidth(i);
        //     double y = hAll_WithUF->GetBinContent(i);
        //     double err = hAll_WithUF->GetBinError(i);

        //     TBox* box = new TBox(xlow, y - err, xup, y + err);
        //     box->SetFillStyle(3004); 
        //     box->SetFillColor(InteractionTypes[0].color);
        //     box->SetLineColor(InteractionTypes[0].color);
        //     box->Draw("SAME");
        // }
        // hAll_WithUF->SetLineColor(InteractionTypes[0].color);
        // hAll_WithUF->SetLineWidth(2);
        // hAll_WithUF->Draw("HIST SAME");

        // all slices with errors
        gStyle->SetHatchesLineWidth(2);
        gStyle->SetHatchesSpacing(0.5);
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
        // title = std::string("plots/") + SelectionPlots[iVar].suffix + std::string(".pdf");
        // c[iVar]->SaveAs(title.c_str());
    }

    FOut.Close();

    return;
}