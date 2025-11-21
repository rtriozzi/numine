#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "TruthCuts.h"

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
#include "TGaxis.h"

using namespace ana;

void Michels() {

    // FNAL development NuMI MC / NG2 filter + NG2 PID
    const std::string TargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco/*nom*_NuGraphReco.flat.caf.root"; ///< nominal flux, mostly NuMu 

    SpectrumLoader NuLoader(TargetFile);

    const unsigned int kNVar = SelectionPlots.size();
    const unsigned int kNSel = ParticleTypes.size();
    Spectrum *spectra[kNVar][kNSel];

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        for(unsigned int jSel = 0; jSel < kNSel; ++jSel){
            spectra[iVar][jSel] = new Spectrum(SelectionPlots[iVar].label, 
                                               SelectionPlots[iVar].bins, 
                                               NuLoader, 
                                               SelectionPlots[iVar].var, 
                                               kNoSpillCut,
                                               kMichelSelection && ParticleTypes[jSel].cut);  
        }
    }

    NuLoader.Go();

    TFile FOut("Michels.root", "recreate");

    TCanvas *c[kNVar];
    TLegend *l[kNVar];
    THStack *hs[kNVar];
    double TargetPOT = 0.;
    std::string title;

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        gStyle->SetCanvasDefW(250); gStyle->SetCanvasDefH(250); 
        c[iVar] = new TCanvas(SelectionPlots[iVar].suffix.c_str(), SelectionPlots[iVar].suffix.c_str(), 300, 300);
        c[iVar]->SetTopMargin(0.05); c[iVar]->SetRightMargin(0.05); c[iVar]->SetBottomMargin(0.175); c[iVar]->SetLeftMargin(0.175);
        hs[iVar] = new THStack(SelectionPlots[iVar].suffix.c_str(), SelectionPlots[iVar].label.c_str());
        l[iVar] = new TLegend(0.625, 0.625, 0.85, 0.925, "NuMI CV");

        TGaxis::SetMaxDigits(3);

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

            h->SetFillColor(ParticleTypes[jSel].color);
            h->SetLineWidth(0);
            h->SetLineColor(kBlack);
            l[iVar]->AddEntry(h, ParticleTypes[jSel].label.c_str(), "f");

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
            box->SetFillColor(ParticleTypes[0].color);
            box->SetLineColor(ParticleTypes[0].color);
            box->Draw("SAME");
        }
        hAll->SetLineColor(ParticleTypes[0].color);
        hAll->SetLineWidth(2);

        l[iVar]->SetTextSize(0.06);
        l[iVar]->Draw();
        c[iVar]->Write();

        gStyle->SetLabelSize(0.07, "XY"); gStyle->SetTitleSize(0.08, "XY");
        gStyle->SetTitleOffset(1.1, "Y"); gStyle->SetTitleOffset(1.1, "X");

        gStyle->SetLineScalePS(5);
        title = std::string("plots/") + SelectionPlots[iVar].suffix + std::string(".pdf");
        c[iVar]->SaveAs(title.c_str());
    }

    FOut.Close();

    return;
}
