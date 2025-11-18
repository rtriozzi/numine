#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "CC1e0piSelection_Cuts.h"
#include "CC1e0piSelection_TruthCuts.h"
#include "CC1e0piSelection_Efficiency.h"
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

void CC1e0piSelection_MultiSample() {

    // FNAL development NuMI MC / standard reconstruction
    // const std::string TargetFile_NuE = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/standard/numinue.flat.caf.root"; ///< NuE
    // const std::string TargetFile_Nom = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/standard/*nom*.flat.caf.root"; ///< nominal flux, mostly NuMu, new production :) 
    // const std::string TargetFile_Nom = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/standard/numinom.flat.caf.root"; ///< nominal flux, mostly NuMu 

    // FNAL development NuMI MC / NG2 filter + NG2 PID
    const std::string TargetFile_NuE = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco/numinue_NuGraphReco.flat.caf.root"; ///< NuE
    const std::string TargetFile_Nom = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco/*nom*_NuGraphReco.flat.caf.root"; ///< nominal flux, mostly NuMu 

    // FNAL development NuMI MC / NG2 filter + NG2 PID / New reco
    // const std::string TargetFile_NuE = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_newvtx/numinue_NuGraphReco_NewVtx.flat.caf.root"; ///< NuE
    // const std::string TargetFile_Nom = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_newvtx/*nom*_NuGraphReco_NewVtx.flat.caf.root"; ///< nominal flux, mostly NuMu 

    // FNAL development NuMI MC / NG2 filter + NG2 PID / Semantic splitting
    // const std::string TargetFile_NuE = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_semsplit/numinue_NuGraphReco_OnlySemSplit.flat.caf.root"; ///< NuE
    // const std::string TargetFile_Nom = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_semsplit/*nom*_NuGraphReco_OnlySemSplit.flat.caf.root"; ///< nominal flux, mostly NuMu 

    // FNAL development NuMI MC / NG2 filter + NG2 PID / Only vertex promotion
    // const std::string TargetFile_NuE = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_onlyvtxprom/numinue_NuGraphReco_OnlyVtxProm.flat.caf.root"; ///< NuE
    // const std::string TargetFile_Nom = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_onlyvtxprom/*nom*_NuGraphReco_OnlyVtxProm.flat.caf.root"; ///< nominal flux, mostly NuMu 

    // FNAL development NuMI MC / NG2 filter + NG2 PID / Only vertex promotion + semantic confidence cut at 3
    // const std::string TargetFile_NuE = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_onlyvtxprom_semconf3/numinue_NuGraphReco_OnlyVtxProm_SemConf3.flat.caf.root"; ///< NuE
    // const std::string TargetFile_Nom = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_onlyvtxprom_semconf3/*nom*_NuGraphReco_OnlyVtxProm_SemConf3.flat.caf.root"; ///< nominal flux, mostly NuMu 

    // FNAL development NuMI MC / NG2 filter + NG2 PID / NG2-based shower growing
    // const std::string TargetFile_NuE = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_ng2showergrowing/numinue_NuGraphReco_NG2ShowerGrowing.flat.caf.root"; ///< NuE
    // const std::string TargetFile_Nom = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_ng2showergrowing/numinom_noyzsim_NuGraphReco_NG2ShowerGrowing.flat.caf.root"; ///< nominal flux, mostly NuMu 

    SpectrumLoader NuLoader_NuE(TargetFile_NuE);
    SpectrumLoader NuLoader_Nom(TargetFile_Nom);

    Spectrum *sMCEventDump = new Spectrum("", Binning::Simple(3, 0, 3), NuLoader_NuE, kMCEventDump_DebugSelection, kNoSpillCut);

    const unsigned int kNVar = SelectionPlots.size();
    const unsigned int kNSel = InteractionTypes.size();
    Spectrum *spectra_NuE[kNVar][kNSel];
    Spectrum *spectra_Nom[kNVar][kNSel];

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        for(unsigned int jSel = 0; jSel < kNSel; ++jSel){
            spectra_NuE[iVar][jSel] = new Spectrum(SelectionPlots[iVar].label, 
                                                SelectionPlots[iVar].bins, 
                                                NuLoader_NuE, 
                                                SelectionPlots[iVar].var, 
                                                kNoSpillCut, // kCRTPMTNeutrino, //kCRTPMTNeutrino && kNoPileUp,
                                                kAutomaticSelection && InteractionTypes[jSel].cut);  ///< change selection here if needed  
            spectra_Nom[iVar][jSel] = new Spectrum(SelectionPlots[iVar].label, 
                                                SelectionPlots[iVar].bins, 
                                                NuLoader_Nom, 
                                                SelectionPlots[iVar].var, 
                                                kNoSpillCut, // kCRTPMTNeutrino, //kCRTPMTNeutrino && kNoPileUp,
                                                kIsNotNue && kAutomaticSelection && InteractionTypes[jSel].cut);  ///< change selection here if needed                                                  
        }
    }

    NuLoader_NuE.Go();
    NuLoader_Nom.Go();

    TFile FOut("CC1e0piSelection.root", "recreate");

    TCanvas *c[kNVar];
    TLegend *l[kNVar];
    THStack *hs[kNVar];
    std::string title;
    double TargetPOT = std::min(spectra_NuE[0][0]->POT(), spectra_Nom[0][0]->POT());

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        gStyle->SetCanvasDefW(250); gStyle->SetCanvasDefH(250); 
        c[iVar] = new TCanvas(SelectionPlots[iVar].suffix.c_str(), SelectionPlots[iVar].suffix.c_str(), 300, 300);
        c[iVar]->SetTopMargin(0.025); c[iVar]->SetRightMargin(0.025); c[iVar]->SetBottomMargin(0.15); c[iVar]->SetLeftMargin(0.15);
        hs[iVar] = new THStack(SelectionPlots[iVar].suffix.c_str(), SelectionPlots[iVar].label.c_str());
        l[iVar] = new TLegend(0.65, 0.45, 0.875, 0.925, "NuMI CV");

        // all slices with margins
        TH1* hAll = spectra_NuE[iVar][0]->ToTH1(TargetPOT);
        TH1* hAll_Nom = spectra_Nom[iVar][0]->ToTH1(TargetPOT);
        hAll->Add(hAll_Nom);

        int yMax = 0;
        for (int i = 1; i <= hAll->GetNbinsX(); ++i) {
            double y = hAll->GetBinContent(i);
            double err = hAll->GetBinError(i);
            if (y + err > yMax) yMax = y + err;
        }
        
        // stack by interaction type
        for(unsigned int jSel = 1; jSel < kNSel; ++jSel) {
            TH1* h = spectra_NuE[iVar][jSel]->ToTH1(TargetPOT);
            TH1* h_Nom = spectra_Nom[iVar][jSel]->ToTH1(TargetPOT);
            h->Add(h_Nom);

            h->SetFillColor(InteractionTypes[jSel].color);
            h->SetFillStyle(1001);
            h->SetLineColor(h->GetFillColor());
            h->SetLineWidth(0);
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

        l[iVar]->SetTextSize(0.06);
        l[iVar]->Draw();
        c[iVar]->Write();

        gStyle->SetLabelSize(0.06, "XY"); gStyle->SetTitleSize(0.07, "XY");
        gStyle->SetTitleOffset(0.9, "Y"); gStyle->SetTitleOffset(0.9, "X");

        gStyle->SetLineScalePS(5);
        title = std::string("plots/") + SelectionPlots[iVar].suffix + std::string(".pdf");
        c[iVar]->SaveAs(title.c_str());
    }

    FOut.Close();

    return;
}

