#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "CC1e0pi0p_Cuts.h"
#include "CC1e0pi0p_TruthCuts.h"

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

void CC1e0pi0p_Data_MultiSample_Offbeam() {

    // FNAL development NuMI prescaled data / NG2 filter + NG2 PID
    const std::string DataTargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_HIPTagger/numi_prescaled_NuGraphReco_HIPTagger.unblind.flat.caf.root";

    SpectrumLoader dataNuLoader(DataTargetFile);

    const unsigned int kNVar = SelectionPlots_LowStat.size();
    Spectrum *dataSpectra[kNVar];

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
            dataSpectra[iVar] = new Spectrum(SelectionPlots_LowStat[iVar].label, 
                                             SelectionPlots_LowStat[iVar].bins, 
                                             dataNuLoader, 
                                             SelectionPlots_LowStat[iVar].var, 
                                             kNoSpillCut,//kCRTPMTNeutrino,
                                             kAutomaticSelection);          
    }
    
    dataNuLoader.Go();

    // FNAL development NuMI off-beam data / NG2 filter + NG2 PID
    const std::string OffBeamDataTargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_HIPTagger/numi_offbeam_NuGraphReco_HIPTagger.unblind.flat.caf.root";

    SpectrumLoader offBeamDataNuLoader(OffBeamDataTargetFile);

    Spectrum *offBeamDataSpectra[kNVar];

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
            offBeamDataSpectra[iVar] = new Spectrum(SelectionPlots_LowStat[iVar].label, 
                                                    SelectionPlots_LowStat[iVar].bins, 
                                                    offBeamDataNuLoader, 
                                                    SelectionPlots_LowStat[iVar].var, 
                                                    kNoSpillCut,//kCRTPMTNeutrino,
                                                    kAutomaticSelection);          
    }
    
    offBeamDataNuLoader.Go();

    // FNAL development NuMI MC / NG2 filter + NG2 PID
    const std::string TargetFile_NuE = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_HIPTagger/numinue_NuGraphReco_HIPTagger.flat.caf.root"; ///< NuE
    const std::string TargetFile_Nom = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_HIPTagger/*nom*_NuGraphReco_HIPTagger.flat.caf.root"; ///< nominal flux, mostly NuMu 

    SpectrumLoader NuLoader_NuE(TargetFile_NuE);
    SpectrumLoader NuLoader_Nom(TargetFile_Nom);

    const unsigned int kNSel = InteractionTypes.size();
    Spectrum *spectra_NuE[kNVar][kNSel];
    Spectrum *spectra_Nom[kNVar][kNSel];

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        for(unsigned int jSel = 0; jSel < kNSel; ++jSel){
            spectra_NuE[iVar][jSel] = new Spectrum(SelectionPlots_LowStat[iVar].label, 
                                                SelectionPlots_LowStat[iVar].bins, 
                                                NuLoader_NuE, 
                                                SelectionPlots_LowStat[iVar].var, 
                                                kNoSpillCut, // kCRTPMTNeutrino, //kCRTPMTNeutrino && kNoPileUp,
                                                kAutomaticSelection && InteractionTypes[jSel].cut);  ///< change selection here if needed  
            spectra_Nom[iVar][jSel] = new Spectrum(SelectionPlots_LowStat[iVar].label, 
                                                SelectionPlots_LowStat[iVar].bins, 
                                                NuLoader_Nom, 
                                                SelectionPlots_LowStat[iVar].var, 
                                                kNoSpillCut, // kCRTPMTNeutrino, //kCRTPMTNeutrino && kNoPileUp,
                                                kIsNotNue && kAutomaticSelection && InteractionTypes[jSel].cut);  ///< change selection here if needed                                                  
        }
    }

    NuLoader_NuE.Go();
    NuLoader_Nom.Go();

    TFile FOut("CC1e0piSelection_Data_Offbeam.root", "recreate");

    TCanvas *c[kNVar];
    TLegend *l[kNVar];
    THStack *hs[kNVar];
    std::string title;

    // POT handling for MC
    // just normalize multiple neutrino samples to the lowest POT value
    double TargetMCPOT = std::min(spectra_NuE[0][0]->POT(), spectra_Nom[0][0]->POT());
    std::cout << "Normalizing MC neutrinos to " << TargetMCPOT << "POT" << std::endl;

    // POT handling for off-beam data
    // TargetOffBeamPOT = Livetime * 6e13 POT gives the number of POT-equivalent OffBeam
    // NeutrinoMCToOffBeamFactor = TargetOffBeamPOT / TargetMCPOT gives the neutrino-to-OffBeam scaling 
    // the TargetOffBeamLivetime will be the Livetime / LivetimeToPOTScaleFactor
    double TargetOffBeamEquivalentPOT = offBeamDataSpectra[0]->Livetime() * 6.e13;
    std::cout << "The OffBeam livetime (" << offBeamDataSpectra[0]->Livetime()  << ") corresponds to " 
              << TargetOffBeamEquivalentPOT << " POT" << std::endl;
    double NeutrinoMCToOffBeamFactor = TargetOffBeamEquivalentPOT / TargetMCPOT;
    std::cout << "The Off-Beam equivalent-POT over MC neutrino POT ratio is " << NeutrinoMCToOffBeamFactor << std::endl;
    double TargetOffBeamLivetime = offBeamDataSpectra[0]->Livetime() / NeutrinoMCToOffBeamFactor;
    std::cout << "The target Off-Beam livetime is " << TargetOffBeamLivetime << std::endl;

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        gStyle->SetCanvasDefW(250); gStyle->SetCanvasDefH(250); 
        c[iVar] = new TCanvas(SelectionPlots_LowStat[iVar].suffix.c_str(), SelectionPlots_LowStat[iVar].suffix.c_str(), 300, 300);
        c[iVar]->SetTopMargin(0.025); c[iVar]->SetRightMargin(0.025); c[iVar]->SetBottomMargin(0.225); c[iVar]->SetLeftMargin(0.225);
        hs[iVar] = new THStack(SelectionPlots_LowStat[iVar].suffix.c_str(), SelectionPlots_LowStat[iVar].label.c_str());
        l[iVar] = new TLegend(0.625, 0.41, 0.85, 0.925, "NuMI CV");

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
        TH1* hAll = spectra_NuE[iVar][0]->ToTH1(TargetMCPOT);
        TH1* hAll_Nom = spectra_Nom[iVar][0]->ToTH1(TargetMCPOT);
        TH1* hOffBeam = offBeamDataSpectra[iVar]->ToTH1(TargetOffBeamLivetime, kLivetime);
        hAll->Add(hAll_Nom);
        hAll->Add(hOffBeam);
        float MCIntegral = hAll->Integral();
        hAll->Scale(1.0 / MCIntegral);

        // stack by interaction type
        for(unsigned int jSel = 1; jSel < kNSel; ++jSel) {
            TH1* h = spectra_NuE[iVar][jSel]->ToTH1(TargetMCPOT);
            TH1* h_Nom = spectra_Nom[iVar][jSel]->ToTH1(TargetMCPOT);
            h->Add(h_Nom);

            h->SetFillColor(InteractionTypes[jSel].color);
            h->SetFillStyle(1001);
            h->SetLineColor(h->GetFillColor());
            h->SetLineWidth(0);
            l[iVar]->AddEntry(h, InteractionTypes[jSel].label.c_str(), "f");

            h->Scale(1.0 / MCIntegral);
            hs[iVar]->Add(h);
        }

        // add off-beam to MC
        hOffBeam->SetFillColor(kAzure-3);
        hOffBeam->SetFillStyle(3005);
        hOffBeam->SetLineColor(kAzure-3);
        hOffBeam->SetLineWidth(1);
        hOffBeam->Scale(1.0 / MCIntegral);
        l[iVar]->AddEntry(hOffBeam, "Off-beam", "f");
        hs[iVar]->Add(hOffBeam);
        
        hs[iVar]->SetMaximum(yMax + 0.1*yMax + 0.1);
        hs[iVar]->Draw("HIST");

        title = std::string(";") + 
                SelectionPlots_LowStat[iVar].label + std::string(";") + 
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
        // hAll->Draw("HIST SAME");

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
        title = std::string("plots/") + SelectionPlots_LowStat[iVar].suffix + std::string(".pdf");
        c[iVar]->SaveAs(title.c_str());
    }

    FOut.Close();

    return;
}

