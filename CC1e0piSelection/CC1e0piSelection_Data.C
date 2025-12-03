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

void CC1e0piSelection_Data() {

    // CNAF Prescaled data
    // const std::string DataTargetFile = "/storage/gpfs_data/icarus/local/users/rtriozzi/concats/NuMI_Prescaled/NuMI_Prescaled_1.root";
    // FNAL Prescaled data
    const std::string DataTargetFile = "/pnfs/sbn/data/sbn_fd/poms_production/data/Run2reprocess/reconstructed/icaruscode_v09_89_01_02p02/numimajority/flatcaf_prescaled/*/*/*.root";

    SpectrumLoader dataNuLoader(DataTargetFile);

    const unsigned int kNVar = SelectionPlots.size();
    Spectrum *dataSpectra[kNVar];

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
            dataSpectra[iVar] = new Spectrum(SelectionPlots[iVar].label, 
                                             SelectionPlots[iVar].bins, 
                                             dataNuLoader, 
                                             SelectionPlots[iVar].var, 
                                             kCRTPMTNeutrino,
                                             kPreSelection);          
    }
    
    // Spectrum *sEventDump = new Spectrum("", Binning::Simple(3, 0, 3), dataNuLoader, kEventDump, kCRTPMTNeutrino); 

    dataNuLoader.Go();

    // CNAF NuMI MC
    // const std::string MCTargetFile = "/storage/gpfs_data/icarus/local/users/cfarnese/NUMI/NUMI_MC/0.root";
    // FNAL NuMI MC
    // const std::string MCTargetFile = "/exp/icarus/data/users/rtriozzi/mc/numi_FRFIX/concat_NuMI_MC_FRFIX_*.root";
    const std::string MCTargetFile = "/pnfs/sbn/data/sbn_fd/poms_production/2025A_icarus_NuMI_MC/FHC_NuMI/mc/reconstructed/icaruscode_v09_89_01_02p02/flatcaf/*/*/*.root";

    SpectrumLoader NuLoader(MCTargetFile);

    const unsigned int kNSel = InteractionTypes.size();
    Spectrum *spectra[kNVar][kNSel];

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        for(unsigned int jSel = 0; jSel < kNSel; ++jSel){
            spectra[iVar][jSel] = new Spectrum(SelectionPlots[iVar].label, 
                                               SelectionPlots[iVar].bins, 
                                               NuLoader, 
                                               SelectionPlots[iVar].var, 
                                               kCRTPMTNeutrino,
                                               kPreSelection && InteractionTypes[jSel].cut);          
        }
    }

    // Spectrum *sMCEventDump = new Spectrum("", Binning::Simple(3, 0, 3), NuLoader, kMCEventDump, kCRTPMTNeutrino); 

    NuLoader.Go();

    TFile FOut("CC1e0piSelection_Data.root", "recreate");

    TCanvas *c[kNVar];
    TLegend *l[kNVar];
    THStack *hs[kNVar];
    double dataPOT = dataSpectra[0]->POT(); ///< get POT from data, scale everything to that
    std::string title;

    for(unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        c[iVar] = new TCanvas(SelectionPlots[iVar].suffix.c_str(), SelectionPlots[iVar].suffix.c_str(), 500, 500);
        hs[iVar] = new THStack(SelectionPlots[iVar].suffix.c_str(), SelectionPlots[iVar].label.c_str());
        l[iVar] = new TLegend(0.65, 0.5, 0.85, 0.85, "NuMI CV");

        // set up data plot
        TH1* hData = dataSpectra[iVar]->ToTH1(dataPOT);
        float yMax = 0;
        for (int i = 1; i <= hData->GetNbinsX(); ++i) {
            double y = hData->GetBinContent(i);
            double err = hData->GetBinError(i);
            if (y + err > yMax) yMax = y + err;
        }

        // all slices with margins
        TH1* hAll = spectra[iVar][0]->ToTH1(dataPOT);
        
        // stack by interaction type
        for(unsigned int jSel = 1; jSel < kNSel; ++jSel) {
            TH1* h = spectra[iVar][jSel]->ToTH1(dataPOT);

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
                Form("Slices / %.1e POT", dataPOT);
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

        // plot data
        hData->SetMarkerStyle(20); 
        hData->SetLineWidth(1);
        hData->SetMarkerSize(0.75);
        hData->SetMarkerColor(kBlack);
        hData->Draw("EX0 SAME");
        l[iVar]->AddEntry(hData, "Data", "f");

        l[iVar]->SetTextSize(0.04);
        l[iVar]->Draw();
        c[iVar]->Write();

        c[iVar]->Update();
        c[iVar]->Modified();

        gStyle->SetLineScalePS(5);
        title = std::string("plots/") + SelectionPlots[iVar].suffix + std::string("_PreSelection_withData.pdf");
        c[iVar]->SaveAs(title.c_str());
    }

    FOut.Close();

    return;
}
