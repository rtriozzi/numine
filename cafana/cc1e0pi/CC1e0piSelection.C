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

void CC1e0piSelection() {

    // CNAF NuMI MC
    // const std::string TargetFile = "/storage/gpfs_data/icarus/local/users/cfarnese/NUMI/NUMI_MC/*.root"; ///< CV
    
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/user/cfarnese/RT_production_NuMIreference_20May25/mc*/caf_here/*.flat.caf.root"; ///< reference
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/user/cfarnese/RT_NUMI_nuonly_May18/run*/cafmakerjob_here_2d_updated/*.flat.caf.root"; ///< new BDT
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/user/cfarnese/RT_production_NuMIcheating_20May25/mc*/caf_here/*.flat.caf.root"; ///< vertex cheated
    // const std::string TargetFile = "/storage/gpfs_data/icarus/plain/user/cfarnese/RT_production_NuMIcheatingnew_21May25/mc*/caf_here/*.flat.caf.root"; ///< BDT vertex closest to truth
    // const std::string TargetFile = "/storage/gpfs_data/icarus/local/users/rtriozzi/concats/NuMI_CV_MopUp_NewBDT/*.root"; ///< Mop-up and then new BDT (trained without mop-up)

    // FNAL NuMI MC
    // const std::string TargetFile = "/exp/icarus/data/users/rtriozzi/mc/numi_FRFIX/concat_NuMI_MC_FRFIX_*.root";

    // FNAL NuGraph MC
    // const std::string TargetFile = "/pnfs/sbn/data/sbn_fd/poms_production/mc/2025A_ICARUS_NuGraph2/AddedNuGraph2_NuMI_sample_29May2025/v10_06_00_01p01/haddedFlatcaf/*.root";

    // FNAL development NuMI MC
    // const std::string TargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/standard/numinue.flat.caf.root"; ///< NuE
    // const std::string TargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/standard/numinom.flat.caf.root"; ///< nominal flux, mostly NuMu 
    // const std::string TargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco/numinue_NuGraphReco.flat.caf.root"; ///< NuE
    // const std::string TargetFile = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco/numinom_NuGraphReco.flat.caf.root"; ///< nominal flux, mostly NuMu 

    // CNAF nue-dis
    const std::string TargetFile = "/storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/concats/cv/cv_run*.flat.caf.root"; ///< nominal flux

    // detector variations...
    // const std::string TargetFile = "/storage/gpfs_data/icarus/local/users/rtriozzi/nuedis/concats/var1_hitcohnoise/var1_run*.flat.caf.root"; ///< high coherent noise

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
                                               kNoSpillCut, // kCRTPMTNeutrino, //kCRTPMTNeutrino && kNoPileUp,
                                               kAutomaticSelection && InteractionTypes[jSel].cut);  ///< change selection here if needed    
        }
    }

    Spectrum *sEventDump = new Spectrum("", Binning::Simple(3, 0, 3), NuLoader, kMCEventDump_DebugSelection, kNoSpillCut); 

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
