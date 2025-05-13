#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <fstream>
#include <vector>
#include <math.h>

#include "TVector3.h"

#include "CC1e0piSelection_Vars.h"

namespace ana {

    // general helper functions
    bool kIsInFV(double x, double y, double z) {  
        if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;

        return (( ( x < -61.94 - 25 && x > -358.49 + 25 ) ||
            ( x >  61.94 + 25 && x <  358.49 - 25 )) &&
            ( ( y > -181.86 + 25 && y < 134.96 - 25 ) &&
            ( z > -894.95 + 30 && z < 894.95 - 50 ) ));
    }

    // pre-selection
    const SpillCut kCRTPMTNeutrino([](const caf::SRSpillProxy* sr) {

        double minTime = 0., maxTime = 0.;

        for ( const auto& match: sr->crtpmt_matches ) {
            if (sr->hdr.ismc) { minTime = 0.0; maxTime = 9.7; }
            if (!sr->hdr.ismc) { minTime = -0.4; maxTime = 10.5; }

            // in-time flash and CRT-PMT match
            if(match.flashGateTime > minTime && 
               match.flashGateTime < maxTime && 
               match.flashClassification == 0)  
                return true;
        }

        return false;
    });

    const Cut kNotClearCosmic([](const caf::SRSliceProxy* slc) { 
        return !(slc->is_clear_cosmic);
    });

    const Cut kVertexInFV([](const caf::SRSliceProxy* slc) { 
        return kIsInFV(slc->vertex.x, slc->vertex.y, slc->vertex.z);
    });

    // light matching
    const Cut kTrigFlashMatch([](const caf::SRSliceProxy* slc) { 
        return (slc->barycenterFM.deltaZ_Trigger >= 0 && 
                slc->barycenterFM.deltaZ_Trigger <= 100);        
    });

    // electron identification
    const Cut kLargestRecoShower_EnergyCut([](const caf::SRSliceProxy* slc) { 
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return false;

        return slc->reco.pfp[largestShwIdx].shw.plane[2].energy > 0.200;
    });

    const Cut kLargestRecoShower_dEdxCut([](const caf::SRSliceProxy* slc) { 
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return false;

        return slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx < 3.5;
    });

    const Cut kLargestRecoShower_OpenAngleCut([](const caf::SRSliceProxy* slc) { 
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return false;

        double openAngle = 180. * slc->reco.pfp[largestShwIdx].shw.open_angle / M_PI;

        return (openAngle < 10);
    });

    const Cut kLargestRecoShower_ConvGapCut([](const caf::SRSliceProxy* slc) { 
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return false;

        return slc->reco.pfp[largestShwIdx].shw.conversion_gap < 5;
    });

    // automatic selection
    const Cut kAutomaticSelection = kNotClearCosmic && kVertexInFV && kTrigFlashMatch
                                    && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut 
                                    && kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut;

    // selections
    struct SelDef {
        std::string suffix = "";
        std::string label = "";
        Cut cut = kNoCut;
        int color = kBlack;
    };

    // std::vector<SelDef> InteractionTypes = {
    //     {"selected", "Selected 1eNp",                           k1eNpAutomaticSelection,  kGray+2},
    //     {"signal", "Signal: True 1eNp",                         k1eNpAutomaticSelection && kTrueSelection1eNpContained,     kRed-7},
    //     {"othernuecc", "(#nu_{e} + anti-#nu_{e})CC",            k1eNpAutomaticSelection && !kTrueSelection1eNpContained && kIsNue && kIsCC && kTrueVertexInFV,   kOrange-3},
    //     {"nuenc", "(#nu_{e} + anti-#nu_{e})NC",                 k1eNpAutomaticSelection && !kTrueSelection1eNpContained && kIsNue && kIsNC && kTrueVertexInFV,   kGreen-2},
    //     {"numucc", "(#nu_{#mu} + anti-#nu_{#mu})CC",            k1eNpAutomaticSelection && !kTrueSelection1eNpContained && kIsNuMu && kIsCC && kTrueVertexInFV,   kMagenta-3},
    //     {"numunc", "(#nu_{#mu} + anti-#nu_{#mu})NC",            k1eNpAutomaticSelection && !kTrueSelection1eNpContained && kIsNuMu && kIsNC && kTrueVertexInFV,   kMagenta-10},
    //     {"oofvnu", "Out-of-FV #nu",                             k1eNpAutomaticSelection && !kTrueSelection1eNpContained && kIsNuOFV,   kCyan-9},
    //     {"ootcosmic", "Out-of-time Cosmic",                     k1eNpAutomaticSelection && !kTrueSelection1eNpContained && kIsCosmic,   kAzure-3}
    // };

    std::vector<SelDef> SelectionSteps = {
        {"presel", "Pre-sel. (NCC, FV)",                 kNotClearCosmic && kVertexInFV,     kRed+2},
        {"flash", "Barycenter FM",                       kNotClearCosmic && kVertexInFV && kTrigFlashMatch,     kRed-7},
        {"shower", "Electron ID, E > 200 MeV, TS < 0.5", kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut,   kOrange-3},
        {"showercuts1", "dE/dx < 3.5 MeV/cm",            kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut,   kGreen-2},
        {"showercuts2", "Open. angle < 10 deg.",         kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut && kLargestRecoShower_OpenAngleCut,   kCyan-7},
        {"showercuts3", "Conv. gap < 5 cm",              kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut && kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut,   kBlue-4},
    };
}