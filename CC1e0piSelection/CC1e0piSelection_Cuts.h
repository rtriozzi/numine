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

    // pre-selection
    const SpillCut kCRTPMTNeutrino([](const caf::SRSpillProxy* sr) {

        double minTime = 0., maxTime = 0.;

        for ( const auto& match: sr->crtpmt_matches ) {
            if (sr->hdr.ismc) { minTime = 0.0; maxTime = 9.8; }
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
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;

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
        if (largestShwIdx == -1) return false;
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[2].energy)) return false;

        return slc->reco.pfp[largestShwIdx].shw.plane[2].energy > 0.200;
    });

    const Cut kLargestRecoShower_dEdxCut([](const caf::SRSliceProxy* slc) { 
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return false;
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx)) return false;

        // return (slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx) > 0 && (slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx < 3.5);
        return slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx < 3.5;
    });

    const Cut kLargestRecoShower_OpenAngleCut([](const caf::SRSliceProxy* slc) { 
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return false;
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.open_angle)) return false;

        double openAngle = 180. * slc->reco.pfp[largestShwIdx].shw.open_angle / M_PI;

        return (openAngle < 10);
    });

    const Cut kLargestRecoShower_ConvGapCut([](const caf::SRSliceProxy* slc) { 
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return false;
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.conversion_gap)) return false;   

        return slc->reco.pfp[largestShwIdx].shw.conversion_gap < 5;
    });

    // Np0π selection
    const Cut kNSelectedProtons([](const caf::SRSliceProxy* slc) { 
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return false;

        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);

        return selectedProtonIdx.size() > 0;
    }); 

    const Cut kNoOtherParticle([](const caf::SRSliceProxy* slc) { 
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return false;

        int NOtherParticles = kNSelectedProtonsIdx_NOtherParticles(slc);
        if (NOtherParticles == -1) return false;

        return NOtherParticles == 0;
    }); 

    // automatic selection
    const Cut kAutomaticSelection = kNotClearCosmic && kVertexInFV && kTrigFlashMatch
                                    && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut 
                                    && kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut
                                    && kNSelectedProtons && kNoOtherParticle;

    // selections
    struct SelDef {
        std::string suffix = "";
        std::string label = "";
        Cut cut = kNoCut;
        int color = kBlack;
    };

    std::vector<SelDef> SelectionSteps = {
        {"presel", "Preselection",            kNotClearCosmic && kVertexInFV,     kRed+2},
        {"flash",  "Barycenter FM",           kNotClearCosmic && kVertexInFV && kTrigFlashMatch,     kRed-7},
        {"shower", "Electron ID",             kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut,   kOrange-3},
        {"showercuts1", "dE/dx",              kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut,   kGreen-2},
        {"showercuts2", "Angle and gap",      kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut && kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut,   kCyan-7},
        {"proton", "Proton ID",               kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut && kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut && kNSelectedProtons,   kBlue-4},
        {"nothingelse", "No other particle",  kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut && kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut && kNSelectedProtons && kNoOtherParticle,   kMagenta-3},
    };
}