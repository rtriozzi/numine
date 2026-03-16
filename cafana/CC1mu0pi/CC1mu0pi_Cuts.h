#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <fstream>
#include <vector>
#include <math.h>
#include "TVector3.h"

#include "CC1mu0pi_Vars.h"

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

    const Cut kFlashMatch([](const caf::SRSliceProxy* slc) { 
        return (slc->barycenterFM.deltaZ >= 0 && 
                slc->barycenterFM.deltaZ <= 100 &&
                slc->barycenterFM.flashTime > -1 &&
                slc->barycenterFM.flashTime < 11);        
    });

    // muon identification
    const Cut kMuon_LengthCut([](const caf::SRSliceProxy* slc) { 
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx < 0) return false;
        if (std::isnan(slc->reco.pfp[muonIdx].trk.len)) return false;

        return slc->reco.pfp[muonIdx].trk.len >= 50;
    });

    const Cut kMuon_Containment([](const caf::SRSliceProxy* slc) { 
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx < 0) return false;
        if (std::isnan(slc->reco.pfp[muonIdx].trk.end.x) || std::isnan(slc->reco.pfp[muonIdx].trk.end.y) || std::isnan(slc->reco.pfp[muonIdx].trk.end.z)) 
            return false;

        return kIsInContained(slc->reco.pfp[muonIdx].trk.end.x, slc->reco.pfp[muonIdx].trk.end.y, slc->reco.pfp[muonIdx].trk.end.z);
    });

    // Np0π selection
    const Cut kNSelectedProtons([](const caf::SRSliceProxy* slc) { 
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx < 0) return false;

        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);

        return selectedProtonIdx.size() > 0;
    }); 

    const Cut kNoOtherParticle([](const caf::SRSliceProxy* slc) { 
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx <= 0) return false;

        int NOtherParticles = kNSelectedProtonsIdx_NOtherParticles(slc);
        if (NOtherParticles < 0) return false;

        return NOtherParticles == 0;
    }); 

    // automatic selection
    const Cut kAutomaticNuMuSelection = kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kMuon_LengthCut && kMuon_Containment &&
                                    kNSelectedProtons && kNoOtherParticle;

    const Cut kAutomaticNuMuSelection_NoTrigger = kNotClearCosmic && kVertexInFV && kFlashMatch && kMuon_LengthCut && kMuon_Containment &&
                                              kNSelectedProtons && kNoOtherParticle;

    const Cut kPreNuMuSelection = kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kMuon_LengthCut && kMuon_Containment;

    const Cut kPreNuMuSelection_NoTrigger = kNotClearCosmic && kVertexInFV && kFlashMatch && kMuon_LengthCut && kMuon_Containment;

    // selections
    struct SelDef {
        std::string suffix = "";
        std::string label = "";
        Cut cut = kNoCut;
        int color = kBlack;
    };

    std::vector<SelDef> SelectionSteps = {
        {"presel", "Presel.",               kNotClearCosmic && kVertexInFV,     kBlack},
        {"flash",  "FM",                    kNotClearCosmic && kVertexInFV && kTrigFlashMatch,     kRed+2},
        {"mulen", "Muon ID",                kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kMuon_LengthCut,     kRed-7},
        {"mucont", "Cont'd",                kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kMuon_LengthCut && kMuon_Containment,     kOrange-3},
        {"proton", "Proton ID",             kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kMuon_LengthCut && kMuon_Containment && kNSelectedProtons,   kGreen-2},
        {"nothingelse", "0#pi",             kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kMuon_LengthCut && kMuon_Containment && kNSelectedProtons && kNoOtherParticle,   kGreen-7},
    };

    std::vector<SelDef> SelectionSteps_NoTrigger = {
        {"presel", "Presel.",               kNotClearCosmic && kVertexInFV,     kBlack},
        {"flash",  "FM",                    kNotClearCosmic && kVertexInFV && kFlashMatch,     kRed+2},
        {"mulen", "Muon ID",                kNotClearCosmic && kVertexInFV && kFlashMatch && kMuon_LengthCut,     kRed-7},
        {"mucont", "Cont'd",                kNotClearCosmic && kVertexInFV && kFlashMatch && kMuon_LengthCut && kMuon_Containment,     kOrange-3},
        {"proton", "Proton ID",             kNotClearCosmic && kVertexInFV && kFlashMatch && kMuon_LengthCut && kMuon_Containment && kNSelectedProtons,   kGreen-2},
        {"nothingelse", "0#pi",             kNotClearCosmic && kVertexInFV && kFlashMatch && kMuon_LengthCut && kMuon_Containment && kNSelectedProtons && kNoOtherParticle,   kGreen-7},
    };
} 