#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <fstream>
#include <vector>
#include <math.h>
#include "TVector3.h"

#include "Vars.h"

namespace ana {

    // event level
    const Cut kNotClearCosmic([](const caf::SRSliceProxy* slc) { 
        return !(slc->is_clear_cosmic);
    });

    const Cut kVertexInFV([](const caf::SRSliceProxy* slc) { 
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;

        return kIsInFV(slc->vertex.x, slc->vertex.y, slc->vertex.z);
    });

    // michel cuts
    const Cut kMichel_EnergyCut([](const caf::SRSliceProxy* slc) { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return false;
        if (std::isnan(slc->reco.pfp[mhlIdx].shw.plane[2].energy)) return false;
        return (slc->reco.pfp[mhlIdx].shw.plane[2].energy >= 0 && 
                slc->reco.pfp[mhlIdx].shw.plane[2].energy <= 0.1);        
    });

    const Cut kMichel_MhlFracCut([](const caf::SRSliceProxy* slc) { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return false;
        if (std::isnan(slc->reco.pfp[mhlIdx].ngscore.mhl_frac)) return false;
        return (slc->reco.pfp[mhlIdx].ngscore.mhl_frac > 0.85);        
    });

    const Cut kMichel_StartInFiducial([](const caf::SRSliceProxy* slc) { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return false;
        if (std::isnan(slc->reco.pfp[mhlIdx].shw.start.x)) return false;

        // starts >10 cm away from boundaries
        return kIsInContainedInXVol(slc->reco.pfp[mhlIdx].shw.start.x, 
                                    slc->reco.pfp[mhlIdx].shw.start.y, 
                                    slc->reco.pfp[mhlIdx].shw.start.z,
                                    10.);        
    });  

    // mip cuts
    const Cut kMichel_CloseMip([](const caf::SRSliceProxy* slc) { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return false;
        const double minMipDist = kMichel_MinMipDistance(slc);
        if (minMipDist == -5) return false;

        return (minMipDist > 0) && (minMipDist < 10);        
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

    const Cut kMichelSelection = kMichel_EnergyCut && kMichel_StartInFiducial && kMichel_CloseMip && kMichel_MhlFracCut; 

    // selections
    struct SelDef {
        std::string suffix = "";
        std::string label = "";
        Cut cut = kNoCut;
        int color = kBlack;
    };
} 