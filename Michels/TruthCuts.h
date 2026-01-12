#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <fstream>
#include <vector>
#include <math.h>

#include "TVector3.h"

#include "Vars.h"
#include "Cuts.h"

namespace ana {

    // truth-level classification of slice
    const Cut kIsCosmic([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.index < 0);
    });

    const Cut kIsNu([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.index >= 0);
    });

    const Cut kIsNue([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.index >= 0) && (abs(slc->truth.pdg) == 12);
    });

    const Cut kIsNotNue([](const caf::SRSliceProxy* slc) { 
        return !(abs(slc->truth.pdg) == 12);
    });

    const Cut kIsNuMu([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.index >= 0) && (abs(slc->truth.pdg) == 14);
    });

    const Cut kIsCC([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.index >= 0) && (slc->truth.iscc);
    });

    const Cut kIsNC([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.index >= 0) && (slc->truth.isnc);
    });


    // truth-level classification of michels
    const Cut kIsMichel_E([](const caf::SRSliceProxy* slc) { 
        const int mhlPDG = kRecoMichelPdg(slc);
        if(mhlPDG == -5) return false;
        return abs(mhlPDG) == 11;
    });    

    const Cut kIsMichel_EPlus([](const caf::SRSliceProxy* slc) { 
        const int mhlPDG = kRecoMichelPdg(slc);
        if(mhlPDG == -5) return false;
        return mhlPDG == -11;
    });    

    const Cut kIsMichel_EMinus([](const caf::SRSliceProxy* slc) { 
        const int mhlPDG = kRecoMichelPdg(slc);
        if(mhlPDG == -5) return false;
        return mhlPDG == 11;
    });    

    const Cut kIsMichel_Ph([](const caf::SRSliceProxy* slc) { 
        const int mhlPDG = kRecoMichelPdg(slc);
        if(mhlPDG == -5) return false;
        return abs(mhlPDG) == 22;
    });   

    const Cut kIsMichel_P([](const caf::SRSliceProxy* slc) { 
        const int mhlPDG = kRecoMichelPdg(slc);
        if(mhlPDG == -5) return false;
        return abs(mhlPDG) == 2212;
    });   

    // selections
    std::vector<SelDef> ParticleTypes = {
        {"selected", "",        kNoCut, kBlack},
        {"electron", "e^{+}",   kIsMichel_EPlus, kMagenta-10},
        {"electron", "e^{-}",   kIsMichel_EMinus, kMagenta-3},
        {"photon", "#gamma",    kIsMichel_Ph, kGreen-2},
        {"other", "Other",      !kIsMichel_E && !kIsMichel_Ph, kPink+1},        
    };        
}