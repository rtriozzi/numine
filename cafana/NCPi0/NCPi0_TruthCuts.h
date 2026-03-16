#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <fstream>
#include <vector>
#include <math.h>

#include "TVector3.h"

#include "NCPi0_Vars.h"
#include "NCPi0_Cuts.h"

namespace ana {

    // signal definitions at the MCNeutrino level
    bool kIsTrueNCPi0(const caf::Proxy<caf::SRTrueInteraction>& nu) {

        bool kTrueNC = (nu.isnc) && kIsInFV(nu.position.x, nu.position.y, nu.position.z);

        int nPi0 = 0;
        int vCryo = nu.position.x < 0 ? 0 : 1;
        int bestPlaneIdx = 2; ///< for now, just stick to Collection for calorimetry

        for (int ip(0); ip < nu.nprim ; ++ip) {
            // pi0
            if (nu.prim[ip].pdg == 111) {
                ++nPi0;
            }
        }

        return kTrueNC && (nPi0 == 1);
    }

    // signal definitions at the slice level
    const Cut kTrueNCPi0([](const caf::SRSliceProxy* slc) { 

        if (std::isnan(slc->truth.position.x) || std::isnan(slc->truth.position.y) || std::isnan(slc->truth.position.z)) return false;

        bool kTrueNC = (slc->truth.index >= 0) && 
                       (slc->truth.isnc) && 
                       kIsInFV(slc->truth.position.x, slc->truth.position.y, slc->truth.position.z);

        int nPi0 = 0;
        int vCryo = slc->truth.position.x < 0 ? 0 : 1;
        int bestPlaneIdx = 2; ///< for now, just stick to Collection for calorimetry

        for (int ip(0); ip < slc->truth.nprim ; ++ip) {
            // pi0
            if (slc->truth.prim[ip].pdg == 111) {
                ++nPi0;
            }
        }

        return kTrueNC && (nPi0 == 1);
    });

    // truth-level classification of slice
    const Cut kIsNuOOFV([](const caf::SRSliceProxy* slc) { 
        if (std::isnan(slc->truth.position.x) || std::isnan(slc->truth.position.y) || std::isnan(slc->truth.position.z)) return false;

        return (slc->truth.index >= 0) &&
            !kIsInFV(slc->truth.position.x, slc->truth.position.y, slc->truth.position.z);
    });

    const Cut kIsNuinFV([](const caf::SRSliceProxy* slc) { 
        if (std::isnan(slc->truth.position.x) || std::isnan(slc->truth.position.y) || std::isnan(slc->truth.position.z)) return false;

        return (slc->truth.index >= 0) &&
            kIsInFV(slc->truth.position.x, slc->truth.position.y, slc->truth.position.z);       
    });

    const Cut kIsCosmic([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.index < 0);
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

    const Cut kIsThereOnePi0([](const caf::SRSliceProxy* slc) { 
        int NPi0 = 0;
        for (int ip(0); ip < slc->truth.nprim ; ++ip) {
            if (slc->truth.prim[ip].pdg == 111) {
                ++NPi0;
            }
        }

        return (slc->truth.index >= 0) && (NPi0 == 1);
    });

    const Cut kIsTherePi0([](const caf::SRSliceProxy* slc) { 
        int NPi0 = 0;
        for (int ip(0); ip < slc->truth.nprim ; ++ip) {
            if (slc->truth.prim[ip].pdg == 111) {
                ++NPi0;
            }
        }

        return (slc->truth.index >= 0) && (NPi0 > 0);
    });

    const Cut kTrueQE([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.genie_mode == caf::kQE);
    });

    const Cut kTrueDIS([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.genie_mode == caf::kDIS);
    });

    const Cut kTrueRes([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.genie_mode == caf::kRes);
    });

    const Cut kTrueCoh([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.genie_mode == caf::kCoh);
    });

    const Cut kTrueMEC([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.genie_mode == caf::kMEC);
    });

   // truth-level classification of reconstructed leading shower
    const Cut kIsLargestShower_E([](const caf::SRSliceProxy* slc) { 
        const int largestShwPDG = kLargestRecoShower_TruePdg(slc);
        if(largestShwPDG == -5) return false;
        return abs(largestShwPDG) == 11;
    });    

    const Cut kIsLargestShower_Mu([](const caf::SRSliceProxy* slc) { 
        const int largestShwPDG = kLargestRecoShower_TruePdg(slc);
        if(largestShwPDG == -5) return false;
        return abs(largestShwPDG) == 13;
    });    

    const Cut kIsLargestShower_Pi([](const caf::SRSliceProxy* slc) { 
        const int largestShwPDG = kLargestRecoShower_TruePdg(slc);
        if(largestShwPDG == -5) return false;
        return abs(largestShwPDG) == 211;
    });   

    const Cut kIsLargestShower_Ph([](const caf::SRSliceProxy* slc) { 
        const int largestShwPDG = kLargestRecoShower_TruePdg(slc);
        if(largestShwPDG == -5) return false;
        return abs(largestShwPDG) == 22;
    });   

    const Cut kIsLargestShower_P([](const caf::SRSliceProxy* slc) { 
        const int largestShwPDG = kLargestRecoShower_TruePdg(slc);
        if(largestShwPDG == -5) return false;
        return abs(largestShwPDG) == 2212;
    });   

    // true vertex was in FV?
    const Cut kTrueVertexInFV([](const caf::SRSliceProxy* slc) { 
        if (std::isnan(slc->truth.position.x) || std::isnan(slc->truth.position.y) || std::isnan(slc->truth.position.z)) return false;

        return (slc->truth.index >= 0) &&
               kIsInFV(slc->truth.position.x, slc->truth.position.y, slc->truth.position.z);
    });

    // selections
    std::vector<SelDef> InteractionTypes = {
        {"selected", "",                    kNoCut,  kBlack},
        // NC
        {"signal", "NC1#pi^{0}",            kTrueNCPi0,     kRed-7},
        {"ncmorepi0", "NC>1#pi^{0}",        !kTrueNCPi0 && kIsNC && !kIsThereOnePi0 && kIsTherePi0 && kTrueVertexInFV,   kOrange-3},
        {"othernc", "Other NC",             !kTrueNCPi0 && kIsNC && !kIsTherePi0 && kTrueVertexInFV,   kGreen-2},
        // CC
        {"ccpi0", "CC#pi^{0}",              !kTrueNCPi0 && kIsCC && kIsCC && kIsTherePi0 && kTrueVertexInFV,   kMagenta-10},
        {"othercc", "Other CC",             !kTrueNCPi0 && kIsCC && !kIsTherePi0 && kTrueVertexInFV,   kMagenta-3},
        {"oofvnu", "OoFV",                  !kTrueNCPi0 && kIsNuOOFV,   kCyan-9},
        {"ootcosmic", "Cosmic",             !kTrueNCPi0 && kIsCosmic,   kAzure-3}
    };

   std::vector<SelDef> LeadingShowerParticleTypes = {
        {"selected", "",        kNoCut,  kBlack},
        {"electron", "e^{#pm}", kIsLargestShower_E,     kRed-7},
        {"photon", "#gamma",    kIsLargestShower_Ph,   kOrange-3},
        {"muon", "#mu^{#pm}",   kIsLargestShower_Mu,   kGreen-2},
        {"pion", "#pi^{#pm}",   kIsLargestShower_Pi,   kMagenta-10},
        {"proton", "p",         kIsLargestShower_P,   kMagenta-3},
        {"other", "Other",      !kIsLargestShower_E && !kIsLargestShower_Ph && !kIsLargestShower_Mu && !kIsLargestShower_Pi && !kIsLargestShower_P,   kPink+1},        
    };        
}