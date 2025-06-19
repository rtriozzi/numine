#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <fstream>
#include <vector>
#include <math.h>

#include "TVector3.h"

#include "CC1e0piSelection_Vars.h"
#include "CC1e0piSelection_Cuts.h"

namespace ana {

    // signal definitions at the MCNeutrino level
    bool kIsTrueCC1e0pi(const caf::Proxy<caf::SRTrueInteraction>& nu) {

        bool kTrueNueCC = (nu.iscc) && 
                          (abs(nu.pdg) == 12) &&
                          kIsInFV(nu.position.x, nu.position.y, nu.position.z);

        int nPrimElectron = 0, nVisProtons = 0, nVisOther = 0;
        int vCryo = nu.position.x < 0 ? 0 : 1;
        int bestPlaneIdx = 2; ///< for now, just stick to Collection for calorimetry

        for (int ip(0); ip < nu.nprim ; ++ip) {

            // int bestPlaneIdx = nu.prim[ip].plane[vCryo][1].nhit > nu.prim[ip].plane[vCryo][2].nhit ? 1 : 2;

            // electron
            if (abs(nu.prim[ip].pdg) == 11) {
                if (nu.prim[ip].plane[vCryo][bestPlaneIdx].visE > 0.2) {
                    ++nPrimElectron;
                }
            }

            // protons
            if (nu.prim[ip].pdg == 2212) {
                if ((nu.prim[ip].plane[vCryo][bestPlaneIdx].visE > VISIBILTY_THRESHOLD_P) &&
                    kIsInContained(nu.prim[ip].end.x, nu.prim[ip].end.y, nu.prim[ip].end.z)) {
                        ++nVisProtons;
                    }
            }

            // pi0
            if (nu.prim[ip].pdg == 111) {
                ++nVisOther;
            }

            // neutrons (any number)
            if (nu.prim[ip].pdg == 2112) {
                continue;
            }

            // other particles (e.g., pions)
            if ((nu.prim[ip].pdg != 2212) && 
                (abs(nu.prim[ip].pdg) != 11) & 
                (nu.prim[ip].pdg != 2112)  && 
                (nu.prim[ip].pdg != 111)) {
                if (nu.prim[ip].plane[vCryo][bestPlaneIdx].visE >= VISIBILTY_THRESHOLD_PI) {
                    ++nVisOther;
                }
            }

        }

        return kTrueNueCC && (nPrimElectron == 1) && (nVisProtons > 0) && (nVisOther == 0);
    }

    // signal definitions at the slice level
    const Cut kTrueCC1e0pi([](const caf::SRSliceProxy* slc) { 

        if (std::isnan(slc->truth.position.x) || std::isnan(slc->truth.position.y) || std::isnan(slc->truth.position.z)) return false;

        bool kTrueNueCC = (slc->truth.index >= 0) && 
                          (slc->truth.iscc) && 
                          (abs(slc->truth.pdg) == 12) &&
                          kIsInFV(slc->truth.position.x, slc->truth.position.y, slc->truth.position.z);

        int nPrimElectron = 0, nVisProtons = 0, nVisOther = 0;
        int vCryo = slc->truth.position.x < 0 ? 0 : 1;
        int bestPlaneIdx = 2; ///< for now, just stick to Collection for calorimetry

        for (int ip(0); ip < slc->truth.nprim ; ++ip) {

            // int bestPlaneIdx = slc->truth.prim[ip].plane[vCryo][1].nhit > slc->truth.prim[ip].plane[vCryo][2].nhit ? 1 : 2;

            // electron
            if (abs(slc->truth.prim[ip].pdg) == 11) {
                if (slc->truth.prim[ip].plane[vCryo][bestPlaneIdx].visE > 0.200) {
                    ++nPrimElectron;
                }
            }

            // proton
            if (slc->truth.prim[ip].pdg == 2212) {
                if ((slc->truth.prim[ip].plane[vCryo][bestPlaneIdx].visE > VISIBILTY_THRESHOLD_P) &&
                    kIsInContained(slc->truth.prim[ip].end.x, slc->truth.prim[ip].end.y, slc->truth.prim[ip].end.z)) {
                       ++nVisProtons; 
                    }
            }

            // pi0
            if (slc->truth.prim[ip].pdg == 111) {
                ++nVisOther;
            }
            
            // neutrons (any number)
            if (slc->truth.prim[ip].pdg == 2112) {
                continue;
            }

            // other particles (e.g., pions)
            if ((slc->truth.prim[ip].pdg != 2212) && 
                (abs(slc->truth.prim[ip].pdg) != 11) && 
                (slc->truth.prim[ip].pdg != 2112) && 
                (slc->truth.prim[ip].pdg != 111)) {

                if (slc->truth.prim[ip].plane[vCryo][bestPlaneIdx].visE >= VISIBILTY_THRESHOLD_PI)
                    ++nVisOther;
            }

        }

        return kTrueNueCC && (nPrimElectron == 1) && (nVisProtons > 0) && (nVisOther == 0);
    });

    // const Cut kIsTrueCC1eXContained([](const caf::SRSliceProxy* slc) { 

    //     bool kTrueNueCC = (slc->truth.index >= 0) && 
    //                       (slc->truth.iscc) && 
    //                       (abs(slc->truth.pdg) == 12) &&
    //                       kIsInFV(slc->truth.position.x, slc->truth.position.y, slc->truth.position.z); 
    // });

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

    const Cut kIsNuMu([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.index >= 0) && (abs(slc->truth.pdg) == 14);
    });

    const Cut kIsCC([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.index >= 0) && (slc->truth.iscc);
    });

    const Cut kIsNC([](const caf::SRSliceProxy* slc) { 
        return (slc->truth.index >= 0) && (slc->truth.isnc);
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

    // true vertex was in FV?
    const Cut kTrueVertexInFV([](const caf::SRSliceProxy* slc) { 
        if (std::isnan(slc->truth.position.x) || std::isnan(slc->truth.position.y) || std::isnan(slc->truth.position.z)) return false;

        return (slc->truth.index >= 0) &&
               kIsInFV(slc->truth.position.x, slc->truth.position.y, slc->truth.position.z);
    });

    // for checking pile-up
    const SpillCut kNoPileUp([](const caf::SRSpillProxy* sr) {
        return (sr->mc.nnu == 1);
    });

    // selections
    std::vector<SelDef> InteractionTypes = {
        {"selected", "",                    kNoCut,  kBlack},
        {"signal", "CC1e^{#pm}0#pi",        kTrueCC1e0pi,     kRed-7},
        {"othernuecc", "#nu_{e}CC",         !kTrueCC1e0pi && kIsNue && kIsCC && kTrueVertexInFV,   kOrange-3},
        {"nuenc", "#nu_{e}NC",              !kTrueCC1e0pi && kIsNue && kIsNC && kTrueVertexInFV,   kGreen-2},
        {"numucc", "#nu_{#mu}CC#pi^{0}",    !kTrueCC1e0pi && kIsNuMu && kIsCC && kIsTherePi0 && kTrueVertexInFV,   kMagenta-10},
        {"numucc", "#nu_{#mu}CC",           !kTrueCC1e0pi && kIsNuMu && kIsCC && !kIsTherePi0 && kTrueVertexInFV,   kMagenta-3},
        {"numunc", "#nu_{#mu}NC",           !kTrueCC1e0pi && kIsNuMu && kIsNC && kTrueVertexInFV,   kPink+1},
        {"oofvnu", "OOFV",                  !kTrueCC1e0pi && kIsNuOOFV,   kCyan-9},
        {"ootcosmic", "Cosmic",             !kTrueCC1e0pi && kIsCosmic,   kAzure-3}
    };
}