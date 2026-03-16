#pragma once

// sbnana stuff
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

// c++ stuff
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>

// root stuff
#include "TVector3.h"

namespace ana {

    const double VISIBILTY_THRESHOLD_P = 0.05;
    const double VISIBILTY_THRESHOLD_PI = 0.025;

    // general helper functions
    bool kIsInFV(double x, double y, double z) {  
        if (std::isnan(x) || std::isnan(y) || std::isnan(z)) return false;

        return (( ( x < -61.94 - 25 && x > -358.49 + 25 ) ||
            ( x >  61.94 + 25 && x <  358.49 - 25 )) &&
            ( ( y > -181.86 + 25 && y < 134.96 - 25 ) &&
            ( z > -894.95 + 30 && z < 894.95 - 50 ) ));
    }

    bool kIsInContained(double ex, double ey, double ez) { 
        if (std::isnan(ex) || std::isnan(ey) || std::isnan(ez)) return false;

        return (( ( ex < -61.94 - 5 && ex > -358.49 + 5 ) ||
            ( ex >  61.94 + 5 && ex <  358.49 - 5 )) &&
            ( ( ey > -181.86 + 5 && ey < 134.96 - 5 ) &&
            ( ez > -894.95 + 5 && ez < 894.95 - 5 ) ));
    }

    // general event variables
    const Var kCounting([](const caf::SRSliceProxy *slc) -> int {
        return 1;
    });

    const Var kVertex_vsTruth([](const caf::SRSliceProxy *slc) -> double {
            TVector3 VertexReco(slc->vertex.x, slc->vertex.y, slc->vertex.z);
            TVector3 VertexTrue(slc->truth.position.x, slc->truth.position.y, slc->truth.position.z);
            return (VertexReco-VertexTrue).Mag();
    });

    const Var kTrue_NVisProtons([](const caf::SRSliceProxy* slc) -> int { 
        int nVisProtons = 0;
        int vCryo = slc->truth.position.x < 0 ? 0 : 1;
        for (int ip(0); ip < slc->truth.nprim ; ++ip) {
            if ((slc->truth.prim[ip].pdg == 2212) &&
                ((slc->truth.prim[ip].startE - slc->truth.prim[ip].endE) > VISIBILTY_THRESHOLD_P)) 
                nVisProtons += 1;
        }
        return nVisProtons;
    });

    const Var kBarycenterFM_DeltaZ([](const caf::SRSliceProxy *slc) -> double {
        return slc->barycenterFM.deltaZ;
    });

    const Var kBarycenterFM_DeltaZ_Trigger([](const caf::SRSliceProxy *slc) -> double {
        return slc->barycenterFM.deltaZ_Trigger;
    });

    const Var kBarycenterFM_FlashTime([](const caf::SRSliceProxy *slc) -> double {
        return slc->barycenterFM.flashTime;
    });

    const Var kNuGraph_Ind1ShowerHits([](const caf::SRSliceProxy *slc) -> int {
        if (std::isnan(slc->ng_plane[0].shr_hits) || (slc->ng_plane[0].shr_hits < 0)) return -5;
        return slc->ng_plane[0].shr_hits;
    });

    const Var kNuGraph_Ind2ShowerHits([](const caf::SRSliceProxy *slc) -> int {
        if (std::isnan(slc->ng_plane[1].shr_hits) || (slc->ng_plane[1].shr_hits < 0)) return -5;
        return slc->ng_plane[1].shr_hits;
    });

    const Var kNuGraph_CollShowerHits([](const caf::SRSliceProxy *slc) -> int {
        if (std::isnan(slc->ng_plane[2].shr_hits) || (slc->ng_plane[2].shr_hits < 0)) return -5;
        return slc->ng_plane[2].shr_hits;
    });

    const Var kNuGraph_Ind1ShowerHits_Unclustered([](const caf::SRSliceProxy *slc) -> int {
        if (std::isnan(slc->ng_plane[0].unclustered_shr_hits) || (slc->ng_plane[0].unclustered_shr_hits < 0)) return -5;
        return slc->ng_plane[0].unclustered_shr_hits;
    });

    const Var kNuGraph_Ind2ShowerHits_Unclustered([](const caf::SRSliceProxy *slc) -> int {
        if (std::isnan(slc->ng_plane[1].unclustered_shr_hits) || (slc->ng_plane[1].unclustered_shr_hits < 0)) return -5;
        return slc->ng_plane[1].unclustered_shr_hits;
    });

    const Var kNuGraph_CollShowerHits_Unclustered([](const caf::SRSliceProxy *slc) -> int {
        if (std::isnan(slc->ng_plane[2].unclustered_shr_hits) || (slc->ng_plane[2].unclustered_shr_hits < 0)) return -5;
        return slc->ng_plane[2].unclustered_shr_hits;
    });

    // pion identification
    bool kIsPFPPionLike(const caf::SRSliceProxy* slc, unsigned int iPFP) {
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].trk.start.x) || std::isnan(slc->reco.pfp[iPFP].trk.start.y) || std::isnan(slc->reco.pfp[iPFP].trk.start.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].trk.end.x) || std::isnan(slc->reco.pfp[iPFP].trk.end.y) || std::isnan(slc->reco.pfp[iPFP].trk.end.z)) return false;

        TVector3 recoVertex(slc->vertex.x, slc->vertex.y, slc->vertex.z); 
        TVector3 recoStart(slc->reco.pfp[iPFP].trk.start.x, slc->reco.pfp[iPFP].trk.start.y, slc->reco.pfp[iPFP].trk.start.z);
        TVector3 startMomentum(slc->reco.pfp[iPFP].trk.dir.x * slc->reco.pfp[iPFP].trk.rangeP.p_pion,
                               slc->reco.pfp[iPFP].trk.dir.y * slc->reco.pfp[iPFP].trk.rangeP.p_pion, 
                               slc->reco.pfp[iPFP].trk.dir.z * slc->reco.pfp[iPFP].trk.rangeP.p_pion); 
        double K = sqrt(pow(0.139570, 2) + pow(startMomentum.Mag(), 2)); ///< GeV

        return ((recoStart - recoVertex).Mag() < 10) &&
               (K >= VISIBILTY_THRESHOLD_PI);
    }

    // generic shower identification
    bool kIsPFPShowerLike(const caf::SRSliceProxy* slc, unsigned int iPFP) {
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].shw.start.x) || std::isnan(slc->reco.pfp[iPFP].shw.start.y) || std::isnan(slc->reco.pfp[iPFP].shw.start.z)) return false;

        TVector3 recoVertex(slc->vertex.x, slc->vertex.y, slc->vertex.z); 
        TVector3 recoStart(slc->reco.pfp[iPFP].shw.start.x, slc->reco.pfp[iPFP].shw.start.y, slc->reco.pfp[iPFP].shw.start.z);

        return ((recoStart - recoVertex).Mag() < 50) &&
               (slc->reco.pfp[iPFP].shw.plane[2].energy >= VISIBILTY_THRESHOLD_PI);
    }

    // proton identification
    bool kIsPFPProtonLike(const caf::SRSliceProxy* slc, unsigned int iPFP) {
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].trk.start.x) || std::isnan(slc->reco.pfp[iPFP].trk.start.y) || std::isnan(slc->reco.pfp[iPFP].trk.start.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].trk.end.x) || std::isnan(slc->reco.pfp[iPFP].trk.end.y) || std::isnan(slc->reco.pfp[iPFP].trk.end.z)) return false;

        TVector3 recoVertex(slc->vertex.x, slc->vertex.y, slc->vertex.z); 
        TVector3 recoStart(slc->reco.pfp[iPFP].trk.start.x, slc->reco.pfp[iPFP].trk.start.y, slc->reco.pfp[iPFP].trk.start.z);
        TVector3 startMomentum(slc->reco.pfp[iPFP].trk.dir.x * slc->reco.pfp[iPFP].trk.rangeP.p_proton,
                               slc->reco.pfp[iPFP].trk.dir.y * slc->reco.pfp[iPFP].trk.rangeP.p_proton, 
                               slc->reco.pfp[iPFP].trk.dir.z * slc->reco.pfp[iPFP].trk.rangeP.p_proton); 
        double K = sqrt(pow(0.9383, 2) + pow(startMomentum.Mag(), 2)); ///< GeV

        return kIsInContained(slc->reco.pfp[iPFP].trk.end.x, slc->reco.pfp[iPFP].trk.end.y, slc->reco.pfp[iPFP].trk.end.z) &&
               ((recoStart - recoVertex).Mag() < 10) &&
               (K >= VISIBILTY_THRESHOLD_P);
    }

    // plotting
    struct PlotDef {
        std::string suffix = "";
        std::string label = "";
        Binning bins = Binning::Simple(3, 0, 3);
        Var var = kCounting;
    };
}
