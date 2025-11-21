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

    const Var kBarycenterFM_DeltaZ([](const caf::SRSliceProxy *slc) -> double {
        return slc->barycenterFM.deltaZ;
    });

    const Var kBarycenterFM_DeltaZ_Trigger([](const caf::SRSliceProxy *slc) -> double {
        return slc->barycenterFM.deltaZ_Trigger;
    });

    const Var kBarycenterFM_FlashTime([](const caf::SRSliceProxy *slc) -> double {
        return slc->barycenterFM.flashTime;
    });

    const Var kNuGraph_FilterFraction([](const caf::SRSliceProxy *slc) -> double {
        if (std::isnan(slc->ng_filt_pass_frac)) return -5.;
        return slc->ng_filt_pass_frac;
    });

    // michel identification
    const Var kMichelIdx([](const caf::SRSliceProxy* slc) -> int { 
        int mhlIdx(-1);
        double maxNHits(-1);

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if ((slc->reco.pfp[i].shw.plane[2].nHits > maxNHits) && (slc->reco.pfp[i].ngscore.sem_cat == 3)) {
                mhlIdx = i;
                maxNHits = slc->reco.pfp[i].shw.plane[2].nHits;
            }
        }

        return mhlIdx;
    });

    const Var kRecoMichelPdg([](const caf::SRSliceProxy* slc) -> int {
        const int mhlIdx = kMichelIdx(slc);
        if(mhlIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[mhlIdx].shw.truth.p.pdg)) 
            return -5;

        return slc->reco.pfp[mhlIdx].shw.truth.p.pdg;
    });

    // michel properties
    const Var kMichel_CollEnergy([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[mhlIdx].shw.plane[2].energy)) 
            return -5.;
            
        return slc->reco.pfp[mhlIdx].shw.plane[2].energy * 1.e3;
    });

    const Var kMichel_MhlFrac([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[mhlIdx].ngscore.mhl_frac)) 
            return -5.;
            
        return slc->reco.pfp[mhlIdx].ngscore.mhl_frac;
    });

    const Var kMichel_MipFrac([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[mhlIdx].ngscore.mhl_frac)) 
            return -5.;
            
        return slc->reco.pfp[mhlIdx].ngscore.mip_frac;
    });

    const Var kMichel_ShrFrac([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[mhlIdx].ngscore.shr_frac)) 
            return -5.;
            
        return slc->reco.pfp[mhlIdx].ngscore.shr_frac;
    });

    const Var kMichel_DifFrac([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[mhlIdx].ngscore.dif_frac)) 
            return -5.;
            
        return slc->reco.pfp[mhlIdx].ngscore.dif_frac;
    });

    const Var kMichel_HipFrac([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[mhlIdx].ngscore.hip_frac)) 
            return -5.;
            
        return slc->reco.pfp[mhlIdx].ngscore.hip_frac;
    });

    const Var kMichel_NCollHits([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[mhlIdx].shw.plane[2].nHits)) 
            return -5.;
            
        return slc->reco.pfp[mhlIdx].shw.plane[2].nHits;
    });

    const SpillMultiVar kMichel_TrueEnergy([](const caf::SRSpillProxy* sr) -> std::vector<double> { 
        for (const auto& p: sr->true_particles) {
            
        }
        for (int ip(0); ip < slc->truth.nprim ; ++ip) {
            if ((abs(slc->truth.prim[ip].pdg) == 11) && (slc->truth.prim[ip].start_process != 0)) {  
                return 1.e3*(slc->truth.prim[ip].startE - slc->truth.prim[ip].endE);
            }
        }
        return -5.;
    });

    // plotting
    struct PlotDef {
        std::string suffix = "";
        std::string label = "";
        Binning bins = Binning::Simple(3, 0, 3);
        Var var = kCounting;
    };

    std::vector<PlotDef> SelectionPlots = {   
        {"count", "Counts [#]",                       Binning::Simple(3, 0, 3), kCounting},
        {"collenergy", "E_{Coll} [MeV]",              Binning::Simple(40, 0, 100), kMichel_CollEnergy}, 
        {"mhlfrac", "NG2 mhl_frac",                   Binning::Simple(25, 0, 1), kMichel_MhlFrac}, 
        {"mipfrac", "NG2 mip_frac",                   Binning::Simple(25, 0, 1), kMichel_MipFrac}, 
        {"shrfrac", "NG2 shr_frac",                   Binning::Simple(25, 0, 1), kMichel_ShrFrac}, 
        {"hipfrac", "NG2 hip_frac",                   Binning::Simple(25, 0, 1), kMichel_HipFrac}, 
        {"diffrac", "NG2 dif_frac",                   Binning::Simple(25, 0, 1), kMichel_DifFrac}, 
        {"collnhits", "Coll. hits [#]",               Binning::Simple(25, 0, 100), kMichel_NCollHits}, 
        {"truemichenergy", "E^{true} [MeV]",          Binning::Simple(40, 0, 100), kMichel_TrueEnergy}
    };

}
