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

    bool kIsInContainedInXVol(double ex, double ey, double ez, double X) { 
        if (std::isnan(ex) || std::isnan(ey) || std::isnan(ez)) return false;

        return (( ( ex < -61.94 - X && ex > -358.49 + X ) ||
            ( ex >  61.94 + X && ex <  358.49 - X )) &&
            ( ( ey > -181.86 + X && ey < 134.96 - X ) &&
            ( ez > -894.95 + X && ez < 894.95 - X ) ));
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

    const Var kMichel_Ind2Energy([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[mhlIdx].shw.plane[1].energy)) 
            return -5.;
            
        return slc->reco.pfp[mhlIdx].shw.plane[1].energy * 1.e3;
    });

    const Var kMichel_Ind1Energy([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[mhlIdx].shw.plane[0].energy)) 
            return -5.;
            
        return slc->reco.pfp[mhlIdx].shw.plane[0].energy * 1.e3;
    });

    const Var kMichel_CollEnergy_VsTruth([](const caf::SRSliceProxy* slc) -> double {
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        double recoEnergy = kMichel_CollEnergy(slc);
        if (recoEnergy == -5) return -5.;

        if (std::isnan(slc->reco.pfp[mhlIdx].shw.truth.p.startE) || std::isnan(slc->reco.pfp[mhlIdx].shw.truth.p.endE)) return -5.;
        double trueEnergy = 1.e3 * (slc->reco.pfp[mhlIdx].shw.truth.p.startE - slc->reco.pfp[mhlIdx].shw.truth.p.endE);

        if (trueEnergy == 0) return -5;
        return (recoEnergy - trueEnergy) / trueEnergy;
    });

    const Var kLargestRecoShower_ColldEdx([](const caf::SRSliceProxy* slc) -> double {
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        if(std::isnan(slc->reco.pfp[mhlIdx].shw.plane[2].dEdx)) return -5;
        return slc->reco.pfp[mhlIdx].shw.plane[2].dEdx;
    });

    const Var kLargestRecoShower_Ind2dEdx([](const caf::SRSliceProxy* slc) -> double {
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        if(std::isnan(slc->reco.pfp[mhlIdx].shw.plane[1].dEdx)) return -5;
        return slc->reco.pfp[mhlIdx].shw.plane[1].dEdx;
    });

    const Var kLargestRecoShower_Ind1dEdx([](const caf::SRSliceProxy* slc) -> double {
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        if(std::isnan(slc->reco.pfp[mhlIdx].shw.plane[0].dEdx)) return -5;
        return slc->reco.pfp[mhlIdx].shw.plane[0].dEdx;
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

    const Var kMichel_NInd2Hits([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[mhlIdx].shw.plane[1].nHits)) 
            return -5.;
            
        return slc->reco.pfp[mhlIdx].shw.plane[1].nHits;
    });

    const Var kMichel_NInd1Hits([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[mhlIdx].shw.plane[0].nHits)) 
            return -5.;
            
        return slc->reco.pfp[mhlIdx].shw.plane[0].nHits;
    });

    // relationship with close particles
    const Var kMichel_MinMipDistance([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[mhlIdx].shw.start.x)) return -5.;
        TVector3 mhlStart(slc->reco.pfp[mhlIdx].shw.start.x, slc->reco.pfp[mhlIdx].shw.start.y, slc->reco.pfp[mhlIdx].shw.start.z);

        double minMipDist(9999.);
        TVector3 pStart, pEnd;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (i == (unsigned int) mhlIdx) continue;
            if (slc->reco.pfp[i].ngscore.sem_cat == 0) {
                if (std::isnan(slc->reco.pfp[i].trk.start.x) || std::isnan(slc->reco.pfp[i].trk.end.x)) return -5.;
                pStart.SetXYZ(slc->reco.pfp[i].trk.start.x, slc->reco.pfp[i].trk.start.y, slc->reco.pfp[i].trk.start.z);
                pEnd.SetXYZ(slc->reco.pfp[i].trk.end.x, slc->reco.pfp[i].trk.end.y, slc->reco.pfp[i].trk.end.z);
                if (std::min((pStart - mhlStart).Mag(), (pEnd - mhlStart).Mag()) < minMipDist) 
                    minMipDist = std::min((pStart - mhlStart).Mag(), (pEnd - mhlStart).Mag()); ///< account of swapped particles
            }
        }

        return minMipDist;
    });

    const Var kMichel_MinMipAngle([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[mhlIdx].shw.start.x) || std::isnan(slc->reco.pfp[mhlIdx].shw.dir.x)) return -5.;
        TVector3 mhlStart(slc->reco.pfp[mhlIdx].shw.start.x, slc->reco.pfp[mhlIdx].shw.start.y, slc->reco.pfp[mhlIdx].shw.start.z);
        TVector3 mhlStartDir(slc->reco.pfp[mhlIdx].shw.dir.x, slc->reco.pfp[mhlIdx].shw.dir.y, slc->reco.pfp[mhlIdx].shw.dir.z);

        double minMipDist(9999.);
        int minMipIdx(-1);
        TVector3 pStart, pEnd;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (i == (unsigned int) mhlIdx) continue;
            if (slc->reco.pfp[i].ngscore.sem_cat == 0) {
                if (std::isnan(slc->reco.pfp[i].trk.start.x) || std::isnan(slc->reco.pfp[i].trk.end.x)) return -5.;
                pStart.SetXYZ(slc->reco.pfp[i].trk.start.x, slc->reco.pfp[i].trk.start.y, slc->reco.pfp[i].trk.start.z);
                pEnd.SetXYZ(slc->reco.pfp[i].trk.end.x, slc->reco.pfp[i].trk.end.y, slc->reco.pfp[i].trk.end.z);
                if (std::min((pStart - mhlStart).Mag(), (pEnd - mhlStart).Mag()) < minMipDist) {
                    minMipDist = std::min((pStart - mhlStart).Mag(), (pEnd - mhlStart).Mag()); ///< account of swapped particles
                    minMipIdx = i;
                }
            }
        }
        
        if (minMipIdx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[mhlIdx].trk.dir_end.x)) return -5.;
        TVector3 pEndDir(slc->reco.pfp[minMipIdx].trk.dir_end.x, slc->reco.pfp[minMipIdx].trk.dir_end.y, slc->reco.pfp[minMipIdx].trk.dir_end.z);

        return mhlStartDir.Dot(pEndDir) / (mhlStartDir.Mag() * pEndDir.Mag());
    });

    const Var kMichel_CloseShoweryBlips([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[mhlIdx].shw.start.x)) return -5.;
        TVector3 mhlStart(slc->reco.pfp[mhlIdx].shw.start.x, slc->reco.pfp[mhlIdx].shw.start.y, slc->reco.pfp[mhlIdx].shw.start.z);

        int maxBlipDist(30); ///< cm
        int nCloseShoweryBlips(0);
        TVector3 pStart, pEnd;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (i == (unsigned int) mhlIdx) continue;
            if ((slc->reco.pfp[i].ngscore.sem_cat == 2) ||
                (slc->reco.pfp[i].ngscore.sem_cat == 3) ||
                (slc->reco.pfp[i].ngscore.sem_cat == 4)) {
                if (std::isnan(slc->reco.pfp[i].shw.start.x) || std::isnan(slc->reco.pfp[i].shw.end.x)) return -5.;
                pStart.SetXYZ(slc->reco.pfp[i].shw.start.x, slc->reco.pfp[i].shw.start.y, slc->reco.pfp[i].shw.start.z);
                pEnd.SetXYZ(slc->reco.pfp[i].shw.end.x, slc->reco.pfp[i].shw.end.y, slc->reco.pfp[i].shw.end.z);
                if (std::min((pStart - mhlStart).Mag(), (pEnd - mhlStart).Mag()) < maxBlipDist) 
                    nCloseShoweryBlips += 1;
            }
        }

        return nCloseShoweryBlips;
    });

    const Var kMichel_CloseMichelHits([](const caf::SRSliceProxy* slc) -> double { 
        const int mhlIdx = kMichelIdx(slc);
        if (mhlIdx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[mhlIdx].shw.start.x)) return -5.;
        TVector3 mhlStart(slc->reco.pfp[mhlIdx].shw.start.x, slc->reco.pfp[mhlIdx].shw.start.y, slc->reco.pfp[mhlIdx].shw.start.z);

        int maxDist(30); ///< cm
        int nCloseMichelHits(0);
        TVector3 pStart, pEnd;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (i == (unsigned int) mhlIdx) continue;
            if ((slc->reco.pfp[i].ngscore.sem_cat == 2) ||
                (slc->reco.pfp[i].ngscore.sem_cat == 3) ||
                (slc->reco.pfp[i].ngscore.sem_cat == 4)) {
                if (std::isnan(slc->reco.pfp[i].shw.start.x) || std::isnan(slc->reco.pfp[i].shw.end.x)) return -5.;
                pStart.SetXYZ(slc->reco.pfp[i].shw.start.x, slc->reco.pfp[i].shw.start.y, slc->reco.pfp[i].shw.start.z);
                pEnd.SetXYZ(slc->reco.pfp[i].shw.end.x, slc->reco.pfp[i].shw.end.y, slc->reco.pfp[i].shw.end.z);
                
                if (std::min((pStart - mhlStart).Mag(), (pEnd - mhlStart).Mag()) < maxDist) 
                    nCloseMichelHits += slc->reco.pfp[i].ngscore.mhl_frac *
                                        (slc->reco.pfp[i].shw.plane[0].nHits + slc->reco.pfp[i].shw.plane[1].nHits + slc->reco.pfp[i].shw.plane[2].nHits);
            }
            else if ((slc->reco.pfp[i].ngscore.sem_cat == 0) ||
                     (slc->reco.pfp[i].ngscore.sem_cat == 1)) {
                if (std::isnan(slc->reco.pfp[i].trk.start.x) || std::isnan(slc->reco.pfp[i].trk.end.x)) return -5.;
                pStart.SetXYZ(slc->reco.pfp[i].trk.start.x, slc->reco.pfp[i].trk.start.y, slc->reco.pfp[i].trk.start.z);
                pEnd.SetXYZ(slc->reco.pfp[i].trk.end.x, slc->reco.pfp[i].trk.end.y, slc->reco.pfp[i].trk.end.z);

                if (std::min((pStart - mhlStart).Mag(), (pEnd - mhlStart).Mag()) < maxDist) 
                    nCloseMichelHits += slc->reco.pfp[i].ngscore.mhl_frac *
                                        (slc->reco.pfp[i].trk.calo[0].nhit + slc->reco.pfp[i].trk.calo[1].nhit + slc->reco.pfp[i].trk.calo[2].nhit);

            }
        }

        return nCloseMichelHits;
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
        {"collenergy", "E_{Coll} [MeV]",              Binning::Simple(40, 0, 80), kMichel_CollEnergy}, 
        {"ind2energy", "E_{I2} [MeV]",                Binning::Simple(40, 0, 80), kMichel_Ind2Energy}, 
        {"ind1energy", "E_{I1} [MeV]",                Binning::Simple(40, 0, 80), kMichel_Ind1Energy}, 
        {"collenergyvstruth", "#DeltaE_{Coll} / E",   Binning::Simple(40, -1, 1), kMichel_CollEnergy_VsTruth},
        {"colldedx", "dE/dx_{Coll} [MeV/cm]",         Binning::Simple(40, 0, 8), kLargestRecoShower_ColldEdx}, 
        {"ind2dedx", "dE/dx_{I2} [MeV/cm]",           Binning::Simple(40, 0, 8), kLargestRecoShower_Ind2dEdx}, 
        {"ind1dedx", "dE/dx_{I1} [MeV/cm]",           Binning::Simple(40, 0, 8), kLargestRecoShower_Ind1dEdx}, 
        {"mhlfrac", "NG2 mhl_frac",                   Binning::Simple(40, 0, 1), kMichel_MhlFrac}, 
        {"mipfrac", "NG2 mip_frac",                   Binning::Simple(40, 0, 1), kMichel_MipFrac}, 
        {"shrfrac", "NG2 shr_frac",                   Binning::Simple(40, 0, 1), kMichel_ShrFrac}, 
        {"hipfrac", "NG2 hip_frac",                   Binning::Simple(40, 0, 1), kMichel_HipFrac}, 
        {"diffrac", "NG2 dif_frac",                   Binning::Simple(40, 0, 1), kMichel_DifFrac}, 
        {"collnhits", "Coll. hits [#]",               Binning::Simple(40, 0, 80), kMichel_NCollHits}, 
        {"ind2nhits", "I2 hits [#]",                  Binning::Simple(40, 0, 80), kMichel_NInd2Hits}, 
        {"ind1nhits", "I1 hits [#]",                  Binning::Simple(40, 0, 80), kMichel_NInd1Hits}, 
        {"minmipdist", "Dist. to closest MIP [cm]",   Binning::Simple(40, 0, 10), kMichel_MinMipDistance},
        {"minmipangle", "Angle to closest MIP",       Binning::Simple(40, -1, 1), kMichel_MinMipAngle},
        {"closeblips", "N. close showery blips [#]",  Binning::Simple(10, 0, 10), kMichel_CloseShoweryBlips},
        {"closemichelhits", "Close Michel hits [#]",  Binning::Simple(40, 0, 80), kMichel_CloseMichelHits},
    };

}
