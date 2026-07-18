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
#include "TRandom3.h"

namespace ana {

    // selection thresholds
    const double VISIBILTY_THRESHOLD_P = 0.05;
    const double VISIBILTY_THRESHOLD_PI = 0.025;
    const double TRACK_SCORE_CUT = 0.5;
    const double PROTON_TRACK_SCORE_CUT = 0.45;
    const double BG_SHW_MAX_DISTANCE = 50;

    // track breaking at the cathode and z=0
    const double CATHODE_ABS_X = 210.215;  
    const double Z0_BREAK_OFFSET = 1.;

    const double N_SELECTED_MC = 9421.;
    const double CATHODE_N_CROSSING_MC = 2959.;
    const double Z0_N_CROSSING_MC = 821.;

    const double CATHODE_ENDX_DENSITY_DATA = 0.009;
    const double CATHODE_ENDX_DENSITY_MC = 0.003;
    const double Z0_ENDZ_DENSITY_DATA = 0.005;
    const double Z0_ENDZ_DENSITY_MC = 0.001;

    const double ENDPOINT_BIN_WIDTH = 2.;
    const double CATHODE_BREAK_PROB = (CATHODE_ENDX_DENSITY_DATA - CATHODE_ENDX_DENSITY_MC) * ENDPOINT_BIN_WIDTH * N_SELECTED_MC / CATHODE_N_CROSSING_MC;
    const double Z0_BREAK_PROB = (Z0_ENDZ_DENSITY_DATA - Z0_ENDZ_DENSITY_MC) * ENDPOINT_BIN_WIDTH * N_SELECTED_MC / Z0_N_CROSSING_MC;

    // general helper functions
    bool kIsInAV(double x, double y, double z) {  
        if (std::isnan(x) || std::isnan(y) || std::isnan(z)) return false;

        return (( ( x < -61.94 && x > -358.49 ) ||
            ( x >  61.94 && x <  358.49 )) &&
            ( ( y > -181.86 && y < 134.96 ) &&
            ( z > -894.95 && z < 894.95 ) ));
    }

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

    // muon identification
    const Var kMuonIdx([](const caf::SRSliceProxy* slc) -> int {
        int muonIdx(-1);
        double highestLength(-1);
        TVector3 recoVertex(slc->vertex.x, slc->vertex.y, slc->vertex.z); 

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {

            if (std::isnan(slc->reco.pfp[i].trk.len))
                continue;

            // check start position compared to vertex
            TVector3 recoStart(slc->reco.pfp[i].trk.start.x, slc->reco.pfp[i].trk.start.y, slc->reco.pfp[i].trk.start.z);

            // muon ID 
            if ((slc->reco.pfp[i].trk.len > highestLength) &&
                ((recoStart - recoVertex).Mag() < 10) &&
                // (slc->reco.pfp[i].ngscore.sem_cat == 0) &&
                (slc->reco.pfp[i].trk.chi2pid[2].chi2_muon < 30) && (slc->reco.pfp[i].trk.chi2pid[2].chi2_proton > 60) && 
                (slc->reco.pfp[i].trackScore >= TRACK_SCORE_CUT) &&
                (slc->reco.pfp[i].parent_is_primary)) {
                muonIdx = i;
                highestLength = slc->reco.pfp[i].trk.len;
            }
        }

        return muonIdx;
    });

    // muon properties
    const Var kMuon_Length([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if(muonIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[muonIdx].trk.len)) return -5;

        return slc->reco.pfp[muonIdx].trk.len;
    });

    const Var kMuon_TrackScore([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if(muonIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[muonIdx].trackScore)) return -5;

        return slc->reco.pfp[muonIdx].trackScore;
    });

    const Var kMuon_TrueLength([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if(muonIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[muonIdx].trk.truth.p.length)) return -5;

        return slc->reco.pfp[muonIdx].trk.truth.p.length;
    });

    const Var kMuon_EndX([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if(muonIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[muonIdx].trk.end.x)) return -5;

        return slc->reco.pfp[muonIdx].trk.end.x;
    });

    const Var kMuon_EndY([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if(muonIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[muonIdx].trk.end.y)) return -5;

        return slc->reco.pfp[muonIdx].trk.end.y;
    });

    const Var kMuon_EndZ([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if(muonIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[muonIdx].trk.end.z)) return -5;

        return slc->reco.pfp[muonIdx].trk.end.z;
    });

    const Var kMuon_Chi2Muon([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if(muonIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[muonIdx].trk.chi2pid[2].chi2_muon)) return -5;

        return slc->reco.pfp[muonIdx].trk.chi2pid[2].chi2_muon;
    });

    const Var kMuon_Chi2Proton([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if(muonIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[muonIdx].trk.chi2pid[2].chi2_proton)) return -5;

        return slc->reco.pfp[muonIdx].trk.chi2pid[2].chi2_proton;
    });

    const Var kMuon_Length_VsTruth([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[muonIdx].trk.len) || std::isnan(slc->reco.pfp[muonIdx].trk.truth.p.length)) 
            return -5;

        double lengthResidual = (slc->reco.pfp[muonIdx].trk.truth.p.length > 0) 
                                ? (slc->reco.pfp[muonIdx].trk.len - slc->reco.pfp[muonIdx].trk.truth.p.length) / slc->reco.pfp[muonIdx].trk.truth.p.length
                                : -5;
        return lengthResidual;
    });

    const Var kMuon_Momentum([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if(muonIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[muonIdx].trk.rangeP.p_muon)) return -5;

        return slc->reco.pfp[muonIdx].trk.rangeP.p_muon;
    });

    const Var kMuon_KE([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if(muonIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[muonIdx].trk.len)) return -5;
        if (std::isnan(slc->reco.pfp[muonIdx].trk.dir.x) || std::isnan(slc->reco.pfp[muonIdx].trk.dir.y) || std::isnan(slc->reco.pfp[muonIdx].trk.dir.z)) return -5;

        TVector3 startMomentum(slc->reco.pfp[muonIdx].trk.dir.x * slc->reco.pfp[muonIdx].trk.rangeP.p_muon,
                               slc->reco.pfp[muonIdx].trk.dir.y * slc->reco.pfp[muonIdx].trk.rangeP.p_muon, 
                               slc->reco.pfp[muonIdx].trk.dir.z * slc->reco.pfp[muonIdx].trk.rangeP.p_muon); 
        double K = sqrt(pow(0.10566, 2) + pow(startMomentum.Mag(), 2)) - 0.10566; ///< GeV

        return K;
    });

    const Var kMuon_TrueKE([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if(muonIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[muonIdx].trk.truth.p.startE) || std::isnan(slc->reco.pfp[muonIdx].trk.truth.p.endE)) return -5;

        return slc->reco.pfp[muonIdx].trk.truth.p.startE - slc->reco.pfp[muonIdx].trk.truth.p.endE;
    });

    const Var kMuon_KE_VsTruth([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[muonIdx].trk.len)) return -5;
        if (std::isnan(slc->reco.pfp[muonIdx].trk.truth.p.startE) || std::isnan(slc->reco.pfp[muonIdx].trk.truth.p.endE)) return -5;
        if (std::isnan(slc->reco.pfp[muonIdx].trk.dir.x) || std::isnan(slc->reco.pfp[muonIdx].trk.dir.y) || std::isnan(slc->reco.pfp[muonIdx].trk.dir.z)) return -5;

        TVector3 startMomentum(slc->reco.pfp[muonIdx].trk.dir.x * slc->reco.pfp[muonIdx].trk.rangeP.p_muon,
                               slc->reco.pfp[muonIdx].trk.dir.y * slc->reco.pfp[muonIdx].trk.rangeP.p_muon, 
                               slc->reco.pfp[muonIdx].trk.dir.z * slc->reco.pfp[muonIdx].trk.rangeP.p_muon); 
        double K = sqrt(pow(0.10566, 2) + pow(startMomentum.Mag(), 2)); ///< GeV
        double trueK = slc->reco.pfp[muonIdx].trk.truth.p.startE - slc->reco.pfp[muonIdx].trk.truth.p.endE;
        return (trueK > 0)
               ? (K - trueK) / trueK
               : -5;
    });

    const Var kMuon_NuGraph_MIPFrac([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[muonIdx].ngscore.mip_frac)) return -5.;

        return slc->reco.pfp[muonIdx].ngscore.mip_frac;
    });

    const Var kMuon_NuGraph_MhlFrac([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[muonIdx].ngscore.mhl_frac)) return -5.;

        return slc->reco.pfp[muonIdx].ngscore.mhl_frac;
    });

    // pion identification (was: MIP / sem_cat == 0)
    bool kIsPFPPionLike(const caf::SRSliceProxy* slc, unsigned int iPFP) {
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].trk.start.x) || std::isnan(slc->reco.pfp[iPFP].trk.start.y) || std::isnan(slc->reco.pfp[iPFP].trk.start.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].trk.end.x) || std::isnan(slc->reco.pfp[iPFP].trk.end.y) || std::isnan(slc->reco.pfp[iPFP].trk.end.z)) return false;

        TVector3 recoVertex(slc->vertex.x, slc->vertex.y, slc->vertex.z); 
        TVector3 recoStart(slc->reco.pfp[iPFP].trk.start.x, slc->reco.pfp[iPFP].trk.start.y, slc->reco.pfp[iPFP].trk.start.z);
        TVector3 startMomentum(slc->reco.pfp[iPFP].trk.dir.x * slc->reco.pfp[iPFP].trk.rangeP.p_pion,
                            slc->reco.pfp[iPFP].trk.dir.y * slc->reco.pfp[iPFP].trk.rangeP.p_pion, 
                            slc->reco.pfp[iPFP].trk.dir.z * slc->reco.pfp[iPFP].trk.rangeP.p_pion); 
        double K = sqrt(pow(0.139570, 2) + pow(startMomentum.Mag(), 2)) - 0.139570; ///< GeV

        return ((recoStart - recoVertex).Mag() < 10) &&
            (K >= VISIBILTY_THRESHOLD_PI) &&
            (slc->reco.pfp[iPFP].trk.chi2pid[2].chi2_muon < 30) &&
            (slc->reco.pfp[iPFP].trk.chi2pid[2].chi2_proton > 60);
    }

    // generic shower identification
    bool kIsPFPShowerLike(const caf::SRSliceProxy* slc, unsigned int iPFP) {
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].shw.start.x) || std::isnan(slc->reco.pfp[iPFP].shw.start.y) || std::isnan(slc->reco.pfp[iPFP].shw.start.z)) return false;

        TVector3 recoVertex(slc->vertex.x, slc->vertex.y, slc->vertex.z); 
        TVector3 recoStart(slc->reco.pfp[iPFP].shw.start.x, slc->reco.pfp[iPFP].shw.start.y, slc->reco.pfp[iPFP].shw.start.z);

        return ((recoStart - recoVertex).Mag() < BG_SHW_MAX_DISTANCE) &&
               (slc->reco.pfp[iPFP].shw.plane[2].energy >= VISIBILTY_THRESHOLD_PI);
    }

    // proton identification (was: HIP / sem_cat == 1)
    bool kIsPFPProtonLike(const caf::SRSliceProxy* slc, unsigned int iPFP) {
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].trk.start.x) || std::isnan(slc->reco.pfp[iPFP].trk.start.y) || std::isnan(slc->reco.pfp[iPFP].trk.start.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].trk.end.x) || std::isnan(slc->reco.pfp[iPFP].trk.end.y) || std::isnan(slc->reco.pfp[iPFP].trk.end.z)) return false;

        TVector3 recoVertex(slc->vertex.x, slc->vertex.y, slc->vertex.z); 
        TVector3 recoStart(slc->reco.pfp[iPFP].trk.start.x, slc->reco.pfp[iPFP].trk.start.y, slc->reco.pfp[iPFP].trk.start.z);
        TVector3 startMomentum(slc->reco.pfp[iPFP].trk.dir.x * slc->reco.pfp[iPFP].trk.rangeP.p_proton,
                            slc->reco.pfp[iPFP].trk.dir.y * slc->reco.pfp[iPFP].trk.rangeP.p_proton, 
                            slc->reco.pfp[iPFP].trk.dir.z * slc->reco.pfp[iPFP].trk.rangeP.p_proton); 
        double K = sqrt(pow(0.9383, 2) + pow(startMomentum.Mag(), 2)) - 0.938272; ///< GeV

        return kIsInContained(slc->reco.pfp[iPFP].trk.end.x, slc->reco.pfp[iPFP].trk.end.y, slc->reco.pfp[iPFP].trk.end.z) &&
            ((recoStart - recoVertex).Mag() < 10) &&
            (slc->reco.pfp[iPFP].parent_is_primary) &&
            (K >= VISIBILTY_THRESHOLD_P) &&
            (slc->reco.pfp[iPFP].trk.chi2pid[2].chi2_proton < 90) &&
            (slc->reco.pfp[iPFP].trk.chi2pid[2].chi2_muon > 30);
    }

    // helper to check for vertex-PFP start proximity
    bool kIsNearVertex(const caf::SRSliceProxy* slc, unsigned int iPFP) {
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].trk.start.x) || std::isnan(slc->reco.pfp[iPFP].trk.start.y) || std::isnan(slc->reco.pfp[iPFP].trk.start.z)) return false;

        TVector3 recoVertex(slc->vertex.x, slc->vertex.y, slc->vertex.z);
        TVector3 recoStart(slc->reco.pfp[iPFP].trk.start.x, slc->reco.pfp[iPFP].trk.start.y, slc->reco.pfp[iPFP].trk.start.z);

        return (recoStart - recoVertex).Mag() < 10;
    }

    // proton selection
    const MultiVar kNSelectedProtonsIdx([](const caf::SRSliceProxy* slc) -> std::vector<double> { 

        std::vector<double> selectedProtonIdx;
        int NOtherParticles(0);

        const int muonIdx = kMuonIdx(slc);
        if(muonIdx == -1) return selectedProtonIdx;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (i == (unsigned int) muonIdx) continue;

            // track-like (was: sem_cat == 0 or 1)
            if (slc->reco.pfp[i].trackScore >= PROTON_TRACK_SCORE_CUT) {
                // exiting track from the vertex: unambiguous background
                if (!kIsInContained(slc->reco.pfp[i].trk.end.x, slc->reco.pfp[i].trk.end.y, slc->reco.pfp[i].trk.end.z) && kIsNearVertex(slc, i)) {
                    NOtherParticles += 1; ///< exiting track
                }
                else if (kIsPFPProtonLike(slc, i)) {
                    selectedProtonIdx.push_back(i); ///< visible protons
                }
                else if (kIsPFPPionLike(slc, i)) {
                    NOtherParticles += 1; ///< visible pions
                }
                else if (slc->reco.pfp[i].trackScore < TRACK_SCORE_CUT) {
                    if (kIsPFPShowerLike(slc, i)) {
                        NOtherParticles += 1; ///< visible shower
                    } 
                }
            }
            // shower-like (was: sem_cat == 2)
            else {
                if (kIsPFPShowerLike(slc, i)) {
                    NOtherParticles += 1; ///< visible shower
                }
            }
        }

        return selectedProtonIdx;
    });

    // complementary var to proton selection
    const Var kNSelectedProtonsIdx_NOtherParticles([](const caf::SRSliceProxy* slc) -> int { 

        std::vector<double> selectedProtonIdx;
        int NOtherParticles(0);

        const int muonIdx = kMuonIdx(slc);
        if(muonIdx == -1) return -1;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (i == (unsigned int) muonIdx) continue;

            if (slc->reco.pfp[i].trackScore >= PROTON_TRACK_SCORE_CUT) {
                // exiting track from the vertex: unambiguous background
                if (!kIsInContained(slc->reco.pfp[i].trk.end.x, slc->reco.pfp[i].trk.end.y, slc->reco.pfp[i].trk.end.z) && kIsNearVertex(slc, i)) {
                    NOtherParticles += 1; ///< exiting track
                }
                else if (kIsPFPProtonLike(slc, i)) {
                    selectedProtonIdx.push_back(i);
                }
                else if (kIsPFPPionLike(slc, i)) {
                    NOtherParticles += 1;
                }
                else if (slc->reco.pfp[i].trackScore < TRACK_SCORE_CUT) {
                    if (kIsPFPShowerLike(slc, i)) {
                        NOtherParticles += 1; ///< visible shower
                    } 
                }
            }
            else {
                if (kIsPFPShowerLike(slc, i)) {
                    NOtherParticles += 1;
                }
            }
        }

        return NOtherParticles;
    });

    const Var kNSelectedProtons_N([](const caf::SRSliceProxy* slc) -> double { 
    
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);

        return selectedProtonIdx.size();
    });

    // proton properties
    const Var kLeadingProtonMomentum([](const caf::SRSliceProxy* slc) -> double { 
    
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        std::vector<double> protonMomenta;

        if (selectedProtonIdx.empty()) return -5.;

        for (auto i : selectedProtonIdx) { 
            if (std::isnan(slc->reco.pfp[i].trk.dir.x) || std::isnan(slc->reco.pfp[i].trk.dir.y) || std::isnan(slc->reco.pfp[i].trk.dir.z)) return -5;
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5;
            TVector3 startMomentum(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                                   slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton, 
                                   slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton); 
            protonMomenta.push_back(startMomentum.Mag());
        }

        std::sort(protonMomenta.begin(), protonMomenta.end(), std::greater<>());

        return protonMomenta[0];
    });

    const Var kSubLeadingProtonMomentum([](const caf::SRSliceProxy* slc) -> double { 
    
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        std::vector<double> protonMomenta;

        if (selectedProtonIdx.empty()) return -5.;
        if (selectedProtonIdx.size() < 2) return -5.;

        for (auto i : selectedProtonIdx) { 
            if (std::isnan(slc->reco.pfp[i].trk.dir.x) || std::isnan(slc->reco.pfp[i].trk.dir.y) || std::isnan(slc->reco.pfp[i].trk.dir.z)) return -5;
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5;
            TVector3 startMomentum(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                                   slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton, 
                                   slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton); 
            protonMomenta.push_back(startMomentum.Mag());
        }

        std::sort(protonMomenta.begin(), protonMomenta.end(), std::greater<>());

        return protonMomenta[1];
    });

    const Var kLeadingProton_Chi2Muon([](const caf::SRSliceProxy* slc) -> double { 
    
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        std::vector<double> protonMomenta;
        std::vector<double> Chi2;

        if (selectedProtonIdx.empty()) return -5.;

        for (auto i : selectedProtonIdx) { 
            if (std::isnan(slc->reco.pfp[i].trk.dir.x) || std::isnan(slc->reco.pfp[i].trk.dir.y) || std::isnan(slc->reco.pfp[i].trk.dir.z)) return -5;
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5;
            TVector3 startMomentum(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                                   slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton, 
                                   slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton); 
            protonMomenta.push_back(startMomentum.Mag());
            Chi2.push_back(slc->reco.pfp[i].trk.chi2pid[2].chi2_muon); 
        }

        std::vector<unsigned int> idx(protonMomenta.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(),
              [&](double i1, double i2) {return protonMomenta[i1] > protonMomenta[i2];});

        return Chi2[idx[0]];
    });

    const Var kLeadingProton_Chi2Proton([](const caf::SRSliceProxy* slc) -> double { 
    
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        std::vector<double> protonMomenta;
        std::vector<double> Chi2;

        if (selectedProtonIdx.empty()) return -5.;

        for (auto i : selectedProtonIdx) { 
            if (std::isnan(slc->reco.pfp[i].trk.dir.x) || std::isnan(slc->reco.pfp[i].trk.dir.y) || std::isnan(slc->reco.pfp[i].trk.dir.z)) return -5;
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5;
            TVector3 startMomentum(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                                   slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton, 
                                   slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton); 
            protonMomenta.push_back(startMomentum.Mag());
            Chi2.push_back(slc->reco.pfp[i].trk.chi2pid[2].chi2_proton); 
        }

        std::vector<unsigned int> idx(protonMomenta.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(),
              [&](double i1, double i2) {return protonMomenta[i1] > protonMomenta[i2];});

        return Chi2[idx[0]];
    });

    const Var kLeadingProton_TrackScore([](const caf::SRSliceProxy* slc) -> double { 
    
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        std::vector<double> protonMomenta;
        std::vector<double> trackScores;

        if (selectedProtonIdx.empty()) return -5.;

        for (auto i : selectedProtonIdx) { 
            if (std::isnan(slc->reco.pfp[i].trk.dir.x) || std::isnan(slc->reco.pfp[i].trk.dir.y) || std::isnan(slc->reco.pfp[i].trk.dir.z)) return -5;
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5;
            TVector3 startMomentum(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                                   slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton, 
                                   slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton); 
            protonMomenta.push_back(startMomentum.Mag());
            trackScores.push_back(slc->reco.pfp[i].trackScore); 
        }

        std::vector<unsigned int> idx(protonMomenta.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(),
              [&](double i1, double i2) {return protonMomenta[i1] > protonMomenta[i2];});

        return trackScores[idx[0]];
    });

    const Var kLeadingProton_NuGraph_HipFrac([](const caf::SRSliceProxy* slc) -> double { 
    
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        std::vector<double> protonMomenta;
        std::vector<double> HIPFracs;

        if (selectedProtonIdx.empty()) return -5.;

        for (auto i : selectedProtonIdx) { 
            if (std::isnan(slc->reco.pfp[i].trk.dir.x) || std::isnan(slc->reco.pfp[i].trk.dir.y) || std::isnan(slc->reco.pfp[i].trk.dir.z)) return -5;
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5;
            TVector3 startMomentum(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                                   slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton, 
                                   slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton); 
            protonMomenta.push_back(startMomentum.Mag());
            HIPFracs.push_back(slc->reco.pfp[i].ngscore.hip_frac);
        }

        std::vector<unsigned int> idx(protonMomenta.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(),
              [&](double i1, double i2) {return protonMomenta[i1] > protonMomenta[i2];});

        return HIPFracs[idx[0]];
    });

    const Var kLeadingProton_NuGraph_MipFrac([](const caf::SRSliceProxy* slc) -> double { 
    
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        std::vector<double> protonMomenta;
        std::vector<double> MIPFracs;

        if (selectedProtonIdx.empty()) return -5.;

        for (auto i : selectedProtonIdx) { 
            if (std::isnan(slc->reco.pfp[i].trk.dir.x) || std::isnan(slc->reco.pfp[i].trk.dir.y) || std::isnan(slc->reco.pfp[i].trk.dir.z)) return -5;
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5;
            TVector3 startMomentum(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                                   slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton, 
                                   slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton); 
            protonMomenta.push_back(startMomentum.Mag());
            MIPFracs.push_back(slc->reco.pfp[i].ngscore.mip_frac);
        }

        std::vector<unsigned int> idx(protonMomenta.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(),
              [&](double i1, double i2) {return protonMomenta[i1] > protonMomenta[i2];});

        return MIPFracs[idx[0]];
    });

    // neutrino properties
    const Var kRecoNeutrino_NuMuCC0piEnergy([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx == -1) return -5.;

        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;

        if (kMuon_KE(slc) < 0) return -5.;
        double E_mu = kMuon_KE(slc);

        double E_p = 0.;
        for (auto i : selectedProtonIdx) { 
            if (std::isnan(slc->reco.pfp[i].trk.dir.x) || std::isnan(slc->reco.pfp[i].trk.dir.y) || std::isnan(slc->reco.pfp[i].trk.dir.z)) return -5.;
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5.;
            TVector3 startMomentum(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                                   slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton, 
                                   slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton); 
            E_p += sqrt(pow(0.9383, 2) + pow(startMomentum.Mag(), 2)) - 0.9383 + 0.0309;
        }

        return E_mu + 0.10566 + E_p;
    });

    const Var kRecoNeutrino_NuMuCC0piEnergy_VsTruth([](const caf::SRSliceProxy* slc) -> double {
        double recoNeutrinoEnergy = kRecoNeutrino_NuMuCC0piEnergy(slc);
        if (recoNeutrinoEnergy < 0) return -5.;

        if (std::isnan(slc->truth.E)) return -5.;
        double trueNeutrinoEnergy = slc->truth.E;

        if (recoNeutrinoEnergy < 0 || std::isnan(slc->truth.E) || trueNeutrinoEnergy < 0) {
            return -5;
        }
        else {
            return (recoNeutrinoEnergy - trueNeutrinoEnergy) / trueNeutrinoEnergy;
        }
    });

    const Var kRecoNeutrino_NuMuCC0piTransverseMomentum([](const caf::SRSliceProxy* slc) -> double {
        // muon
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[muonIdx].trk.len)) return -5.;
        TVector3 startMomentum(slc->reco.pfp[muonIdx].trk.dir.x * slc->reco.pfp[muonIdx].trk.rangeP.p_muon,
                               slc->reco.pfp[muonIdx].trk.dir.y * slc->reco.pfp[muonIdx].trk.rangeP.p_muon, 
                               slc->reco.pfp[muonIdx].trk.dir.z * slc->reco.pfp[muonIdx].trk.rangeP.p_muon); 

        // protons
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;

        TVector3 startMomentumP(0., 0., 0.);
        for (auto i : selectedProtonIdx) { 
            if (std::isnan(slc->reco.pfp[i].trk.dir.x) || std::isnan(slc->reco.pfp[i].trk.dir.y) || std::isnan(slc->reco.pfp[i].trk.dir.z)) return -5.;
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5.;
            startMomentumP.SetXYZ(
                startMomentumP.X() + slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                startMomentumP.Y() + slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                startMomentumP.Z() + slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton
            );
        }

        TVector3 totalMomentum = startMomentum + startMomentumP;

        return sqrt(pow(totalMomentum.X(), 2) + pow(totalMomentum.Y(), 2));
    });

    const Var kRecoNeutrino_NuMuCC0piTransverseMomentum_NuMI([](const caf::SRSliceProxy* slc) -> double {
        // muon
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[muonIdx].trk.len)) return -5.;
        TVector3 startMomentum(slc->reco.pfp[muonIdx].trk.dir.x * slc->reco.pfp[muonIdx].trk.rangeP.p_muon,
                               slc->reco.pfp[muonIdx].trk.dir.y * slc->reco.pfp[muonIdx].trk.rangeP.p_muon, 
                               slc->reco.pfp[muonIdx].trk.dir.z * slc->reco.pfp[muonIdx].trk.rangeP.p_muon); 

        // protons
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;

        TVector3 startMomentumP(0., 0., 0.);
        for (auto i : selectedProtonIdx) { 
            if (std::isnan(slc->reco.pfp[i].trk.dir.x) || std::isnan(slc->reco.pfp[i].trk.dir.y) || std::isnan(slc->reco.pfp[i].trk.dir.z)) return -5.;
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5.;
            startMomentumP.SetXYZ(
                startMomentumP.X() + slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                startMomentumP.Y() + slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                startMomentumP.Z() + slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton
            );
        }

        TVector3 totalMomentum = startMomentum + startMomentumP;

        // NuMI beam direction unit vector (already normalised)
        const TVector3 numiDir(3.94583e-01, 4.26067e-02, 9.17677e-01);

        // Longitudinal component along beam: p_L = (p · n_hat) n_hat
        TVector3 pLongitudinal = totalMomentum.Dot(numiDir) * numiDir;

        // Transverse component: p_T = p - p_L
        TVector3 pTransverse = totalMomentum - pLongitudinal;

        return pTransverse.Mag();
    });

    // track breaking systematic knob creation

    // fixed seed for reproducible aggregate behaviour
    double kBreakUniform() {
        static TRandom3 rndmBreak(42);
        return rndmBreak.Rndm();
    }

    // muon momentum [GeV/c] from CSDA range, as larreco TrackMomentumCalculator (PDG table),
    // we need this to recompute the spectrum after track breaking
    double MuonMomentumFromRange(double len) {
        static const std::vector<double> rangeGcm2 = {
            9.833e-1, 1.786e0, 3.321e0, 6.598e0, 1.058e1, 3.084e1, 4.250e1, 6.732e1,
            1.063e2, 1.725e2, 2.385e2, 4.934e2, 6.163e2, 8.552e2, 1.202e3, 1.758e3,
            2.297e3, 4.359e3, 5.354e3, 7.298e3, 1.013e4, 1.469e4, 1.910e4 };
        static const std::vector<double> keMeV = {
            10., 14., 20., 30., 40., 80., 100., 140.,
            200., 300., 400., 800., 1000., 1400., 2000., 3000.,
            4000., 8000., 10000., 14000., 20000., 30000., 40000. };

        double r = len * 1.396; ///< g/cm2, LAr density as in larreco
        if (r <= rangeGcm2.front()) return 0.;
        if (r >= rangeGcm2.back()) return -5.;

        auto it = std::upper_bound(rangeGcm2.begin(), rangeGcm2.end(), r);
        size_t i = std::distance(rangeGcm2.begin(), it);
        double KE = keMeV[i-1] + (keMeV[i] - keMeV[i-1]) * (r - rangeGcm2[i-1]) / (rangeGcm2[i] - rangeGcm2[i-1]);

        return sqrt(KE * (KE + 2. * 105.658)) / 1000.;
    }

    // track breaking at z = 0
    const Var kRecoNeutrino_NuMuCC0piEnergy_BrokenAtZ0([](const caf::SRSliceProxy* slc) -> double {
        double E_nom = kRecoNeutrino_NuMuCC0piEnergy(slc);
        if (E_nom < 0) return -5.;

        const int muonIdx = kMuonIdx(slc);

        double sz = slc->reco.pfp[muonIdx].trk.start.z;
        double ez = slc->reco.pfp[muonIdx].trk.end.z;
        if (std::isnan(sz) || std::isnan(ez)) return -5.;

        // no z=0 crossing: energy unchanged
        if (sz * ez > 0) return E_nom;

        // break decision first: unbroken tracks are untouched
        if (kBreakUniform() >= Z0_BREAK_PROB) return E_nom;

        // stub ends just before the gap, on the entry side
        double breakZ = (sz < 0) ? -Z0_BREAK_OFFSET : Z0_BREAK_OFFSET;

        double newLen = slc->reco.pfp[muonIdx].trk.len * fabs(sz - breakZ) / fabs(ez - sz);
        if (newLen < 50.) return -5.; ///< broken track leaves the selection

        double newP = MuonMomentumFromRange(newLen);
        if (newP <= 0) return -5.;
        double newKE = sqrt(pow(0.10566, 2) + pow(newP, 2)) - 0.10566;

        return E_nom - kMuon_KE(slc) + newKE;
    });

    const Var kMuon_EndZ_BrokenAtZ0([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx == -1) return -9999.;

        double sz = slc->reco.pfp[muonIdx].trk.start.z;
        double ez = slc->reco.pfp[muonIdx].trk.end.z;
        if (std::isnan(sz) || std::isnan(ez)) return -9999.;

        if (sz * ez > 0) return ez;

        if (kBreakUniform() >= Z0_BREAK_PROB) return ez;

        double breakZ = (sz < 0) ? -Z0_BREAK_OFFSET : Z0_BREAK_OFFSET;

        double newLen = slc->reco.pfp[muonIdx].trk.len * fabs(sz - breakZ) / fabs(ez - sz);
        if (newLen < 50.) return -9999.;

        return breakZ;
    });

    // track breaking at cathodes
    const Var kRecoNeutrino_NuMuCC0piEnergy_BrokenAtCathode([](const caf::SRSliceProxy* slc) -> double {
        double E_nom = kRecoNeutrino_NuMuCC0piEnergy(slc);
        if (E_nom < 0) return -5.;

        const int muonIdx = kMuonIdx(slc);

        double sx = slc->reco.pfp[muonIdx].trk.start.x;
        double ex = slc->reco.pfp[muonIdx].trk.end.x;
        if (std::isnan(sx) || std::isnan(ex)) return -5.;

        // no cathode crossing: energy unchanged
        if ((fabs(sx) - CATHODE_ABS_X) * (fabs(ex) - CATHODE_ABS_X) > 0) return E_nom;

        // break decision first: unbroken tracks are untouched
        if (kBreakUniform() >= CATHODE_BREAK_PROB) return E_nom;

        double newLen = slc->reco.pfp[muonIdx].trk.len * fabs(fabs(sx) - CATHODE_ABS_X) / fabs(fabs(ex) - fabs(sx));
        if (newLen < 50.) return -5.; ///< broken track leaves the selection

        double newP = MuonMomentumFromRange(newLen);
        if (newP <= 0.) return -5.;
        double newKE = sqrt(pow(0.10566, 2) + pow(newP, 2)) - 0.10566;

        return E_nom - kMuon_KE(slc) + newKE;
    });

    const Var kMuon_EndX_BrokenAtCathode([](const caf::SRSliceProxy* slc) -> double {
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx == -1) return -9999.;

        double sx = slc->reco.pfp[muonIdx].trk.start.x;
        double ex = slc->reco.pfp[muonIdx].trk.end.x;
        if (std::isnan(sx) || std::isnan(ex)) return -9999.;

        if ((fabs(sx) - CATHODE_ABS_X) * (fabs(ex) - CATHODE_ABS_X) > 0) return ex;

        if (kBreakUniform() >= CATHODE_BREAK_PROB) return ex;

        double newLen = slc->reco.pfp[muonIdx].trk.len * fabs(fabs(sx) - CATHODE_ABS_X) / fabs(fabs(ex) - fabs(sx));
        if (newLen < 50.) return -9999.;

        // new endpoint on the cathode, signed by the cryostat side
        return (sx < 0) ? -CATHODE_ABS_X : CATHODE_ABS_X;
    });

    // plotting
    struct PlotDef {
        std::string suffix = "";
        std::string label = "";
        Binning bins = Binning::Simple(3, 0, 3);
        Var var = kCounting;
    };

    std::vector<PlotDef> NuMuSelectionPlots = {   
        {"count", "Counts [#]",                                             Binning::Simple(3, 0, 3), kCounting},
        
        // muon variables
        {"mulen", "Muon length [cm]",                                       Binning::Simple(30, 0, 600), kMuon_Length}, 
        {"muke", "Muon KE [GeV]",                                           Binning::Simple(30, 0, 3), kMuon_KE}, 
        {"mup", "Muon P [GeV]",                                             Binning::Simple(30, 0, 3), kMuon_Momentum}, 
        {"mulenres", "(L^{reco}_{#mu} - L^{true}_{#mu}) / L^{true}_{#mu}",  Binning::Simple(50, -1, 1), kMuon_Length_VsTruth}, 
        {"mukeres", "(KE^{reco}_{#mu} - KE^{true}_{#mu}) / KE^{true}_{#mu}",Binning::Simple(50, -1, 1), kMuon_KE_VsTruth}, 
        {"ngmipfrac", "Muon mip_frac",                                      Binning::Simple(25, 0, 1), kMuon_NuGraph_MIPFrac}, 
        {"ngmhlfrac", "Muon mhl_frac",                                      Binning::Simple(25, 0, 1), kMuon_NuGraph_MhlFrac}, 
  
        // proton variables
        {"leadproton", "P_{p_{1}} [GeV/c]",                                 Binning::Simple(30, 0, 2), kLeadingProtonMomentum},
        {"subleadproton", "P_{p_{2}} [GeV/c]",                              Binning::Simple(30, 0, 2), kSubLeadingProtonMomentum},
        {"slphipfrac", "P_{1} hip_frac",                                    Binning::Simple(25, 0, 1), kLeadingProton_NuGraph_HipFrac},
        {"slpmipfrac", "P_{1} mip_frac",                                    Binning::Simple(25, 0, 1), kLeadingProton_NuGraph_MipFrac},

        // neutrino variables
        {"reconuenergy", "E^{reco}_{#nu} [GeV]",                                  Binning::Simple(30, 0, 3), kRecoNeutrino_NuMuCC0piEnergy},
        {"nuenergyres", "(E^{reco}_{#nu} - E^{true}_{#nu}) / E^{true}_{#nu}",     Binning::Simple(50, -1, 1), kRecoNeutrino_NuMuCC0piEnergy_VsTruth},   
        {"tranvmomentum", "P_{T} [GeV/c]",                                        Binning::Simple(30, 0, 3), kRecoNeutrino_NuMuCC0piTransverseMomentum},     

        // light information
        {"barycenterfmdeltaztr", "Barycenter-FM #DeltaZ (trigger) [cm]",    Binning::Simple(40, 0, 150), kBarycenterFM_DeltaZ_Trigger},
        {"barycenterfmdeltaz", "Barycenter-FM #DeltaZ [cm]",                Binning::Simple(15, 0, 150), kBarycenterFM_DeltaZ},
        {"barycenterfmtime", "Barycenter-FM time [#mus]",                   Binning::Simple(40, -1, 14), kBarycenterFM_FlashTime},         
        
        // NuGraph2 variables
        {"shrhitsind1", "Ind1 shr_hits [#]",                                Binning::Simple(50, 0, 100), kNuGraph_Ind1ShowerHits}, 
        {"shrhitsind2", "Ind2 shr_hits [#]",                                Binning::Simple(50, 0, 100), kNuGraph_Ind2ShowerHits}, 
        {"shrhitscoll", "Coll shr_hits [#]",                                Binning::Simple(50, 0, 100), kNuGraph_CollShowerHits}, 
        {"unclshrhitsind1", "Ind1 unclustered_shr_hits [#]",                Binning::Simple(25, 0, 50), kNuGraph_Ind1ShowerHits_Unclustered}, 
        {"unclshrhitsind2", "Ind2 unclustered_shr_hits [#]",                Binning::Simple(25, 0, 50), kNuGraph_Ind2ShowerHits_Unclustered}, 
        {"unclshrhitscoll", "Coll unclustered_shr_hits [#]",                Binning::Simple(25, 0, 50), kNuGraph_CollShowerHits_Unclustered}, 
    };
}
