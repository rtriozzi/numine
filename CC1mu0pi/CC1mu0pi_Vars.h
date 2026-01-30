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

    // muon identification
    const Var kMuonIdx([](const caf::SRSliceProxy* slc) -> int {
        int muonIdx(-1);
        double highestLength(-1);

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {

            if (std::isnan(slc->reco.pfp[i].trk.len))
                continue;

            // muon ID 
            if ((slc->reco.pfp[i].trk.len > highestLength) &&
                (slc->reco.pfp[i].ngscore.sem_cat == 0) &&
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
        double K = sqrt(pow(0.10566, 2) + pow(startMomentum.Mag(), 2)); ///< GeV

        return K;
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

    // proton selection
    const MultiVar kNSelectedProtonsIdx([](const caf::SRSliceProxy* slc) -> std::vector<double> { 

        std::vector<double> selectedProtonIdx;
        int NOtherParticles(0);

        const int muonIdx = kMuonIdx(slc);
        if(muonIdx == -1) return selectedProtonIdx;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (i == (unsigned int) muonIdx) continue;

            // MIPs
            if (slc->reco.pfp[i].ngscore.sem_cat == 0) {
                if (kIsPFPPionLike(slc, i)) {
                    NOtherParticles += 1; ///< visible pions
                }
            }
            // HIPs
            else if (slc->reco.pfp[i].ngscore.sem_cat == 1) {
                if (kIsPFPProtonLike(slc, i)) {
                    selectedProtonIdx.push_back(i); ///< visible protons
                }
            }
            // showers
            else if (slc->reco.pfp[i].ngscore.sem_cat == 2) {
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

            // MIPs
            if (slc->reco.pfp[i].ngscore.sem_cat == 0) {
                if (kIsPFPPionLike(slc, i)) {
                    NOtherParticles += 1; ///< visible pions
                }
            }
            // HIPs
            else if (slc->reco.pfp[i].ngscore.sem_cat == 1) {
                if (kIsPFPProtonLike(slc, i)) {
                    selectedProtonIdx.push_back(i); ///< visible protons
                }
            }
            // showers
            else if (slc->reco.pfp[i].ngscore.sem_cat == 2) {
                if (kIsPFPShowerLike(slc, i)) {
                    NOtherParticles += 1; ///< visible shower
                }
            }
        }

        return NOtherParticles;
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

        return E_mu + E_p;
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
