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

    const double VISIBILTY_THRESHOLD = 0.05;

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

    const Var kBarycenterFM_DeltaZ([](const caf::SRSliceProxy *slc) -> double {
        return slc->barycenterFM.deltaZ;
    });

    const Var kBarycenterFM_DeltaZ_Trigger([](const caf::SRSliceProxy *slc) -> double {
        return slc->barycenterFM.deltaZ_Trigger;
    });

    const Var kBarycenterFM_FlashTime([](const caf::SRSliceProxy *slc) -> double {
        return slc->barycenterFM.flashTime;
    });

    // muon rejection
    const Var kHaveMuonCandidate([](const caf::SRSliceProxy *slc) -> bool {
        bool haveMuonCandidate = false;
        double highestLength = 0.;
        
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
        TVector3 recoVertex(slc->vertex.x, slc->vertex.y, slc->vertex.z);
        TVector3 recoStart;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (std::isnan(slc->reco.pfp[i].trk.start.x) || std::isnan(slc->reco.pfp[i].trk.len)) {
                continue;
            }
            recoStart.SetXYZ(slc->reco.pfp[i].trk.start.x, slc->reco.pfp[i].trk.start.y, slc->reco.pfp[i].trk.start.z);

            if (slc->reco.pfp[i].trackScore < 0.5) {
                continue;
            }

            if ((slc->reco.pfp[i].trk.len > highestLength) &&
                // (slc->reco.pfp[i].trk.len > 50) &&
                (slc->reco.pfp[i].trk.chi2pid[2].chi2_muon < 30.) &&
                (slc->reco.pfp[i].trk.chi2pid[2].chi2_proton > 60.) &&
                // ((recoVertex - recoStart).Mag() < 10) &&
                // (kIsInContained(slc->reco.pfp[i].trk.end.x, slc->reco.pfp[i].trk.end.y, slc->reco.pfp[i].trk.end.z)) &&
                (slc->reco.pfp[i].trk.end.x * slc->vertex.x > 0.) &&
                (slc->reco.pfp[i].parent_is_primary)) {
                highestLength = slc->reco.pfp[i].trk.len;
                haveMuonCandidate = true;
            }
        }

        return haveMuonCandidate;
    }); 

    // electron identification
    const Var kLargestRecoShowerIdx([](const caf::SRSliceProxy* slc) -> int {
        int electronIdx(-1);
        double maxNHits(-1);

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {

            // get the hits on the best-plane, Collection or Induction-2
            auto const& shw = slc->reco.pfp[i].shw;
            int bestPlaneIdx = shw.plane[1].nHits > shw.plane[2].nHits ? 1 : 2;

            // electron ID
            if ((shw.plane[bestPlaneIdx].nHits > maxNHits) &&
                (slc->reco.pfp[i].trackScore <= 0.5) &&
                (slc->reco.pfp[i].parent_is_primary)) {
                electronIdx = i;
                maxNHits = shw.plane[bestPlaneIdx].nHits;
            }
        }

        return electronIdx;
    });

    // electron properties
    const Var kLargestRecoShower_CollEnergy([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[2].energy)) return -5;

        return slc->reco.pfp[largestShwIdx].shw.plane[2].energy;
    });

    const Var kLargestRecoShower_TruePdg([](const caf::SRSliceProxy* slc) -> int {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[largestShwIdx].shw.truth.p.pdg)) return -5;

        return slc->reco.pfp[largestShwIdx].shw.truth.p.pdg;
    });

    const Var kLargestRecoShower_ColldEdx([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx)) return -5;

        return slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx;
    });

    const Var kLargestRecoShower_AvailabledEdx([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return -5;

        if(!std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx) && (slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx > 0)) {   
            return slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx;
        }
        else if (!std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[1].dEdx) && (slc->reco.pfp[largestShwIdx].shw.plane[1].dEdx > 0)) {
            return slc->reco.pfp[largestShwIdx].shw.plane[1].dEdx;
        }
        else if (!std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[0].dEdx) && (slc->reco.pfp[largestShwIdx].shw.plane[0].dEdx > 0)) {
            return slc->reco.pfp[largestShwIdx].shw.plane[0].dEdx;
        }
        else {
            return -5.;
        }
    });

    const Var kLargestRecoShower_TrackScore([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return -5;
        if(std::isnan(slc->reco.pfp[largestShwIdx].trackScore)) return -5;

        return slc->reco.pfp[largestShwIdx].trackScore;
    });

    const Var kLargestRecoShower_OpenAngle([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return -5;

        return 180. * slc->reco.pfp[largestShwIdx].shw.open_angle / M_PI;
    });

    const Var kLargestRecoShower_ConvGap([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return -5;

        return slc->reco.pfp[largestShwIdx].shw.conversion_gap;
    });

    const Var kLargestRecoShower_PCA2Ratio([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return -5;

        return slc->reco.pfp[largestShwIdx].pfochar.pca2ratio;
    });

    const Var kLargestRecoShower_PCA3Ratio([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return -5;

        return slc->reco.pfp[largestShwIdx].pfochar.pca3ratio;
    });

    const Var kLargestRecoShower_BestPlaneShowerHitShare([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return -5;
        int NShowerHits = 0;
        int bestPlaneIdx = 2;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            bestPlaneIdx = slc->reco.pfp[i].shw.plane[2].nHits > slc->reco.pfp[i].shw.plane[1].nHits ? 2 : 1;
            NShowerHits += slc->reco.pfp[i].shw.plane[bestPlaneIdx].nHits;
        }

        if (NShowerHits > 0) {
            return (float) slc->reco.pfp[largestShwIdx].shw.plane[bestPlaneIdx].nHits / NShowerHits;
        }
        else {
            return -5;
        }
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

        return (slc->reco.pfp[iPFP].trk.chi2pid[2].chi2_proton > 100) &&
                kIsInContained(slc->reco.pfp[iPFP].trk.end.x, slc->reco.pfp[iPFP].trk.end.y, slc->reco.pfp[iPFP].trk.end.z) &&
                ((recoStart - recoVertex).Mag() < 10) &&
                (K >= VISIBILTY_THRESHOLD);
    }

    // generic shower identification
    bool kIsPFPShowerLike(const caf::SRSliceProxy* slc, unsigned int iPFP) {
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].shw.start.x) || std::isnan(slc->reco.pfp[iPFP].shw.start.y) || std::isnan(slc->reco.pfp[iPFP].shw.start.z)) return false;

        TVector3 recoVertex(slc->vertex.x, slc->vertex.y, slc->vertex.z); 
        TVector3 recoStart(slc->reco.pfp[iPFP].shw.start.x, slc->reco.pfp[iPFP].shw.start.y, slc->reco.pfp[iPFP].shw.start.z);

        return ((recoStart - recoVertex).Mag() < 50) &&
               (slc->reco.pfp[iPFP].shw.plane[2].energy >= VISIBILTY_THRESHOLD);
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

        return (slc->reco.pfp[iPFP].trk.chi2pid[2].chi2_proton <= 100) &&
                kIsInContained(slc->reco.pfp[iPFP].trk.end.x, slc->reco.pfp[iPFP].trk.end.y, slc->reco.pfp[iPFP].trk.end.z) &&
                ((recoStart - recoVertex).Mag() < 10) &&
                (K >= VISIBILTY_THRESHOLD);
    }

    // proton selection
    const MultiVar kNSelectedProtonsIdx([](const caf::SRSliceProxy* slc) -> std::vector<double> { 

        std::vector<double> selectedProtonIdx;
        int NOtherParticles(0);

        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return selectedProtonIdx;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (i == (unsigned int) largestShwIdx) continue;

            if (slc->reco.pfp[i].trackScore > 0.5) {
                if (kIsPFPPionLike(slc, i)) {
                    NOtherParticles += 1; ///< visible pions
                }
                else if (kIsPFPProtonLike(slc, i)) {
                    selectedProtonIdx.push_back(i); ///< visible protons
                }
            }
            else {
                if (slc->reco.pfp[i].trackScore > 0.45) {
                    if (kIsPFPProtonLike(slc, i)) {
                        selectedProtonIdx.push_back(i); ///< shower-like visible protons
                    }
                    else if (kIsPFPShowerLike(slc, i)) {
                        NOtherParticles += 1; ///< visible shower
                    }
                }
                else {
                    if (kIsPFPShowerLike(slc, i)) {
                        NOtherParticles += 1; ///< visible shower
                    }
                }
            }
        }

        return selectedProtonIdx;
    });

    // complementary var to proton selection
    const Var kNSelectedProtonsIdx_NOtherParticles([](const caf::SRSliceProxy* slc) -> int { 

        std::vector<int> selectedProtonIdx;
        int NOtherParticles(0);

        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return -1;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (i == (unsigned int) largestShwIdx) continue;

            if (slc->reco.pfp[i].trackScore > 0.5) {
                if (kIsPFPPionLike(slc, i)) {
                    NOtherParticles += 1; ///< visible pions
                }
                else if (kIsPFPProtonLike(slc, i)) {
                    selectedProtonIdx.push_back(i); ///< visible protons
                }
            }
            else {
                if (slc->reco.pfp[i].trackScore > 0.45) {
                    if (kIsPFPProtonLike(slc, i)) {
                        selectedProtonIdx.push_back(i); ///< shower-like visible protons
                    }
                    else if (kIsPFPShowerLike(slc, i)) {
                        NOtherParticles += 1; ///< visible shower
                    }
                }
                else {
                    if (kIsPFPShowerLike(slc, i)) {
                        NOtherParticles += 1; ///< visible shower
                    }
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

    const Var kSubLeadingProton_TrackScore([](const caf::SRSliceProxy* slc) -> double { 
    
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        std::vector<double> protonMomenta;
        std::vector<double> trackScores;

        if (selectedProtonIdx.empty()) return -5.;
        if (selectedProtonIdx.size() < 2) return -5.;

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

        return trackScores[idx[1]];
    });

    // neutrino properties
    const Var kRecoNeutrino_CC0piEnergy([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;

        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;

        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[2].energy)) return -5.;
        double E_e = slc->reco.pfp[largestShwIdx].shw.plane[2].energy;

        double E_p = 0.;
        for (auto i : selectedProtonIdx) { 
            if (std::isnan(slc->reco.pfp[i].trk.dir.x) || std::isnan(slc->reco.pfp[i].trk.dir.y) || std::isnan(slc->reco.pfp[i].trk.dir.z)) return -5.;
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5.;
            TVector3 startMomentum(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                                   slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton, 
                                   slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton); 
            E_p += sqrt(pow(0.9383, 2) + pow(startMomentum.Mag(), 2)) - 0.9383;
        }

        return E_e + E_p;
    });

    const Var kRecoNeutrino_CC0piEnergy_VsTruth([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;

        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;

        double recoNeutrinoEnergy = kRecoNeutrino_CC0piEnergy(slc);

        if (std::isnan(slc->truth.E)) return -5.;
        double trueNeutrinoEnergy = slc->truth.E;

        return (recoNeutrinoEnergy - trueNeutrinoEnergy) / trueNeutrinoEnergy;
    });

    const Var kRecoNeutrino_CC0piHadronEnergy([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;

        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;

        double E_p = 0.;
        for (auto i : selectedProtonIdx) { 
            if (std::isnan(slc->reco.pfp[i].trk.dir.x) || std::isnan(slc->reco.pfp[i].trk.dir.y) || std::isnan(slc->reco.pfp[i].trk.dir.z)) return -5.;
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5.;
            TVector3 startMomentum(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                                   slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton, 
                                   slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton); 
            E_p += sqrt(pow(0.9383, 2) + pow(startMomentum.Mag(), 2)) - 0.9383;
        }

        return E_p;
    });

    const Var kRecoNeutrino_CC0piInelasticity([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;

        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;

        return kRecoNeutrino_CC0piHadronEnergy(slc) / kRecoNeutrino_CC0piEnergy(slc);
    });

    const Var kRecoNeutrino_CC0piTransverseMomentum([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;

        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;

        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[2].energy)) return -5.;
        double P_e = sqrt(pow(slc->reco.pfp[largestShwIdx].shw.plane[2].energy, 2) - pow(0.510998e-3, 2));
        TVector3 startMomentumE(P_e * slc->reco.pfp[largestShwIdx].shw.dir.x,
                                P_e * slc->reco.pfp[largestShwIdx].shw.dir.y,
                                P_e * slc->reco.pfp[largestShwIdx].shw.dir.z
        );

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

        TVector3 totalMomentum = startMomentumE + startMomentumP;

        return sqrt(pow(totalMomentum.X(), 2) + pow(totalMomentum.Y(), 2));
    });

    // plotting
    struct PlotDef {
        std::string suffix = "";
        std::string label = "";
        Binning bins = Binning::Simple(3, 0, 3);
        Var var = kCounting;
    };

    std::vector<PlotDef> SelectionPlots = {   
        {"count", "Counts [#]",                                             Binning::Simple(3, 0, 3), kCounting},
        {"barycenterfmdeltaztr", "Barycenter-FM #DeltaZ (trigger) [cm]",    Binning::Simple(40, 0, 150), kBarycenterFM_DeltaZ_Trigger},
        {"barycenterfmdeltaz", "Barycenter-FM #DeltaZ [cm]",                Binning::Simple(15, 0, 150), kBarycenterFM_DeltaZ},
        {"barycenterfmtime", "Barycenter-FM time [#mus]",                   Binning::Simple(40, -1, 14), kBarycenterFM_FlashTime},
        {"collenergy", "E_{Coll} [GeV]",                                    Binning::Simple(40, 0, 3), kLargestRecoShower_CollEnergy},
        {"colldedx", "dE/dx_{Coll} [MeV/cm]",                               Binning::Simple(40, 0, 9), kLargestRecoShower_ColldEdx},
        {"availdedx", "dE/dx_{Coll, Ind} [MeV/cm]",                         Binning::Simple(40, 0, 9), kLargestRecoShower_AvailabledEdx},
        {"trackscore", "Track score",                                       Binning::Simple(50, 0, 1), kLargestRecoShower_TrackScore},
        {"openangle", "Opening angle [deg.]",                               Binning::Simple(40, 0, 30), kLargestRecoShower_OpenAngle},
        {"convgap", "Conversion gap [cm]",                                  Binning::Simple(40, 0, 10), kLargestRecoShower_ConvGap},
        {"hitshare", "Hit share",                                           Binning::Simple(40, 0, 1), kLargestRecoShower_BestPlaneShowerHitShare},
        {"muonrej", "#mu veto",                                             Binning::Simple(2, 0, 2), kHaveMuonCandidate},
        {"leadproton", "P_{p_{1}} [GeV/c]",                                 Binning::Simple(30, 0, 2), kLeadingProtonMomentum},
        {"subleadproton", "P_{p_{2}} [GeV/c]",                              Binning::Simple(30, 0, 2), kSubLeadingProtonMomentum},
        {"leadpts", "Track score (p_{1})",                                  Binning::Simple(50, 0, 1), kLeadingProton_TrackScore},
        {"lsubeadpts", "Track score (p_{2})",                               Binning::Simple(50, 0, 1), kSubLeadingProton_TrackScore},
        {"reconuenergy", "E^{reco}_{#nu} [GeV]",                            Binning::Simple(30, 0, 3), kRecoNeutrino_CC0piEnergy},
        {"nuenergyres", "(E^{reco}_{#nu} - E^{true}_{#nu}) / E^{true}_{#nu}",     Binning::Simple(40, -1, 0.5), kRecoNeutrino_CC0piEnergy_VsTruth},    
        {"inelasticity", "y = E^{reco}_{had.} / E^{reco}_{#nu}",                  Binning::Simple(40, 0, 1), kRecoNeutrino_CC0piInelasticity},       
        {"tranvmomentum", "P_{T} [GeV/c]",                                        Binning::Simple(30, 0, 3), kRecoNeutrino_CC0piTransverseMomentum},       
    };

    // meant for data-MC plots with final selection
    std::vector<PlotDef> SelectionPlots_LowStat = {   
        {"count", "Counts [#]",                                             Binning::Simple(3, 0, 3), kCounting},
        {"barycenterfmdeltaztr", "Barycenter-FM #DeltaZ (trigger) [cm]",    Binning::Simple(20, 0, 140), kBarycenterFM_DeltaZ_Trigger},
        {"collenergy", "E_{Coll} [GeV]",                                    Binning::Simple(25, 0, 2.5), kLargestRecoShower_CollEnergy},
        {"colldedx", "dE/dx_{Coll} [MeV/cm]",                               Binning::Simple(30, 0, 7), kLargestRecoShower_ColldEdx},
        {"availdedx", "dE/dx_{Coll, Ind} [MeV/cm]",                         Binning::Simple(20, 0, 7), kLargestRecoShower_AvailabledEdx},
        {"trackscore", "Track score",                                       Binning::Simple(40, 0, 1), kLargestRecoShower_TrackScore},
        {"openangle", "Opening angle [deg.]",                               Binning::Simple(25, 0, 20), kLargestRecoShower_OpenAngle},
        {"convgap", "Conversion gap [cm]",                                  Binning::Simple(25, 0, 8), kLargestRecoShower_ConvGap},
        {"leadproton", "P_{p_{1}} [GeV/c]",                                 Binning::Simple(15, 0, 1.5), kLeadingProtonMomentum},
        {"subleadproton", "P_{p_{2}} [GeV/c]",                              Binning::Simple(15, 0, 1.5), kSubLeadingProtonMomentum},
        {"reconuenergy", "E^{reco}_{#nu} [GeV]",                            Binning::Simple(20, 0, 3), kRecoNeutrino_CC0piEnergy},
        {"inelasticity", "y = E^{reco}_{had.} / E^{reco}_{#nu}",            Binning::Simple(20, 0, 1), kRecoNeutrino_CC0piInelasticity},     
        {"tranvmomentum", "P_{T} [GeV/c]",                                  Binning::Simple(20, 0, 2), kRecoNeutrino_CC0piTransverseMomentum},       
    };

}
