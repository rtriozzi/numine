#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <fstream>
#include <vector>
#include <math.h>

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

    const Var kBarycenterFM_DeltaZ([](const caf::SRSliceProxy *slc) {
        return slc->barycenterFM.deltaZ_Trigger;
    });

    const Var kBarycenterFM_FlashTime([](const caf::SRSliceProxy *slc) {
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
        TVector3 startMomentum(slc->reco.pfp[iPFP].trk.dir.x * slc->reco.pfp[iPFP].trk.rangeP.p_proton,
                               slc->reco.pfp[iPFP].trk.dir.y * slc->reco.pfp[iPFP].trk.rangeP.p_proton, 
                               slc->reco.pfp[iPFP].trk.dir.z * slc->reco.pfp[iPFP].trk.rangeP.p_proton); 
        double K = sqrt(pow(0.139570, 2) + pow(startMomentum.Mag(), 2)); ///< GeV

        return (slc->reco.pfp[iPFP].trk.chi2pid[2].chi2_proton > 100) &&
                kIsInContained(slc->reco.pfp[iPFP].trk.end.x, slc->reco.pfp[iPFP].trk.end.y, slc->reco.pfp[iPFP].trk.end.z) &&
                ((recoStart - recoVertex).Mag() < 10) &&
                (K >= VISIBILTY_THRESHOLD);
    }

    // shower identification
    bool kIsPFPShowerLike(const caf::SRSliceProxy* slc, unsigned int iPFP) {
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].shw.start.x) || std::isnan(slc->reco.pfp[iPFP].shw.start.y) || std::isnan(slc->reco.pfp[iPFP].shw.start.z)) return false;

        TVector3 recoVertex(slc->vertex.x, slc->vertex.y, slc->vertex.z); 
        TVector3 recoStart(slc->reco.pfp[iPFP].shw.start.x, slc->reco.pfp[iPFP].shw.start.y, slc->reco.pfp[iPFP].shw.start.z);

        return ((recoStart - recoVertex).Mag() < 10) &&
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

        for (auto i : selectedProtonIdx) { 
            if (std::isnan(slc->reco.pfp[i].trk.dir.x) || std::isnan(slc->reco.pfp[i].trk.dir.y) || std::isnan(slc->reco.pfp[i].trk.dir.z)) return -5;
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5;
            TVector3 startMomentum(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                                   slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton, 
                                   slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton); 
            protonMomenta.push_back(startMomentum.Mag());
        }

        std::sort(protonMomenta.begin(), protonMomenta.end(), std::greater<>());

        if (selectedProtonIdx.size() < 2) return -5.;

        return protonMomenta[1];
    });

    // plotting
    struct PlotDef {
        std::string suffix = "";
        std::string label = "";
        Binning bins = Binning::Simple(3, 0, 3);
        Var var = kCounting;
    };

    std::vector<PlotDef> SelectionPlots = {   
        {"count", "Counts [#]",                                     Binning::Simple(3, 0, 3), kCounting},
        {"barycenterfmdeltaz", "Barycenter-FM #DeltaZ [cm]",        Binning::Simple(40, -1, 200), kBarycenterFM_DeltaZ},
        {"barycenterfmtime", "Barycenter-FM time [#mus]",           Binning::Simple(40, -1, 12), kBarycenterFM_FlashTime},
        {"collenergy", "E_{Coll} [GeV]",                            Binning::Simple(40, -0.1, 3), kLargestRecoShower_CollEnergy},
        {"colldedx", "dE/dx_{Coll} [MeV/cm]",                       Binning::Simple(40, 0, 10), kLargestRecoShower_ColldEdx},
        {"trackscore", "Track score",                               Binning::Simple(70, 0, 1), kLargestRecoShower_TrackScore},
        {"openangle", "Opening angle [deg.]",                       Binning::Simple(40, -1, 30), kLargestRecoShower_OpenAngle},
        {"convgap", "Conversion gap [cm]",                          Binning::Simple(40, -1, 10), kLargestRecoShower_ConvGap},
        {"hitshare", "Hit share",                                   Binning::Simple(40, -0.05, 1), kLargestRecoShower_BestPlaneShowerHitShare},
        {"pca2ratio", "PCA #lambda_{2} / #lambda_{1}",              Binning::Simple(40, -0.02, 0.15), kLargestRecoShower_PCA2Ratio},
        {"pca3ratio", "PCA #lambda_{3} / #lambda_{1}",              Binning::Simple(40, -0.01, 0.05), kLargestRecoShower_PCA3Ratio},
        {"muonrej", "#mu veto",                                     Binning::Simple(2, 0, 2), kHaveMuonCandidate},
        {"leadproton", "P_{p_{1}} [GeV/c]",                         Binning::Simple(40, 0, 2), kLeadingProtonMomentum},
        {"subleadproton", "P_{p_{2}} [GeV/c]",                      Binning::Simple(40, 0, 2), kSubLeadingProtonMomentum},
    };

}