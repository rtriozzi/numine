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
    bool kIsInContained(double ex, double ey, double ez) { 
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

        return slc->reco.pfp[largestShwIdx].shw.plane[2].energy;
    });

    const Var kLargestRecoShower_TruePdg([](const caf::SRSliceProxy* slc) -> int {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return -5;

        return slc->reco.pfp[largestShwIdx].shw.truth.p.pdg;
    });

    const Var kLargestRecoShower_ColldEdx([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1 || slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx < 0) return -5;

        return slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx;
    });

    const Var kLargestRecoShower_TrackScore([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return -5;

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

    // proton identification
    const Var kNProtons_Reconstructed([](const caf::SRSliceProxy* slc) -> int { 

        int NElse(-1);
        int NProtons(-1);
        TVector3 RecoVtx;
        RecoVtx.SetXYZ(slc->vertex.x, slc->vertex.y, slc->vertex.z);
        TVector3 RecoStart;
        TVector3 RecoEnd;
        double minDistToVtx;
        int bestPlaneIdx;
        TVector3 StartMomentum;

        // identify the electron PFP
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if(largestShwIdx == -1) return -1;

        // loop over the PFPs
        for (unsigned int i = 0; i < slc->reco.npfp; i++) {

            // ignore the electron
            if (i == (unsigned int) largestShwIdx)
                continue;

            // PFP has to be primary
            if (!slc->reco.pfp[i].parent_is_primary)
                continue;

            // tracks
            if (slc->reco.pfp[i].trackScore > 0.5) {
                const caf::SRTrackProxy *trk = &slc->reco.pfp[i].trk;
                if (trk) {
                    RecoStart.SetXYZ(trk->start.x, trk->start.y, trk->start.z);
                    RecoEnd.SetXYZ(trk->end.x, trk->end.y, trk->end.z);

                    // accounts for flipped tracks
                    minDistToVtx = (RecoStart-RecoVtx).Mag() < (RecoEnd-RecoVtx).Mag() ? (RecoStart-RecoVtx).Mag() : (RecoEnd-RecoVtx).Mag();
                    if (minDistToVtx < 50 && 
                        kIsInContained(trk->end.x, trk->end.y, trk->end.z)) {

                        // pions
                        if (trk->chi2pid[2].chi2_proton >= 100) {
                            StartMomentum.SetXYZ(trk->rangeP.p_pion*trk->dir.x, trk->rangeP.p_pion*trk->dir.y, trk->rangeP.p_pion*trk->dir.z);
                            
                            if ( (sqrt(pow(139.570,2) + pow(StartMomentum.Mag()*1.e3, 2))-139.570) >= VISIBILTY_THRESHOLD*1.e3 )
                                ++NElse;
                        }

                        // protons
                        if (trk->chi2pid[2].chi2_proton < 100) {
                            StartMomentum.SetXYZ(trk->rangeP.p_proton*trk->dir.x, trk->rangeP.p_proton*trk->dir.y, trk->rangeP.p_proton*trk->dir.z);

                            if ((sqrt(pow(938.3,2) + pow(StartMomentum.Mag()*1.e3, 2))-938.3) >= VISIBILTY_THRESHOLD*1.e3 &&
                                (minDistToVtx < 10) &&
                                kIsInContained(trk->end.x, trk->end.y, trk->end.z))
                                    ++NProtons;
                        }
                    }
                }
            } // end of track case

            // showers
            if (slc->reco.pfp[i].trackScore <= 0.5) {
                const caf::SRShowerProxy *shw = &slc->reco.pfp[i].shw;
                if (shw) {
                    RecoStart.SetXYZ(shw->start.x, shw->start.y, shw->start.z);
                    RecoEnd.SetXYZ(shw->end.x, shw->end.y, shw->end.z);
                    minDistToVtx = (RecoStart-RecoVtx).Mag() < (RecoEnd-RecoVtx).Mag() ? (RecoStart-RecoVtx).Mag() : (RecoEnd-RecoVtx).Mag();

                    if (minDistToVtx < 50 &&
                        shw->bestplane_energy > VISIBILTY_THRESHOLD*1.e3) 
                            ++NElse;
                }
            } // end of shower case

        } // end of loop over the PFPs

        return NProtons;

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
        {"pca2ratio", "PCA #lambda_2 / #lambda_1",                  Binning::Simple(40, -0.02, 0.15), kLargestRecoShower_PCA2Ratio},
        {"pca3ratio", "PCA #lambda_3 / #lambda_1",                  Binning::Simple(40, -0.01, 0.05), kLargestRecoShower_PCA3Ratio},
        {"muonrej", "#mu veto",                                     Binning::Simple(2, 0, 2), kHaveMuonCandidate},
    };

}