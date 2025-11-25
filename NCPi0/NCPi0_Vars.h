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

    const Var kNuGraph_FilterFraction([](const caf::SRSliceProxy *slc) -> double {
        if (std::isnan(slc->ng_filt_pass_frac)) return -5.;
        return slc->ng_filt_pass_frac;
    });

    const Var kNuGraph_NShowerPFPs([](const caf::SRSliceProxy *slc) -> int {
        int kNPFPs(0);

        // sem_cat == 2 is for showers!
        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (!std::isnan(slc->reco.pfp[i].ngscore.sem_cat) 
                && (slc->reco.pfp[i].ngscore.sem_cat == 2)
                && (slc->reco.pfp[i].shw.plane[2].nHits > 5))
                kNPFPs += 1;
        }

        return kNPFPs;
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
                // (slc->reco.pfp[i].trk.chi2pid[2].chi2_muon < 30.) &&
                // (slc->reco.pfp[i].trk.chi2pid[2].chi2_proton > 60.) &&
                (slc->reco.pfp[i].ngscore.sem_cat == 0) &&
                (slc->reco.pfp[i].trk.end.x * slc->vertex.x > 0.) &&
                (slc->reco.pfp[i].parent_is_primary)) {
                highestLength = slc->reco.pfp[i].trk.len;
                haveMuonCandidate = true;
            }
        }

        return haveMuonCandidate;
    }); 

    // photon identification
    const MultiVar kOrderedShowerIdxs([](const caf::SRSliceProxy* slc) -> std::vector<double> { 

        std::vector<std::pair<int, double>> showerHitsIdx;  ///< (PFP index, number of best-plane hits)
        int bestPlaneIdx;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            auto const& pfp = slc->reco.pfp[i];

            // shower-like
            if (!(slc->reco.pfp[i].ngscore.sem_cat == 2)) continue;

            if (std::isnan(slc->reco.pfp[i].shw.plane[1].nHits)) { bestPlaneIdx = 2; }
            else if (std::isnan(slc->reco.pfp[i].shw.plane[2].nHits)) { bestPlaneIdx = 1; }
            else {
                bestPlaneIdx = (slc->reco.pfp[i].shw.plane[1].nHits > slc->reco.pfp[i].shw.plane[2].nHits) 
                                ? 1 
                                : 2;
            }

            showerHitsIdx.emplace_back(i, slc->reco.pfp[i].shw.plane[bestPlaneIdx].nHits);
        }

        std::sort(showerHitsIdx.begin(), showerHitsIdx.end(),
                [](auto const& a, auto const& b) {
                    return a.second > b.second; });

        std::vector<double> orderedShowerIdx;
        orderedShowerIdx.reserve(showerHitsIdx.size());
        for (auto const& p : showerHitsIdx)
            orderedShowerIdx.push_back(p.first);

        return orderedShowerIdx;
    });

    const Var kLargestRecoShowerIdx([](const caf::SRSliceProxy* slc) -> int {
        std::vector<double> showerIdxs = kOrderedShowerIdxs(slc);
        return showerIdxs.size() > 0 ? (int) showerIdxs[0] : -1;
    });

    const Var kSubleadingRecoShowerIdx([](const caf::SRSliceProxy* slc) -> int {
        std::vector<double> showerIdxs = kOrderedShowerIdxs(slc);
        return showerIdxs.size() > 0 ? (int) showerIdxs[1] : -1;
    });

    // largest shower properties
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

    const Var kLargestRecoShower_CollEnergy_VsTruth([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;

        double recoEnergy = kLargestRecoShower_CollEnergy(slc);

        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.truth.p.startE) || std::isnan(slc->reco.pfp[largestShwIdx].shw.truth.p.endE)) return -5.;
        double trueEnergy = slc->reco.pfp[largestShwIdx].shw.truth.p.startE - slc->reco.pfp[largestShwIdx].shw.truth.p.endE;

        return (recoEnergy - trueEnergy) / trueEnergy;
    });

    const Var kLargestRecoShower_NuGraph_ShowerFrac([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[largestShwIdx].ngscore.shr_frac)) return -5.;

        return slc->reco.pfp[largestShwIdx].ngscore.shr_frac;
    });

    const Var kLargestRecoShower_NuGraph_HipFrac([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[largestShwIdx].ngscore.hip_frac)) return -5.;

        return slc->reco.pfp[largestShwIdx].ngscore.hip_frac;
    });

    const Var kLargestRecoShower_NuGraph_MipFrac([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[largestShwIdx].ngscore.mip_frac)) return -5.;

        return slc->reco.pfp[largestShwIdx].ngscore.mip_frac;
    });

    const Var kLargestRecoShower_NuGraph_MhlFrac([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[largestShwIdx].ngscore.mhl_frac)) return -5.;

        return slc->reco.pfp[largestShwIdx].ngscore.mhl_frac;
    });

    const Var kLargestRecoShower_NuGraph_DifFrac([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;

        if (std::isnan(slc->reco.pfp[largestShwIdx].ngscore.dif_frac)) return -5.;

        return slc->reco.pfp[largestShwIdx].ngscore.dif_frac;
    });

    // subleading shower properties
    const Var kSubleadRecoShower_CollEnergy([](const caf::SRSliceProxy* slc) -> double {
        const int subleadShwIdx = kSubleadingRecoShowerIdx(slc);
        if (subleadShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[subleadShwIdx].shw.plane[2].energy)) return -5;

        return slc->reco.pfp[subleadShwIdx].shw.plane[2].energy;
    });

    const Var kSubleadRecoShower_TruePdg([](const caf::SRSliceProxy* slc) -> int {
        const int subleadShwIdx = kSubleadingRecoShowerIdx(slc);
        if (subleadShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[subleadShwIdx].shw.truth.p.pdg)) return -5;

        return slc->reco.pfp[subleadShwIdx].shw.truth.p.pdg;
    });

    const Var kSubleadRecoShower_ColldEdx([](const caf::SRSliceProxy* slc) -> double {
        const int subleadShwIdx = kSubleadingRecoShowerIdx(slc);
        if (subleadShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[subleadShwIdx].shw.plane[2].dEdx)) return -5;

        return slc->reco.pfp[subleadShwIdx].shw.plane[2].dEdx;
    });

    const Var kSubleadRecoShower_OpenAngle([](const caf::SRSliceProxy* slc) -> double {
        const int subleadShwIdx = kSubleadingRecoShowerIdx(slc);
        if (subleadShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[subleadShwIdx].shw.open_angle)) return -5;

        return 180. * slc->reco.pfp[subleadShwIdx].shw.open_angle / M_PI;
    });

    const Var kSubleadRecoShower_ConvGap([](const caf::SRSliceProxy* slc) -> double {
        const int subleadShwIdx = kSubleadingRecoShowerIdx(slc);
        if (subleadShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[subleadShwIdx].shw.conversion_gap)) return -5;

        return slc->reco.pfp[subleadShwIdx].shw.conversion_gap;
    });
    const Var kSubleadRecoShower_NuGraph_ShowerFrac([](const caf::SRSliceProxy* slc) -> double {
        const int subleadShwIdx = kSubleadingRecoShowerIdx(slc);
        if (subleadShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[subleadShwIdx].ngscore.shr_frac)) return -5.;

        return slc->reco.pfp[subleadShwIdx].ngscore.shr_frac;
    });

    const Var kSubleadRecoShower_NuGraph_HipFrac([](const caf::SRSliceProxy* slc) -> double {
        const int subleadShwIdx = kSubleadingRecoShowerIdx(slc);
        if (subleadShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[subleadShwIdx].ngscore.hip_frac)) return -5.;

        return slc->reco.pfp[subleadShwIdx].ngscore.hip_frac;
    });

    const Var kSubleadRecoShower_NuGraph_MipFrac([](const caf::SRSliceProxy* slc) -> double {
        const int subleadShwIdx = kSubleadingRecoShowerIdx(slc);
        if (subleadShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[subleadShwIdx].ngscore.mip_frac)) return -5.;

        return slc->reco.pfp[subleadShwIdx].ngscore.mip_frac;
    });

    const Var kSubleadRecoShower_NuGraph_MhlFrac([](const caf::SRSliceProxy* slc) -> double {
        const int subleadShwIdx = kSubleadingRecoShowerIdx(slc);
        if (subleadShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[subleadShwIdx].ngscore.mhl_frac)) return -5.;

        return slc->reco.pfp[subleadShwIdx].ngscore.mhl_frac;
    });

    const Var kSubleadRecoShower_NuGraph_DifFrac([](const caf::SRSliceProxy* slc) -> double {
        const int subleadShwIdx = kSubleadingRecoShowerIdx(slc);
        if (subleadShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[subleadShwIdx].ngscore.dif_frac)) return -5.;

        return slc->reco.pfp[subleadShwIdx].ngscore.dif_frac;
    });

    // Pi0 variables
    const Var kPi0_CosPhotonOpenAngle([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[2].energy)) return -5;
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.dir.x)) return -5;

        const int subleadShwIdx = kSubleadingRecoShowerIdx(slc);
        if (subleadShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[subleadShwIdx].shw.plane[2].energy)) return -5;
        if (std::isnan(slc->reco.pfp[subleadShwIdx].shw.dir.x)) return -5;

        TVector3 largestShwDir(slc->reco.pfp[largestShwIdx].shw.dir.x, 
                            slc->reco.pfp[largestShwIdx].shw.dir.y, 
                            slc->reco.pfp[largestShwIdx].shw.dir.z);
        TVector3 subleadShwDir(slc->reco.pfp[subleadShwIdx].shw.dir.x, 
                            slc->reco.pfp[subleadShwIdx].shw.dir.y, 
                            slc->reco.pfp[subleadShwIdx].shw.dir.z);

        return largestShwDir.Dot(subleadShwDir) / (largestShwDir.Mag() * subleadShwDir.Mag());
    });

    const Var kPi0_InvariantMass([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[2].energy)) return -5;
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.dir.x)) return -5;

        const int subleadShwIdx = kSubleadingRecoShowerIdx(slc);
        if (subleadShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[subleadShwIdx].shw.plane[2].energy)) return -5;
        if (std::isnan(slc->reco.pfp[subleadShwIdx].shw.dir.x)) return -5;

        const double cosPhotonOpenAngle = kPi0_CosPhotonOpenAngle(slc);

        return std::sqrt(2 
                * 1.e3 * slc->reco.pfp[largestShwIdx].shw.plane[2].energy 
                * 1.e3 * slc->reco.pfp[subleadShwIdx].shw.plane[2].energy
                * (1 - cosPhotonOpenAngle));
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
        
        // leading shower variables
        {"collenergy", "#gamma_{1} E_{Coll} [GeV]",                           Binning::Simple(25, 0, 1), kLargestRecoShower_CollEnergy}, 
        {"colldedx", "#gamma_{1} dE/dx_{Coll} [MeV/cm]",                      Binning::Simple(25, 0, 10), kLargestRecoShower_ColldEdx},
        {"openangle", "#gamma_{1} Opening angle [deg.]",                      Binning::Simple(30, 0, 30), kLargestRecoShower_OpenAngle},
        {"convgap", "#gamma_{1} Conversion gap [cm]",                         Binning::Simple(25, 0, 10), kLargestRecoShower_ConvGap},
        
        // subleading shower variables
        {"slcollenergy", "#gamma_{2} E_{Coll} [GeV]",                         Binning::Simple(25, 0, 1), kSubleadRecoShower_CollEnergy}, 
        {"slcolldedx", "#gamma_{2} dE/dx_{Coll} [MeV/cm]",                    Binning::Simple(25, 0, 10), kSubleadRecoShower_ColldEdx},
        {"slopenangle", "#gamma_{2} Opening angle [deg.]",                    Binning::Simple(30, 0, 30), kSubleadRecoShower_OpenAngle},
        {"slconvgap", "#gamma_{2} Conversion gap [cm]",                       Binning::Simple(25, 0, 10), kSubleadRecoShower_ConvGap},
             
        // Pi0 variables
        {"cosphopenangle", "cos(#theta_{#gamma#gamma})",                      Binning::Simple(25, -1, 1), kPi0_CosPhotonOpenAngle},
        {"pi0invmass", "M_{#pi^{0}} [MeV]",                                   Binning::Simple(25, 0, 400), kPi0_InvariantMass},

        // // neutrino variables
        // {"reconuenergy", "E^{reco}_{#nu} [GeV]",                            Binning::Simple(30, 0, 3), kRecoNeutrino_CC0piEnergy},
        // {"nuenergyres", "(E^{reco}_{#nu} - E^{true}_{#nu}) / E^{true}_{#nu}",     Binning::Simple(40, -1, 0.5), kRecoNeutrino_CC0piEnergy_VsTruth},   
        // {"collenergyres", "(E^{reco}_{e} - E^{true}_{e}) / E^{true}_{e}",         Binning::Simple(40, -1, 0.5), kLargestRecoShower_CollEnergy_VsTruth},   
        // {"inelasticity", "y = E^{reco}_{had.} / E^{reco}_{#nu}",                  Binning::Simple(40, 0, 1), kRecoNeutrino_CC0piInelasticity},       
        // {"tranvmomentum", "P_{T} [GeV/c]",                                        Binning::Simple(30, 0, 3), kRecoNeutrino_CC0piTransverseMomentum},     

        // // light information
        // {"barycenterfmdeltaztr", "Barycenter-FM #DeltaZ (trigger) [cm]",    Binning::Simple(40, 0, 150), kBarycenterFM_DeltaZ_Trigger},
        // {"barycenterfmdeltaz", "Barycenter-FM #DeltaZ [cm]",                Binning::Simple(15, 0, 150), kBarycenterFM_DeltaZ},
        // {"barycenterfmtime", "Barycenter-FM time [#mus]",                   Binning::Simple(40, -1, 14), kBarycenterFM_FlashTime},         
        
        // // NuGraph2 variables
        {"ngtaggedshws", "NG2-tagged showers [#]",                                Binning::Simple(7, 0, 7), kNuGraph_NShowerPFPs},
        {"ngshwfrac", "NG2 #gamma_{1} shw_frac",                                  Binning::Simple(25, 0, 1), kLargestRecoShower_NuGraph_ShowerFrac},
        {"nghipfrac", "NG2 #gamma_{1} hip_frac",                                  Binning::Simple(25, 0, 1), kLargestRecoShower_NuGraph_HipFrac},
        {"ngmipfrac", "NG2 #gamma_{1} mip_frac",                                  Binning::Simple(25, 0, 1), kLargestRecoShower_NuGraph_MipFrac},
        {"ngmhlfrac", "NG2 #gamma_{1} mhl_frac",                                  Binning::Simple(25, 0, 1), kLargestRecoShower_NuGraph_MhlFrac},
        {"ngdiffrac", "NG2 #gamma_{1} dif_frac",                                  Binning::Simple(25, 0, 1), kLargestRecoShower_NuGraph_DifFrac}, 
        {"slngshwfrac", "NG2 #gamma_{2} shw_frac",                                Binning::Simple(25, 0, 1), kSubleadRecoShower_NuGraph_ShowerFrac},
        {"slnghipfrac", "NG2 #gamma_{2} hip_frac",                                Binning::Simple(25, 0, 1), kSubleadRecoShower_NuGraph_HipFrac},
        {"slngmipfrac", "NG2 #gamma_{2} mip_frac",                                Binning::Simple(25, 0, 1), kSubleadRecoShower_NuGraph_MipFrac},
        {"slngmhlfrac", "NG2 #gamma_{2} mhl_frac",                                Binning::Simple(25, 0, 1), kSubleadRecoShower_NuGraph_MhlFrac},
        {"slngdiffrac", "NG2 #gamma_{2} dif_frac",                                Binning::Simple(25, 0, 1), kSubleadRecoShower_NuGraph_DifFrac}, 
    };

    // meant for data-MC plots with final selection
    std::vector<PlotDef> SelectionPlots_LowStat = {   
        {"count", "Counts [#]",                                             Binning::Simple(3, 0, 3), kCounting},
        
        // electron variables
        {"collenergy", "E_{Coll} [GeV]",                                    Binning::Simple(25, 0, 2.5), kLargestRecoShower_CollEnergy}, 
        {"colldedx", "dE/dx_{Coll} [MeV/cm]",                               Binning::Simple(40, 0, 11), kLargestRecoShower_ColldEdx},
        {"openangle", "Opening angle [deg.]",                               Binning::Simple(25, 0, 20), kLargestRecoShower_OpenAngle},
        {"convgap", "Conversion gap [cm]",                                  Binning::Simple(25, 0, 8), kLargestRecoShower_ConvGap}, 
        
        // light information
        {"barycenterfmdeltaztr", "Barycenter-FM #DeltaZ (trigger) [cm]",    Binning::Simple(15, 0, 150), kBarycenterFM_DeltaZ_Trigger},
        {"barycenterfmdeltaz", "Barycenter-FM #DeltaZ [cm]",                Binning::Simple(15, 0, 150), kBarycenterFM_DeltaZ},
        {"barycenterfmtime", "Barycenter-FM time [#mus]",                   Binning::Simple(40, -4, 14), kBarycenterFM_FlashTime},         

        // NuGraph2 variables
        {"ngshwfrac", "NG2 e^{#pm} shw_frac",                               Binning::Simple(20, 0, 1), kLargestRecoShower_NuGraph_ShowerFrac},
        {"nghipfrac", "NG2 e^{#pm} hip_frac",                               Binning::Simple(20, 0, 1), kLargestRecoShower_NuGraph_HipFrac},
        {"ngmipfrac", "NG2 e^{#pm} mip_frac",                               Binning::Simple(20, 0, 1), kLargestRecoShower_NuGraph_MipFrac},
        {"ngmhlfrac", "NG2 e^{#pm} mhl_frac",                               Binning::Simple(20, 0, 1), kLargestRecoShower_NuGraph_MhlFrac},
        {"ngdiffrac", "NG2 e^{#pm} dif_frac",                               Binning::Simple(20, 0, 1), kLargestRecoShower_NuGraph_DifFrac},     };

}
