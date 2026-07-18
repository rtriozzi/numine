#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <fstream>
#include <vector>
#include <math.h>
#include "TVector3.h"

#include "CC1mu0pi_Vars.h"   
#include "CC1mu0pi_Cuts.h"
#include "CC1mu0pi_TruthCuts.h"

namespace ana {

    // generic debugger
    const SpillMultiVar kDebugger([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::ofstream myOut("debug/CC1mu0pi_Debug_MC.txt", std::ios::app); 
        std::vector<double> tempSpillVar;
        std::string SourceName = sr->hdr.sourceName;

        // loop over slices
        for (auto const &islc : sr->slc) {
            // automatic event selection
            if (kAutomaticNuMuSelection(&islc)) {

                // header stuff for the selected slice
                myOut << SourceName << "\t" << sr->hdr.run << "\t" << sr->hdr.evt << std::endl;
                myOut << kTrueCC1mu0pi(&islc) << "\t" << kIsNuMu(&islc) << "\t" 
                        << kIsCC(&islc) << "\t" << kTrueVertexInFV(&islc) << std::endl;

                const int muonIdx = kMuonIdx(&islc);
                const double muonKE = kMuon_KE(&islc);
                std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(&islc);
                const double recoNeutrinoEnergy = kRecoNeutrino_NuMuCC0piEnergy(&islc);
                const double trueNeutrinoEnergy = islc.truth.E;

                // NuMuCC other background
                if (!kTrueCC1mu0pi(&islc) && kIsNuMu(&islc) && kIsCC(&islc) && kTrueVertexInFV(&islc)) {
                    myOut << "***********\nTRUE NuMuCC OTHER\n";
                    myOut << "TRUE PARTICLES\n";
                    for (int ip(0); ip < islc.truth.nprim ; ++ip) {
                        myOut << islc.truth.prim[ip].pdg << "\t" << islc.truth.prim[ip].length << "\t";
                        myOut << islc.truth.prim[ip].startE - islc.truth.prim[ip].endE << "\t";
                        myOut << kIsInContained(islc.truth.prim[ip].end.x, islc.truth.prim[ip].end.y, islc.truth.prim[ip].end.z) << "\t";
                        myOut << std::endl;
                    }
                    myOut << "RECO PARTICLES\n";
                    for (unsigned int i = 0; i < islc.reco.npfp; i++) {
                        myOut << islc.reco.pfp[i].trk.truth.p.pdg << "\t" << islc.reco.pfp[i].trk.len << "\t" 
                              << islc.reco.pfp[i].ngscore.sem_cat << "\t";
                        myOut << kIsInContained(islc.reco.pfp[i].trk.end.x, islc.reco.pfp[i].trk.end.y, islc.reco.pfp[i].trk.end.z) << "\t";
                        myOut << std::endl;
                    }
                }
            }
        }
        myOut.close();

        return tempSpillVar;
    });

    // inspect data
    const SpillMultiVar kDataDump([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::ofstream myOut("CC1mu0pi_DataDump.txt", std::ios::app); 
        std::vector<double> tempSpillVar;
        std::string SourceName = sr->hdr.sourceName;

        // loop over slices
        for (auto const &islc : sr->slc) {
            // automatic event selection
            if (kAutomaticNuMuSelection(&islc)) {
                myOut << SourceName << "\t" 
                      << sr->hdr.run << "\t" 
                      << sr->hdr.evt << "\t" 
                      << kMuon_KE(&islc) << "\t" 
                      << kRecoNeutrino_NuMuCC0piEnergy(&islc) << "\t"
                      << kNSelectedProtons_N(&islc) << "\t"
                      << std::endl;
            }
        }
        myOut.close();

        return tempSpillVar;
    });

    // truth dumping
    const SpillMultiVar kEventTruthDump([](const caf::SRSpillProxy* sr) -> std::vector<double>
    {
        
        std::vector<double> tempSpillVar;
        std::string SourceName = sr->hdr.sourceName;

        std::ofstream myOut("debug/NuMI_Truth_Dump.txt", std::ios::app);
        for (auto const& nu : sr->mc.nu) {

            // count stuff in the active volume...
            if (kIsInAV(nu.position.x, nu.position.y, nu.position.z)) {

                // generic neutrino
                myOut << nu.index << "\t" 
                    << nu.pdg << "\t" << nu.iscc << "\t" << nu.genie_mode << "\t" << kIsTrueCC1mu0pi(nu) << "\t"
                    << nu.position.x << "\t" << nu.position.y << "\t" << nu.position.z << "\t"
                    << nu.momentum.x << "\t" << nu.momentum.y << "\t" << nu.momentum.z << "\t"
                    << nu.E << "\t" << nu.baseline << "\t" << nu.xsec << "\t";

                // leading muon, if present
                double LeadingMuonEnergy = -1.;
                double LeadingMuonLength = -1;
                double LeadingMuonPX = -1.;
                double LeadingMuonPY = -1.;
                double LeadingMuonPZ = -1.;
                for (int ip(0); ip < nu.nprim ; ++ip) {
                    if (abs(nu.prim[ip].pdg) == 13) {
                        LeadingMuonEnergy = nu.prim[ip].startE - nu.prim[ip].endE;
                        LeadingMuonLength = nu.prim[ip].length;
                        LeadingMuonPX = nu.prim[ip].startp.x;
                        LeadingMuonPY = nu.prim[ip].startp.y;
                        LeadingMuonPZ = nu.prim[ip].startp.z;
                    }
                }
                myOut << LeadingMuonEnergy << "\t" << LeadingMuonLength << "\t"
                      << LeadingMuonPX << "\t" << LeadingMuonPY << "\t" << LeadingMuonPZ << "\t";
                myOut << std::endl;
            }
        }
        myOut.close();

        return tempSpillVar;
    });

    // count cathode- and z=0-crossing muon candidates among selected slices
    const SpillMultiVar kMuonCrossingDump([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> crossingCode;

        for (auto const &islc : sr->slc) {
            if (!kAutomaticNuMuSelection(&islc)) continue;

            const int muonIdx = kMuonIdx(&islc);
            if (muonIdx == -1) continue; ///< cannot happen post-selection, kept for safety

            double sx = islc.reco.pfp[muonIdx].trk.start.x;
            double ex = islc.reco.pfp[muonIdx].trk.end.x;
            double sz = islc.reco.pfp[muonIdx].trk.start.z;
            double ez = islc.reco.pfp[muonIdx].trk.end.z;
            if (std::isnan(sx) || std::isnan(ex) || std::isnan(sz) || std::isnan(ez)) continue;

            bool crossesCathode = ((fabs(sx) - CATHODE_ABS_X) * (fabs(ex) - CATHODE_ABS_X)) < 0;
            bool crossesZ0 = (sz * ez) < 0;

            // 0: no crossing, 1: all cathode crossers; 2: all z0 crossers, 3: both crossing points
            if (!crossesCathode && !crossesZ0)  crossingCode.push_back(0.);
            if (crossesCathode)                 crossingCode.push_back(1.);
            if (crossesZ0)                      crossingCode.push_back(2.);
            if (crossesCathode && crossesZ0)    crossingCode.push_back(3.);

        }

        return crossingCode;
    });

}