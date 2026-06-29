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
                // // signal investigation
                // if (kTrueCC1mu0pi(&islc)) {
                //     myOut << "***********\nTRUE SIGNAL\n";
                //     myOut << "TRUE PARTICLES\n";
                //     for (int ip(0); ip < islc.truth.nprim ; ++ip) {
                //         myOut << islc.truth.prim[ip].pdg << "\t" << islc.truth.prim[ip].length << "\t";
                //         myOut << islc.truth.prim[ip].startE - islc.truth.prim[ip].endE << "\t";
                //         myOut << std::endl;
                //     }
                //     myOut << "RECO PARTICLES\n";
                //     for (unsigned int i = 0; i < islc.reco.npfp; i++) {
                //         myOut << islc.reco.pfp[i].trk.truth.p.pdg << "\t" << islc.reco.pfp[i].trk.len << "\t" 
                //               << islc.reco.pfp[i].ngscore.sem_cat << "\t";
                //         myOut << std::endl;
                //     }
                //     myOut << "RECO SIGNAL ENERGY\n";
                //     myOut << kRecoNeutrino_CC0piEnergy_VsTruth(&islc) << "\t" 
                //           << kRecoNeutrino_CC0piEnergy(&islc) << "\t" 
                //           << trueNeutrinoEnergy << std::endl;

                //     // muon
                //     TVector3 startMomentumMu(islc.reco.pfp[muonIdx].trk.dir.x * islc.reco.pfp[muonIdx].trk.rangeP.p_muon,
                //                              islc.reco.pfp[muonIdx].trk.dir.y * islc.reco.pfp[muonIdx].trk.rangeP.p_muon, 
                //                              islc.reco.pfp[muonIdx].trk.dir.z * islc.reco.pfp[muonIdx].trk.rangeP.p_muon); 
                //     double K = sqrt(pow(0.10566, 2) + pow(startMomentumMu.Mag(), 2)); ///< GeV
                //     myOut << K << "\t" << (islc.reco.pfp[muonIdx].trk.truth.p.startE - islc.reco.pfp[muonIdx].trk.truth.p.endE) << "\t"
                //           << islc.reco.pfp[muonIdx].trk.truth.p.pdg << "\t";
                //     myOut << std::endl;

                //     // proton
                //     for (auto i : selectedProtonIdx) {
                //         TVector3 startMomentum(islc.reco.pfp[i].trk.dir.x * islc.reco.pfp[i].trk.rangeP.p_proton,
                //                                islc.reco.pfp[i].trk.dir.y * islc.reco.pfp[i].trk.rangeP.p_proton, 
                //                                islc.reco.pfp[i].trk.dir.z * islc.reco.pfp[i].trk.rangeP.p_proton); 
                //         myOut << sqrt(pow(0.9383, 2) + pow(startMomentum.Mag(), 2)) - 0.9383 << "\t" << (islc.reco.pfp[i].trk.truth.p.startE - islc.reco.pfp[i].trk.truth.p.endE) << "\t"
                //               << islc.reco.pfp[i].trk.truth.p.pdg << "\t";
                //         myOut << std::endl;
                //     }
                // }
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

}