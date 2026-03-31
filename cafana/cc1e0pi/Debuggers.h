#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <fstream>
#include <vector>
#include <math.h>
#include "TVector3.h"

#include "CC1e0piSelection_Vars.h"   
#include "CC1e0piSelection_Cuts.h"
#include "CC1e0piSelection_TruthCuts.h"

namespace ana {

    // event dumping
    const SpillMultiVar kEventDump([](const caf::SRSpillProxy* sr) -> std::vector<double>
    {
        
        std::vector<double> tempSpillVar;
        std::string SourceName = sr->hdr.sourceName;

        std::ofstream myOut("NuMI_Prescaled_CC1e0piSelection.txt", std::ios::app);
        for (auto const &islc : sr->slc) {
            if (kAutomaticSelection(&islc)) {
                // reco information
                myOut << sr->hdr.run << "\t" << sr->hdr.evt << "\t" << SourceName << "\t"
                      << islc.vertex.x << "\t" << kLargestRecoShower_CollEnergy(&islc) << "\t" 
                      << kRecoNeutrino_CC0piEnergy(&islc) << "\t";
                myOut << std::endl;
            }
        }
        myOut.close();

        return tempSpillVar;
    });

    // event dumping
    const SpillMultiVar kEventDump_PreSelection([](const caf::SRSpillProxy* sr) -> std::vector<double>
    {
        
        std::vector<double> tempSpillVar;
        std::string SourceName = sr->hdr.sourceName;

        std::ofstream myOut("NuMI_Prescaled_ShowerPreSelection.txt", std::ios::app);
        for (auto const &islc : sr->slc) {
            if (kPreSelection(&islc)) {
                // reco information
                myOut << sr->hdr.run << "\t" << sr->hdr.evt << "\t" << SourceName << "\t"
                      << islc.vertex.x << "\t" << kLargestRecoShower_CollEnergy(&islc) << "\t" 
                      << kRecoNeutrino_CC0piEnergy(&islc) << "\t" << kLargestRecoShower_OpenAngle(&islc) << "\t" 
                      << kLargestRecoShower_ConvGap(&islc) << "\t";
                myOut << std::endl;
            }
        }
        myOut.close();

        return tempSpillVar;
    });

    const SpillMultiVar kMCEventDump_PreSelection([](const caf::SRSpillProxy* sr) -> std::vector<double>
    {
        
        std::vector<double> tempSpillVar;
        std::string SourceName = sr->hdr.sourceName;
        int NPi0;

        std::ofstream myOut("NuMI_MC_ShowerPreSelection.txt", std::ios::app);
        for (auto const &islc : sr->slc) {

            NPi0 = 0;

            if (kPreSelection(&islc)) {
                // reco information
                myOut << sr->hdr.run << "\t" << sr->hdr.evt << "\t" << SourceName << "\t"
                      << islc.vertex.x << "\t" << kLargestRecoShower_CollEnergy(&islc) << "\t" 
                      << kRecoNeutrino_CC0piEnergy(&islc) << "\t" << kLargestRecoShower_OpenAngle(&islc) << "\t" 
                      << kLargestRecoShower_ConvGap(&islc) << "\t";
                // truth information (discard for data)
                myOut << islc.truth.index << "\t" << islc.truth.pdg << "\t" 
                      << islc.truth.iscc << "\t";
                // Pi0
                for (int ip(0); ip < islc.truth.nprim ; ++ip) {
                    if (islc.truth.prim[ip].pdg == 111) {
                        ++NPi0;
                    }
                }
                myOut << NPi0 << "\t";
                myOut << std::endl;
            }
        }
        myOut.close();

        return tempSpillVar;
    });

    const SpillMultiVar kDebugSelection([](const caf::SRSpillProxy* sr) -> std::vector<double> {

        bool trueNeutrinoWasCounted;
        std::vector<double> selectedEnergies;
        std::string SourceName = sr->hdr.sourceName;

        // std::cout << sr->hdr.run << "\t" << sr->hdr.evt << std::endl;
        // std::cout << "Entering spill with " << sr->mc.nnu << " neutrinos." << std::endl;

        for (auto const& nu : sr->mc.nu) {
            if (kIsTrueCC1e0pi(nu)) {
                trueNeutrinoWasCounted = false;

                std::cout << SourceName << std::endl;
                std::cout << sr->hdr.run << "\t" << sr->hdr.evt << std::endl;
                std::cout << "True neutrino with index " << nu.index << std::endl;
                std::cout << "Positions: " << nu.position.x << "\t" << nu.position.y << "\t" << nu.position.z << std::endl;

                for (auto const& islc : sr->slc) {
                    std::cout << " ->Slice with index " << islc.truth.index << std::endl;
                    std::cout << " ->Reconstructed slice is signal? " << kTrueCC1e0pi(&islc) << std::endl;
                    //std::cout << " ->Reconstructed slice passes cut? " << cut(&islc) << std::endl;
                    std::cout << " ->True neutrino was already counted? " << trueNeutrinoWasCounted << std::endl;
                    std::cout << " ->Reco vertex: " << islc.vertex.x << "\t" << islc.vertex.y << "\t" << islc.vertex.z << std::endl;
                    if ((islc.truth.index == nu.index) &&  ///< same index, to account for pile-up
                        (kTrueCC1e0pi(&islc)) &&           ///< signal slice
                        //(cut(&islc)) &&                    ///< apply reconstruction cut step
                        (!trueNeutrinoWasCounted)) {       ///< make sure to count only one slice per true neutrino

                        selectedEnergies.push_back(islc.truth.E);
                        trueNeutrinoWasCounted = true;
                    }
                }
            }
        }

        return selectedEnergies;
    });

    const SpillMultiVar kMCEventDump_DebugSelection([](const caf::SRSpillProxy* sr) -> std::vector<double>
    {
        std::vector<double> tempSpillVar;
        std::string SourceName = sr->hdr.sourceName;

        std::ofstream myOut("debug/numinue_CV.txt", std::ios::app);
        // std::ofstream myOutSlim("debug/numinue_StandardReco_Slim.txt", std::ios::app);     

        for (auto const &islc : sr->slc) {
            if (kTrueCC1e0pi(&islc)) {
                myOut << SourceName << "\t" << sr->hdr.run << "\t" << sr->hdr.evt << "\t" << islc.tmatch.eff << std::endl;
                myOut << "Properties..." << std::endl;
                myOut << kLargestRecoShower_TrueEnergy(&islc) << "\t" << kLargestRecoShower_CollEnergy(&islc) << "\t"
                      << kLargestRecoShower_TrueLength(&islc) << "\t" << kLargestRecoShower_Length(&islc) << "\t"
                      << kLargestRecoShower_EndZ(&islc) << std::endl;
                myOut << "Cuts..." << std::endl;
                myOut << "NCC: " << "\t" << kNotClearCosmic(&islc) << std::endl;
                myOut << "Vtx in FV: " << "\t" << kVertexInFV(&islc) << std::endl;
                myOut << "Trig. flash: " << "\t" << kTrigFlashMatch(&islc) << std::endl;
                myOut << "Shw > 200 MeV: " << "\t" << kLargestRecoShower_EnergyCut(&islc) << std::endl;
                myOut << "dEdx cut: " << "\t" << kLargestRecoShower_dEdxCut(&islc) << std::endl;
                myOut << "angle cut: " << "\t" << kLargestRecoShower_OpenAngleCut(&islc) << std::endl;
                myOut << "gap cut: " << "\t" << kLargestRecoShower_ConvGapCut(&islc) << std::endl;
                myOut << "protons: " << "\t" << kNSelectedProtons(&islc) << std::endl;
                myOut << "muon veto: " << "\t" << kMuonVeto(&islc) << std::endl;

            }
        }
        myOut.close();

        return tempSpillVar;
    });   

}