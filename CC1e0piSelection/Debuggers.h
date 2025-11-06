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

    const SpillMultiVar kMCEventDump_DebugSelection([](const caf::SRSpillProxy* sr) -> std::vector<double>
    {
        std::vector<double> tempSpillVar;
        std::string SourceName = sr->hdr.sourceName;

        std::ofstream myOut("NuMI_MC_DebugSelection.txt", std::ios::app);
        for (auto const &islc : sr->slc) {
            if (islc.truth.index >= 0) {
                myOut << SourceName << "\t" << sr->hdr.run << "\t" << sr->hdr.evt << std::endl;
                myOut << "Neutrino slice of PDG " << islc.truth.pdg << " and CC/NC " << islc.truth.iscc << std::endl;
                
                // myOut << "Primaries [pdg, 2.visE, startE, endE, startE-endE]" << std::endl;
                // int vCryo = islc.truth.position.x < 0 ? 0 : 1;
                // for (int ip(0); ip < islc.truth.nprim ; ++ip) {
                //     myOut << islc.truth.prim[ip].pdg << "\t" << islc.truth.prim[ip].plane[vCryo][2].visE << "\t"
                //           << islc.truth.prim[ip].startE << "\t" << islc.truth.prim[ip].endE << "\t"
                //           << islc.truth.prim[ip].startE - islc.truth.prim[ip].endE << std::endl;
                // }

                myOut << "Is this CC1e0pi in truth: " << kTrueCC1e0pi(&islc) << std::endl;
                myOut << "Is this selected as reco'd CC1e0pi: " << kAutomaticSelection_NoTrigger(&islc) << std::endl;
                myOut << "Is this flash-matched: " << kFlashMatch(&islc) << std::endl;
                myOut << "Is this trigger-flash-matched: " << kTrigFlashMatch(&islc) << std::endl;
                myOut << "Flash-matching info [DZ, T, DZ_Trig]: " << islc.barycenterFM.deltaZ << "\t" << islc.barycenterFM.flashTime << "\t" << islc.barycenterFM.deltaZ_Trigger << std::endl;
            }
        }
        myOut.close();

        return tempSpillVar;
    });   

}