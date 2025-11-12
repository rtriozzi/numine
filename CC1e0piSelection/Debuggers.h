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

        std::ofstream myOut("debug/numinue_NG2Filter.txt", std::ios::app);
        std::ofstream myOutSlim("debug/numinue_NG2Filter_Slim.txt", std::ios::app);     

        for (auto const &islc : sr->slc) {
            if (kTrueCC1e0pi(&islc)) {
                myOut << SourceName << "\t" << sr->hdr.run << "\t" << sr->hdr.evt << "\t" << islc.tmatch.eff << std::endl;
                myOut << "Is this selected as reco'd CC1e0pi: " << kAutomaticSelection_NoTrigger(&islc) << std::endl;
                myOut << "Gap, dE/dx: " << kLargestRecoShower_ConvGap(&islc) << "\t" << kLargestRecoShower_ColldEdx(&islc) << std::endl;
                myOut << "Vtx x reco, truth, DeltaVtx: " << islc.vertex.z << "\t" << islc.truth.position.z << "\t" << kVertex_vsTruth(&islc) << std::endl;
                myOut << "Flash-matching info [DZ, T, DZ_Trig]: " << islc.barycenterFM.deltaZ << "\t" << islc.barycenterFM.flashTime << "\t" << islc.barycenterFM.deltaZ_Trigger << std::endl;
                myOut << "Flash matchinf info charge Z, flash Z: " << islc.barycenterFM.chargeCenter.z << "\t" << islc.barycenterFM.flashCenter.z << std::endl;
                myOut << "Is this flash-matched: " << kFlashMatch(&islc) << std::endl;
                myOut << "Is this trigger-flash-matched: " << kTrigFlashMatch(&islc) << std::endl;

                // slimmed version for further analysis, for selected events
                myOutSlim << SourceName << "\t" << sr->hdr.run << "\t" << sr->hdr.evt << "\t" << kAutomaticSelection_NoTrigger(&islc) << "\t" << islc.tmatch.eff << "\t"
                          << islc.truth.position.z << "\t" << kVertex_vsTruth(&islc) << "\t" << kTrue_NVisProtons(&islc) << "\t"
                          << kFlashMatch(&islc) << "\t" << kTrigFlashMatch(&islc) << "\t" << islc.barycenterFM.chargeCenter.z << "\t" << islc.barycenterFM.flashCenter.z << "\t"
                          << islc.barycenterFM.deltaZ << "\t" << islc.barycenterFM.flashTime << "\t" << islc.barycenterFM.deltaZ_Trigger << std::endl;
            }
        }
        myOut.close();

        return tempSpillVar;
    });   

}