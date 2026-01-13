#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <fstream>
#include <vector>
#include <math.h>
#include "TVector3.h"

#include "CC1e0pi0p_Vars.h"   
#include "CC1e0pi0p_Cuts.h"
#include "CC1e0pi0p_TruthCuts.h"

namespace ana {

    // event dumping
    const SpillMultiVar kEventDump([](const caf::SRSpillProxy* sr) -> std::vector<double>
    {
        
        std::vector<double> tempSpillVar;
        std::string SourceName = sr->hdr.sourceName;

        std::ofstream myOut("debug/CC1e0pi0p_Data.txt", std::ios::app);
        for (auto const &islc : sr->slc) {
            if (kAutomaticSelection(&islc)) {
                // reco information
                myOut << SourceName << "\t" << sr->hdr.run << "\t" << sr->hdr.evt << "\t"
                      << islc.vertex.x << "\t" 
                      << kLargestRecoShower_CollEnergy(&islc) << "\t" << kProton_NuGraph_MaxHIPTag(&islc) << "\t"
                      << kNuGraph_NShrPFPs(&islc) << "\t";
                myOut << std::endl;
            }
        }
        myOut.close();

        return tempSpillVar;
    });

    const SpillMultiVar kMCEventDump([](const caf::SRSpillProxy* sr) -> std::vector<double>
    {
        std::vector<double> tempSpillVar;
        std::string SourceName = sr->hdr.sourceName;

        std::ofstream myOut("debug/CC1e0pi0p_MC_Nom.txt", std::ios::app);   
        for (auto const &islc : sr->slc) {
            // look at selected slices
            if (kAutomaticSelection(&islc)) {
                myOut << SourceName << "\t" << sr->hdr.run << "\t" << sr->hdr.evt << "\t" 
                      << islc.truth.index << "\t" << islc.truth.iscc << "\t" << islc.truth.pdg << "\t" << islc.truth.E << "\t"
                      << islc.truth.position.x << "\t" 
                      << kLargestRecoShower_CollEnergy(&islc) << "\t" << kProton_NuGraph_MaxHIPTag(&islc) << "\t"
                      << kNuGraph_NShrPFPs(&islc) << "\t" << islc.reco.npfp << "\t";
                myOut << std::endl;
            }
        }
        myOut.close();

        return tempSpillVar;
    });   

}