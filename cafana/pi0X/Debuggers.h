#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <fstream>
#include <vector>
#include <math.h>
#include "TVector3.h"

#include "Pi0X_Vars.h"   
#include "Pi0X_Cuts.h"
#include "Pi0X_TruthCuts.h"

namespace ana {

    // event dumping
    const SpillMultiVar kEventDump([](const caf::SRSpillProxy* sr) -> std::vector<double>
    {
        
        std::vector<double> tempSpillVar;
        std::string SourceName = sr->hdr.sourceName;

        std::ofstream myOut("NuMI_Prescaled_Pi0XSelection.txt", std::ios::app);
        for (auto const &islc : sr->slc) {
            if (kPi0Selection(&islc)) {
                // reco information
                myOut << sr->hdr.run << "\t" << sr->hdr.evt << "\t" << SourceName << "\t"
                      << islc.vertex.x << "\t" 
                      << kLargestRecoShower_CollEnergy(&islc) << "\t" 
                      << kSubleadRecoShower_CollEnergy(&islc) << "\t"
                      << kPi0_InvariantMass(&islc) << "\t";
                myOut << std::endl;
            }
        }
        myOut.close();

        return tempSpillVar;
    });
}