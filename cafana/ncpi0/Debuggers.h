#pragma once

#include <fstream>
#include <vector>
#include <math.h>

// explicit SBNAna includes
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

// my stuff
#include "NCPi0_Vars.h"   
#include "NCPi0_Cuts.h"
#include "NCPi0_TruthCuts.h"

// root stuff
#include "TVector3.h"

namespace ana {

    // event dumping
    const SpillMultiVar kEventDump([](const caf::SRSpillProxy* sr) -> std::vector<double>
    {
        
        std::vector<double> tempSpillVar;
        std::string SourceName = sr->hdr.sourceName;

        std::ofstream myOut("NuMI_Prescaled_NCPi0.txt", std::ios::app);
        for (auto const &islc : sr->slc) {
            if (kAutomaticSelection(&islc)) {
                // reco information
                myOut << sr->hdr.run << "\t" << sr->hdr.evt << "\t" << SourceName << "\t"
                      << islc.vertex.x << "\t" << kPi0_InvariantMass(&islc) << "\t";
                myOut << std::endl;
            }
        }
        myOut.close();

        return tempSpillVar;
    });

}