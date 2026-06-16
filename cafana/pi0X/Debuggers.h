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


    /*
     *  `SpillMultiVar`s for reconstruction/selection efficiency studies,
     *  here with a more flexible output being a text file with the 
     *  truth variable of interest (e.g., true neutrino energy), followed
     *  by a bool signaling whether the interaction passed the sequential
     *  selection cut.
     *  Selection cuts are defined in `CC1e0piSelection_Cuts.h`.
     */

    const SpillMultiVar kPi0X_MakeSelectionEfficiencyOutput([](const caf::SRSpillProxy* sr) -> std::vector<double> {

        // prepare output
        std::ofstream fOut("Pi0X_Efficiency.txt", std::ios::app);

        for (auto const& nu : sr->mc.nu) {
            if (!kIsTrueNuPi0(nu)) continue;
 
            // truth neutrino energy from matched slice
            double truthNuE = -999.;
            for (auto const& islc : sr->slc) {
                if (islc.truth.index == nu.index && kTrueNuPi0(&islc)) {
                    truthNuE = islc.truth.E;
                    break;
                }
            }
            if (truthNuE < -998.) continue;

            // find pi0 among primaries
            int    pi0Idx  = -1;
            double pi0P    = -999.;  ///< pi0 momentum magnitude [GeV]
            double pi0E    = -999.;  ///< pi0 total energy [GeV]
            for (int ip = 0; ip < nu.nprim; ++ip) {
                if (nu.prim[ip].pdg == 111) {
                    TVector3 pi0Mom(nu.prim[ip].startp.x,
                                    nu.prim[ip].startp.y,
                                    nu.prim[ip].startp.z);
                    double p = pi0Mom.Mag();
                    if (p > pi0P) { // take highest-momentum pi0 if more than one
                        pi0P   = p;
                        pi0E   = nu.prim[ip].startE;
                        pi0Idx = ip;
                    }
                }
            }

            // cos(theta_gamma_gamma) minimum opening angle from pi0 kinematics
            // for a symmetric decay: cos(theta_min) = 1 - 2/(gamma+1), gamma = E/m
            const double mPi0 = 0.13498; // GeV
            double cosGammaGamma = -999.;
            if (pi0Idx >= 0 && pi0E > mPi0) {
                double gamma = pi0E / mPi0;
                cosGammaGamma = 1. - 2. / (gamma + 1.);
            }

            // cos(pi0, NuMI beam direction)
            double cosPi0NuMI = -999.;
            double cosPi0Z    = -999.;
            if (pi0Idx >= 0) {
                TVector3 pi0Dir(nu.prim[pi0Idx].startp.x,
                                nu.prim[pi0Idx].startp.y,
                                nu.prim[pi0Idx].startp.z);
                if (pi0Dir.Mag() > 0) {
                    pi0Dir   = pi0Dir.Unit();
                    cosPi0Z    = pi0Dir.z();
                    cosPi0NuMI = pi0Dir.x() * 3.94583e-01 +
                                pi0Dir.y() * 4.26067e-02 +
                                pi0Dir.z() * 9.17677e-01;
                }
            }

            // dump to file
            fOut << truthNuE    << "\t"
                << pi0P        << "\t"
                << cosGammaGamma << "\t"
                << cosPi0Z     << "\t"
                << cosPi0NuMI;
 
            // apply cuts
            for (auto const& step : SelectionSteps) {
                bool passed = false;
                for (auto const& islc : sr->slc) {
                    if (islc.truth.index == nu.index &&
                        kTrueNuPi0(&islc)          &&
                        step.cut(&islc)) {
                        passed = true;
                        break;
                    }
                }
                fOut << "\t" << (passed ? 1 : 0);
            }
            fOut << "\n";
        }

        fOut.close();

        return {};
    });

}