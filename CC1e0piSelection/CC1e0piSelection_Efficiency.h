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

    // true energy of signal true neutrinos
    const SpillMultiVar kCC1e0p1Signal_TrueNeutrinoEnergy([](const caf::SRSpillProxy* sr)-> std::vector<double> {

        std::vector<double> trueEnergies;

        for (auto const& nu : sr->mc.nu) { 
            if (kIsTrueCC1e0pi(nu)) {
                trueEnergies.push_back(nu.E);    
            }
        }
        
        return trueEnergies;
    });

    // true energy of signal true neutrinos
    const SpillMultiVar kNuCC_TrueNeutrinoEnergy([](const caf::SRSpillProxy* sr)-> std::vector<double> {

        std::vector<double> trueEnergies;

        for (auto const& nu : sr->mc.nu) { 
            if ((nu.iscc) &&
                (kIsInFV(nu.position.x, nu.position.y, nu.position.z)) &&
                (sr->mc.nnu > 0)) {
                trueEnergies.push_back(nu.E);    
            }
        }
        
        return trueEnergies;
    });

    const SpillMultiVar kCC1e0p1Signal_NoPileup_TrueNeutrinoEnergy([](const caf::SRSpillProxy* sr)-> std::vector<double> {

        std::vector<double> trueEnergies;

        for (auto const& nu : sr->mc.nu) { 
            if (kIsTrueCC1e0pi(nu) && (sr->mc.nnu == 1)) {
                trueEnergies.push_back(nu.E);    
            }
        }
        
        return trueEnergies;
    });

    // factory of true energies for reconstructed and selected neutrinos matched to the truth
    SpillMultiVar kCC1e0p1Signal_TrueNeutrinoEnergy_MakeSelectionStep(const Cut& cut)
    {
        return SpillMultiVar([cut](const caf::SRSpillProxy* sr) -> std::vector<double> {
            bool trueNeutrinoWasCounted;
            std::vector<double> selectedEnergies;

            for (auto const& nu : sr->mc.nu) {
                if (kIsTrueCC1e0pi(nu)) {
                    trueNeutrinoWasCounted = false;

                    for (auto const& islc : sr->slc) {
                        if ((islc.truth.index == nu.index) &&  ///< same index, to account for pile-up
                            (kTrueCC1e0pi(&islc)) &&           ///< signal slice
                            (cut(&islc)) &&                    ///< apply reconstruction cut step
                            (!trueNeutrinoWasCounted)) {       ///< make sure to count only one slice per true neutrino

                            selectedEnergies.push_back(islc.truth.E);
                            trueNeutrinoWasCounted = true;
                        }
                    }
                }
            }

            return selectedEnergies;
        });
    }

    // factory of true energies for reconstructed and selected neutrinos matched to the truth
    SpillMultiVar kNuCC_TrueNeutrinoEnergy_MakeSelectionStep(const Cut& cut)
    {
        return SpillMultiVar([cut](const caf::SRSpillProxy* sr) -> std::vector<double> {
            bool trueNeutrinoWasCounted;
            std::vector<double> selectedEnergies;

            for (auto const& nu : sr->mc.nu) {
                if ((nu.iscc) &&
                    (kIsInFV(nu.position.x, nu.position.y, nu.position.z)) &&
                    (sr->mc.nnu > 0)) {
                    trueNeutrinoWasCounted = false;

                    for (auto const& islc : sr->slc) {
                        if ((islc.truth.index == nu.index) &&  ///< same index, to account for pile-up
                            (kIsCC(&islc)) &&                  ///< signal slice
                            (kIsNuinFV(&islc)) &&
                            (cut(&islc)) &&                    ///< apply reconstruction cut step
                            (!trueNeutrinoWasCounted)) {       ///< make sure to count only one slice per true neutrino

                            selectedEnergies.push_back(islc.truth.E);
                            trueNeutrinoWasCounted = true;
                        }
                    }
                }
            }

            return selectedEnergies;
        });
    }

    SpillMultiVar kCC1e0p1Signal_NoPileup_TrueNeutrinoEnergy_MakeSelectionStep(const Cut& cut)
    {
        return SpillMultiVar([cut](const caf::SRSpillProxy* sr) -> std::vector<double> {
            bool trueNeutrinoWasCounted;
            std::vector<double> selectedEnergies;

            for (auto const& nu : sr->mc.nu) {
                if (kIsTrueCC1e0pi(nu) && (sr->mc.nnu == 1)) {
                    trueNeutrinoWasCounted = false;

                    for (auto const& islc : sr->slc) {
                        if ((islc.truth.index == nu.index) &&  ///< same index, to account for pile-up
                            (kTrueCC1e0pi(&islc)) &&           ///< signal slice
                            (cut(&islc)) &&                    ///< apply reconstruction cut step
                            (!trueNeutrinoWasCounted)) {       ///< make sure to count only one slice per true neutrino

                            selectedEnergies.push_back(islc.truth.E);
                            trueNeutrinoWasCounted = true;
                        }
                    }
                }
            }

            return selectedEnergies;
        });
    }

    // debug selection
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

}