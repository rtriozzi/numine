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

    /*
     *  Truth-level variables for signal neutrinos.
     */

    const SpillMultiVar kCC1e0p1Signal_TrueNeutrinoEnergy([](const caf::SRSpillProxy* sr)-> std::vector<double> {

        std::vector<double> trueEnergies;

        for (auto const& nu : sr->mc.nu) { 
            if (kIsTrueCC1e0pi(nu)) {
                trueEnergies.push_back(nu.E);    
            }
        }
        
        return trueEnergies;
    });

    const SpillMultiVar kCC1e0p1Signal_TrueElectronEnergy([](const caf::SRSpillProxy* sr)-> std::vector<double> {

        std::vector<double> trueEnergies;

        for (auto const& nu : sr->mc.nu) { 
            if (kIsTrueCC1e0pi(nu)) {

                // get the true electron
                for (int ip(0); ip < nu.nprim; ++ip) {
                    if (abs(nu.prim[ip].pdg) == 11) {
                        trueEnergies.push_back(nu.prim[ip].startE - nu.prim[ip].endE);   
                        break;
                    }
                } 
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

    /*
     *  Factory of `SpillMultiVar`s for reconstruction/selection efficiency studies,
     *  having as output the list of selected events in the truth variable of interest,
     *  for plotting the corresponding "selected" histogram.
     *  A factory is needed to handle, e.g., pile-up of multiple true neutrino
     *  interactions in the beam spill, which happens somewhat frequently in NuMI.
     *  The factory wants the corresponding `Cut` as input.
     */

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

    SpillMultiVar kCC1e0p1Signal_TrueElectronEnergy_MakeSelectionStep(const Cut& cut)
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

                            // get the true electron
                            for (int ip(0); ip < nu.nprim; ++ip) {
                                if (abs(nu.prim[ip].pdg) == 11) {
                                    selectedEnergies.push_back(nu.prim[ip].startE - nu.prim[ip].endE);
                                    trueNeutrinoWasCounted = true;
                                    break;
                                }
                            } 
                        }
                    }
                }
            }

            return selectedEnergies;
        });
    }

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

    SpillMultiVar kCC1e0piSignal_NoPileup_TrueNeutrinoEnergy_MakeSelectionStep(const Cut& cut)
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

    /*
    *  `SpillMultiVar`s for reconstruction/selection efficiency studies,
    *  here with a more flexible output being a text file with the 
    *  truth variable of interest (e.g., true neutrino energy), followed
    *  by a bool signaling whether the interaction passed the sequential
    *  selection cut.
    *  Selection cuts are defined in `CC1e0piSelection_Cuts.h`.
    */

    const SpillMultiVar kCC1e0piSignal_MakeSelectionEfficiencyOutput([](const caf::SRSpillProxy* sr) -> std::vector<double> {

        // prepare output
        std::ofstream fOut("CC1e0piSelection_Efficiency.txt", std::ios::app);
        // fOut << "true_nu_energy_GeV" << "\t" << "true_elec_energy_GeV" << "\t" << "true_proton_KE_GeV" << "\t" << "true_ep_cosAngle_3D" << "\t" << "true_elec_cosAngle_NuMI";
        // for (auto const& step : SelectionSteps)
        //     fOut << "\t" << step.suffix;
        // fOut << "\n";

        for (auto const& nu : sr->mc.nu) {
            if (!kIsTrueCC1e0pi(nu)) continue;

            // truth neutrino energy from matched slice
            double truthNuE = -999.;
            for (auto const& islc : sr->slc) {
                if (islc.truth.index == nu.index && kTrueCC1e0pi(&islc)) {
                    truthNuE = islc.truth.E;
                    break;
                }
            }
            if (truthNuE < -998.) continue;

            // find leading electron and leading proton among primaries
            int    elecIdx   = -1;
            int    protonIdx = -1;
            double elecE     = -999.;
            double protonKE  = -999.;
            for (int ip = 0; ip < nu.nprim; ++ip) {
                if (std::abs(nu.prim[ip].pdg) == 11 && elecIdx == -1) {
                    elecIdx = ip;
                    elecE = nu.prim[elecIdx].startE;
                }
                if (nu.prim[ip].pdg == 2212) {
                    double ke = nu.prim[ip].startE - nu.prim[ip].endE;
                    if (ke > protonKE) { protonKE = ke; protonIdx = ip; }
                }
            }

            // cos(electron, NuMI) 
            double cosElecNuMI = -999.;
            double cosElecZ = -999.;
            if (elecIdx >= 0) {
                TVector3 elecDir(nu.prim[elecIdx].startp.x,
                                 nu.prim[elecIdx].startp.y,
                                 nu.prim[elecIdx].startp.z);
                if (elecDir.Mag() > 0)
                    elecDir = elecDir.Unit();
                    cosElecZ = elecDir.z();
                    cosElecNuMI = elecDir.x() * 3.94583e-01 +
                        elecDir.y() * 4.26067e-02 +
                        elecDir.z() * 9.17677e-01;
            }

            // cos(electron, proton) 3D opening angle
            double cosEP = -999.;
            if (elecIdx >= 0 && protonIdx >= 0) {
                TVector3 elecDir(nu.prim[elecIdx].startp.x,
                                 nu.prim[elecIdx].startp.y,
                                 nu.prim[elecIdx].startp.z);
                TVector3 protonDir(nu.prim[protonIdx].startp.x,
                                   nu.prim[protonIdx].startp.y,
                                   nu.prim[protonIdx].startp.z);
                if (elecDir.Mag() > 0 && protonDir.Mag() > 0)
                    cosEP = std::cos(elecDir.Angle(protonDir));
            }

            // dump to file
            fOut << truthNuE << "\t" << elecE << "\t" << protonKE << "\t" << cosElecZ << "\t" << cosElecNuMI << "\t" << cosEP;

            // apply cuts
            for (auto const& step : SelectionSteps) {
                bool passed = false;
                for (auto const& islc : sr->slc) {
                    if (islc.truth.index == nu.index &&
                        kTrueCC1e0pi(&islc)          &&
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