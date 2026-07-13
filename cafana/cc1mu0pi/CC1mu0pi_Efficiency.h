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

    // true energy of signal true neutrinos
    const SpillMultiVar kCC1mu0p1Signal_TrueNeutrinoEnergy([](const caf::SRSpillProxy* sr)-> std::vector<double> {

        std::vector<double> trueEnergies;

        for (auto const& nu : sr->mc.nu) { 
            if (kIsTrueCC1mu0pi(nu)) {
                trueEnergies.push_back(nu.E);    
            }
        }
        
        return trueEnergies;
    });

    const SpillMultiVar kCC1mu0p1Signal_TrueMuonEnergy([](const caf::SRSpillProxy* sr)-> std::vector<double> {

        std::vector<double> trueEnergies;

        for (auto const& nu : sr->mc.nu) { 
            if (kIsTrueCC1mu0pi(nu)) {

                // get the true electron
                for (int ip(0); ip < nu.nprim; ++ip) {
                    if (abs(nu.prim[ip].pdg) == 13) {
                        trueEnergies.push_back(nu.prim[ip].startE - nu.prim[ip].endE);   
                        break;
                    }
                } 
            }
        }
        
        return trueEnergies;
    });

    const SpillMultiVar kCC1mu0p1Signal_NoPileup_TrueNeutrinoEnergy([](const caf::SRSpillProxy* sr)-> std::vector<double> {

        std::vector<double> trueEnergies;

        for (auto const& nu : sr->mc.nu) { 
            if (kIsTrueCC1mu0pi(nu) && (sr->mc.nnu == 1)) {
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

    // factory of true energies for reconstructed and selected neutrinos matched to the truth, at the selection step <cut>
    SpillMultiVar kCC1mu0p1Signal_TrueNeutrinoEnergy_MakeSelectionStep(const Cut& cut)
    {
        return SpillMultiVar([cut](const caf::SRSpillProxy* sr) -> std::vector<double> {
            bool trueNeutrinoWasCounted;
            std::vector<double> selectedEnergies;

            for (auto const& nu : sr->mc.nu) {
                if (kIsTrueCC1mu0pi(nu)) {
                    trueNeutrinoWasCounted = false;

                    for (auto const& islc : sr->slc) {
                        if ((islc.truth.index == nu.index) &&  ///< same index, to account for pile-up
                            (kTrueCC1mu0pi(&islc)) &&          ///< signal slice
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

    SpillMultiVar kCC1mu0p1Signal_TrueMuonEnergy_MakeSelectionStep(const Cut& cut)
    {
        return SpillMultiVar([cut](const caf::SRSpillProxy* sr) -> std::vector<double> {
            bool trueNeutrinoWasCounted;
            std::vector<double> selectedEnergies;

            for (auto const& nu : sr->mc.nu) {
                if (kIsTrueCC1mu0pi(nu)) {
                    trueNeutrinoWasCounted = false;

                    for (auto const& islc : sr->slc) {
                        if ((islc.truth.index == nu.index) &&  ///< same index, to account for pile-up
                            (kTrueCC1mu0pi(&islc)) &&           ///< signal slice
                            (cut(&islc)) &&                    ///< apply reconstruction cut step
                            (!trueNeutrinoWasCounted)) {       ///< make sure to count only one slice per true neutrino

                            // get the true muon
                            for (int ip(0); ip < nu.nprim; ++ip) {
                                if (abs(nu.prim[ip].pdg) == 13) {
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

    SpillMultiVar kCC1mu0p1Signal_NoPileup_TrueNeutrinoEnergy_MakeSelectionStep(const Cut& cut)
    {
        return SpillMultiVar([cut](const caf::SRSpillProxy* sr) -> std::vector<double> {
            bool trueNeutrinoWasCounted;
            std::vector<double> selectedEnergies;

            for (auto const& nu : sr->mc.nu) {
                if (kIsTrueCC1mu0pi(nu) && (sr->mc.nnu == 1)) {
                    trueNeutrinoWasCounted = false;

                    for (auto const& islc : sr->slc) {
                        if ((islc.truth.index == nu.index) &&  ///< same index, to account for pile-up
                            (kTrueCC1mu0pi(&islc)) &&          ///< signal slice
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

    const SpillMultiVar kCC1mu0piSignal_MakeSelectionEfficiencyOutput([](const caf::SRSpillProxy* sr) -> std::vector<double> {

        // prepare output
        std::ofstream fOut("CC1mu0piSelection_Efficiency.txt", std::ios::app);
        // fOut << "true_nu_energy_GeV" << "\t" << "true_elec_energy_GeV" << "\t" << "true_proton_KE_GeV" << "\t" << "true_ep_cosAngle_3D" << "\t" << "true_elec_cosAngle_NuMI";
        // for (auto const& step : SelectionSteps)
        //     fOut << "\t" << step.suffix;
        // fOut << "\n";

        for (auto const& nu : sr->mc.nu) {
            if (!kIsTrueCC1mu0pi(nu)) continue;

            // truth neutrino energy from matched slice
            double truthNuE = -999.;
            for (auto const& islc : sr->slc) {
                if (islc.truth.index == nu.index && kTrueCC1mu0pi(&islc)) {
                    truthNuE = islc.truth.E;
                    break;
                }
            }
            if (truthNuE < -998.) continue;

            // find leading electron and leading proton among primaries
            int    muIdx   = -1;
            int    protonIdx = -1;
            double muLen     = -999.;
            double protonKE  = -999.;
            for (int ip = 0; ip < nu.nprim; ++ip) {
                if (std::abs(nu.prim[ip].pdg) == 13 && muIdx == -1) {
                    muIdx = ip;
                    muLen = nu.prim[muIdx].length;
                }
                if (nu.prim[ip].pdg == 2212) {
                    double ke = nu.prim[ip].startE - nu.prim[ip].endE;
                    if (ke > protonKE) { protonKE = ke; protonIdx = ip; }
                }
            }

            // cos(electron, NuMI) 
            double cosMuNuMI = -999.;
            double cosMuZ = -999.;
            if (muIdx >= 0) {
                TVector3 muDir(nu.prim[muIdx].startp.x,
                                nu.prim[muIdx].startp.y,
                                nu.prim[muIdx].startp.z);
                if (muDir.Mag() > 0)
                    muDir = muDir.Unit();
                    cosMuZ = muDir.z();
                    cosMuNuMI = muDir.x() * 3.94583e-01 +
                        muDir.y() * 4.26067e-02 +
                        muDir.z() * 9.17677e-01;
            }

            // cos(muon, proton) 3D opening angle
            double cosMP = -999.;
            if (muIdx >= 0 && protonIdx >= 0) {
                TVector3 muDir(nu.prim[muIdx].startp.x,
                                 nu.prim[muIdx].startp.y,
                                 nu.prim[muIdx].startp.z);
                TVector3 protonDir(nu.prim[protonIdx].startp.x,
                                   nu.prim[protonIdx].startp.y,
                                   nu.prim[protonIdx].startp.z);
                if (muDir.Mag() > 0 && protonDir.Mag() > 0)
                    cosMP = std::cos(muDir.Angle(protonDir));
            }

            // dump to file
            fOut << truthNuE << "\t" << muLen << "\t" << protonKE << "\t" << cosMuZ << "\t" << cosMuNuMI << "\t" << cosMP;

            // apply cuts
            for (auto const& step : SelectionSteps) {
                bool passed = false;
                for (auto const& islc : sr->slc) {
                    if (islc.truth.index == nu.index &&
                        kTrueCC1mu0pi(&islc)          &&
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