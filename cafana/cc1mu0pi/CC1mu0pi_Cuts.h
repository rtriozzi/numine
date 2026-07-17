#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <fstream>
#include <vector>
#include <math.h>
#include "TVector3.h"
#include "TGraph.h"
#include "TF1.h"

#include "CC1mu0pi_Vars.h"

namespace ana {

    // pre-selection
    const SpillCut kCRTPMTNeutrino([](const caf::SRSpillProxy* sr) {

        double minTime = 0., maxTime = 0.;

        for ( const auto& match: sr->crtpmt_matches ) {
            if (sr->hdr.ismc) { minTime = 0.0; maxTime = 9.8; }
            if (!sr->hdr.ismc) { minTime = -0.4; maxTime = 10.5; }

            // in-time flash and CRT-PMT match
            if(match.flashGateTime > minTime && 
               match.flashGateTime < maxTime && 
               match.flashClassification == 0)  
                return true;
        }

        return false;
    });

    const Cut kNotClearCosmic([](const caf::SRSliceProxy* slc) { 
        return !(slc->is_clear_cosmic);
    });

    const Cut kVertexInFV([](const caf::SRSliceProxy* slc) { 
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;

        return kIsInFV(slc->vertex.x, slc->vertex.y, slc->vertex.z);
    });

    // light matching
    const Cut kTrigFlashMatch([](const caf::SRSliceProxy* slc) { 
        return (slc->barycenterFM.deltaZ_Trigger >= 0 && 
                slc->barycenterFM.deltaZ_Trigger <= 100);        
    });

    const Cut kFlashMatch([](const caf::SRSliceProxy* slc) { 
        return (slc->barycenterFM.deltaZ >= 0 && 
                slc->barycenterFM.deltaZ <= 100 &&
                slc->barycenterFM.flashTime > -1 &&
                slc->barycenterFM.flashTime < 11);        
    });

    // muon identification
    const Cut kMuon_LengthCut([](const caf::SRSliceProxy* slc) { 
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx < 0) return false;
        if (std::isnan(slc->reco.pfp[muonIdx].trk.len)) return false;

        return slc->reco.pfp[muonIdx].trk.len >= 50;
    });

    const Cut kMuon_Containment([](const caf::SRSliceProxy* slc) { 
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx < 0) return false;
        if (std::isnan(slc->reco.pfp[muonIdx].trk.end.x) || std::isnan(slc->reco.pfp[muonIdx].trk.end.y) || std::isnan(slc->reco.pfp[muonIdx].trk.end.z)) 
            return false;

        return kIsInContained(slc->reco.pfp[muonIdx].trk.end.x, slc->reco.pfp[muonIdx].trk.end.y, slc->reco.pfp[muonIdx].trk.end.z);
    });

    // Np0π selection
    const Cut kNSelectedProtons([](const caf::SRSliceProxy* slc) { 
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx < 0) return false;

        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);

        return selectedProtonIdx.size() > 0;
    }); 

    const Cut kNoOtherParticle([](const caf::SRSliceProxy* slc) { 
        const int muonIdx = kMuonIdx(slc);
        if (muonIdx <= 0) return false;

        int NOtherParticles = kNSelectedProtonsIdx_NOtherParticles(slc);
        if (NOtherParticles < 0) return false;

        return NOtherParticles == 0;
    }); 

    // good NuMI Run1/Run2 runs
    bool IsNuMIRun1Run2GoodRun(int run) {
    static const std::vector<int> kGoodRuns = {
        8461,  8462,  8468,  8469,  8470,  8471,
        8505,  8506,  8507,  8513,  8514,  8515,
        8521,  8522,  8527,  8528,  8529,  8530,
        8531,  8552,  8553,  9593,  9594,  9595,
        9597,  9599,  9602,  9610,  9642,  9646,
        9648,  9649,  9688,  9690,  9691,  9692,
        9693,  9694,  9695,  9696,  9699,  9700,
        9704,  9715,  9716,  9717,  9721,  9723,
        9725,  9726,  9728,  9729,  9730,  9731,
        9732,  9733,  9735,  9743,  9744,  9745,
        9746,  9747,  9750,  9752,  9753,  9755,
        9762,  9763,  9764,  9765,  9781,  9791,
        9792,  9794,  9796,  9807,  9834,  9835,
        9837,  9838,  9840,  9844,  9847,  9849,
        9851,  9854,  9855,  9860,  9862,  9867,
        9869,  9870,  9892,  9894,  9896,  9897,
        9914,  9919,  9921,  9922,  9924,  9925,
        9926,  9944,  9945,  9949,  9950,  9951,
        9953,  9954,  9956,  9959,  9960,  9972,
        9974,  9977,  9979,  9981,  10054, 10059,
        10060, 10061, 10064, 10065, 10066, 10084,
        10085, 10096, 10097
    };
        return std::find(kGoodRuns.begin(), kGoodRuns.end(), run) != kGoodRuns.end();
    }

    const SpillCut kGoodRunCut([](const caf::SRSpillProxy* sr) {
        return IsNuMIRun1Run2GoodRun(sr->hdr.run);
    });

    // beam quality cuts
    const SpillCut kNuMISpillHornCurrentCut ( [](const caf::SRSpillProxy *sr)
    {
        double iHorn = ((sr->hdr.spillnumiinfo.NSLINA - 0.01)/0.9951) 
            + ((sr->hdr.spillnumiinfo.NSLINB - (-0.14))/0.9957) 
            + ((sr->hdr.spillnumiinfo.NSLINC - (-0.05))/0.9965) 
            + ((sr->hdr.spillnumiinfo.NSLIND - (-0.07))/0.9945);
        return (iHorn > -202. && iHorn <= -196.4);
    });

    const SpillCut kNuMISpillPOTCut ( [](const caf::SRSpillProxy *sr)
    {
        
        double pot = sr->hdr.spillnumiinfo.TRTGTD;
        if ( sr->hdr.spillnumiinfo.TRTGTD < 0.02 ) {
            pot = sr->hdr.spillnumiinfo.TR101D; // fall back to TR101D
            if ( pot < 0. ) pot = 0.; // clip
        }

        return pot > 2.0e12;
    });

    // Beam Pos at target
    //////////////////////////////////////////////////////////////
    double extrapToLoc ( const double var1, const double loc1, const double var2, const double loc2, const double loc3 )
    {
        double slope = (var2-var1)/(loc2-loc1);
        double proj = var1 + ((loc3-loc1)*slope);
        return proj;
    }

    // NOTE: this function does averages THEN projects, I believe NOvA projects then averages. Maybe that's better...
    std::pair<double,double> BeamPositionAtTargetVal ( const caf::SRNuMIInfoProxy& info, const unsigned int runNumber )
    {
        double mmPerFoot = 12.*2.54*10.; // inch per foot * cm per inch * mm per cm

        // From NOvA, converted to mm
        double z_hp121 = -68.04458 * mmPerFoot;
        double z_vp121 = -66.99283 * mmPerFoot;
        double z_hptgt = -31.25508 * mmPerFoot;
        double z_vptgt = -30.16533 * mmPerFoot;
        double z_targ  = 0.; 

        // If any monitor doesn't have 7 then return -999, -999
        if ( info.HP121.size() != 7 || info.VP121.size() != 7 || info.HPTGT.size() != 7 ||
            info.HITGT.size() != 7 || info.VPTGT.size() != 7 || info.VITGT.size() != 7 ) {
        std::pair<double,double> ret = std::make_pair<double,double>(-999.,-999.);
        return ret;
        }

        double ave_HP121(0.), ave_VP121(0.), ave_HPTGT(0.), ave_VPTGT(0.), wgt_HPTGT(0.), wgt_VPTGT(0.);
        for ( unsigned int idxSpillPart=0; idxSpillPart<7; ++idxSpillPart ) {
        ave_HP121+=info.HP121[idxSpillPart];
        ave_VP121+=info.VP121[idxSpillPart];
        ave_HPTGT+=info.HPTGT[idxSpillPart]*info.HITGT[idxSpillPart];
        wgt_HPTGT+=info.HITGT[idxSpillPart];
        ave_VPTGT+=info.VPTGT[idxSpillPart]*info.VITGT[idxSpillPart];
        wgt_VPTGT+=info.VITGT[idxSpillPart];
        }

        ave_HP121/=7.;
        ave_VP121/=7.;
        if ( wgt_HPTGT < std::numeric_limits<double>::epsilon() ||
            wgt_VPTGT < std::numeric_limits<double>::epsilon() ) {
        std::pair<double,double> ret = std::make_pair<double,double>(-999.,-999.);
        return ret;
        }
        ave_HPTGT/=wgt_HPTGT;
        ave_VPTGT/=wgt_VPTGT;

        // TODO: double check that these offsets make sense...
        // Based on values provided by NuMI Target group. Thanks!
        double xCorr(-9999.);
        double yCorr(-9999.);
        // Run 1 offsets -- for now apply same cuts to everything before Run1 as well
        if ( runNumber <= 8553 ) {
        double x_nom121(0.), x_nomtgt(0.), y_nom121(-1.5), y_nomtgt(-1.5), x_nomTarg(0.398), y_nomTarg(-0.39);
        
        xCorr = extrapToLoc(x_nom121,z_hp121,x_nomtgt,z_hptgt,z_targ) - x_nomTarg;
        yCorr = extrapToLoc(y_nom121,z_hp121,y_nomtgt,z_hptgt,z_targ) - y_nomTarg;
        }
        // Run 2 offsets -- will have to edit this to give the max run at the end of Run 2, or when conditions change.
        else if ( runNumber > 8553 ) {
        double x_nom121(1.2), x_nomtgt(0.7), y_nom121(-0.4), y_nomtgt(-0.94), x_nomTarg(0.03), y_nomTarg(-0.59);
        
        xCorr = extrapToLoc(x_nom121,z_hp121,x_nomtgt,z_hptgt,z_targ) - x_nomTarg;
        yCorr = extrapToLoc(y_nom121,z_hp121,y_nomtgt,z_hptgt,z_targ) - y_nomTarg;
        }

        std::pair<double,double> ret = std::make_pair<double,double>( extrapToLoc(ave_HP121,z_hp121,ave_HPTGT,z_hptgt,z_targ) - xCorr,
                                                                    extrapToLoc(ave_VP121,z_vp121,ave_VPTGT,z_vptgt,z_targ) - yCorr );

        return ret;
    }

    const SpillCut kNuMISpillBeamPosCut ( [](const caf::SRSpillProxy *sr)
    {
        std::pair<double,double> beamPos = BeamPositionAtTargetVal(
            sr->hdr.spillnumiinfo, sr->hdr.run
        );

        return fabs(beamPos.first) <= 1. && fabs(beamPos.second) <= 1.;  
    });

    // Beam Widths
    // This also motivated by/borrows what is done in the NOvA modules
    // is it [103,103+48) or (103,103+48] ? -- I think the first?
    //////////////////////////////////////////////////////////////
    std::pair<double,double> BeamWidthVal ( const caf::SRNuMIInfoProxy& info )
    {
        if ( info.MTGTDS.size() < 199 ) {
        std::pair<double,double> ret = std::make_pair<double,double>(-999.,-999.);
        return ret;
        }

        double x_pos[48];
        double h_chan[48];
        double v_chan[48];

        for ( unsigned int idxX=0; idxX<48; ++idxX ) {
        x_pos[ idxX ] = 0.5 * (idxX+1-24.5);
        }

        double h_sum = 0.;
        double h_min = 99999;
        for ( unsigned int idxH=0; idxH<48; ++idxH ) {
        h_chan[ idxH ] = -1.*info.MTGTDS.at(103+idxH);
        h_sum += h_chan[ idxH ];
        if(h_chan[ idxH ]<h_min){
            h_min = h_chan[ idxH ];
        }
        }

        double v_sum = 0.;
        double v_min = 99999;
        for ( unsigned int idxV=0; idxV<48; ++idxV ) {
        v_chan[ idxV ] = -1.*info.MTGTDS.at(151+idxV);
        v_sum += v_chan[ idxV ];
        if(v_chan[ idxV ]<v_min){
            v_min = v_chan[ idxV ];
        }
        }

        // H width
        TGraph *hGraph = new TGraph(48, x_pos, h_chan);
        TF1 *hGauss = new TF1("hGauss","[0] + [1]*TMath::Exp(-((x-[2])*(x-[2]))/(2*[3]*[3]))");
        hGauss->SetParameter(0,-1.*h_sum/48.);
        hGauss->SetParameter(1,-1.*h_min);
        hGauss->SetParameter(2,0.);
        hGauss->SetParameter(3,2.);
        hGauss->SetParLimits(3, 0., 10.);
        /*
            // Remove 0s
            std::vector<int> to_remove_h;
            for ( int idx=0; idx<hGraph->GetN(); ++idx ) {
            double x, y;
            hGraph->GetPoint( idx, x, y );
            if ( fabs(y) < std::numeric_limits<double>::epsilon() ) to_remove_h.push_back( idx );
            }
            for ( unsigned int idx=0; idx<to_remove_h.size(); ++idx ) hGraph->RemovePoint( to_remove_h.at( to_remove_h.size()-1-idx ) );
        */
        hGraph->Fit(hGauss, "Q");

        // V width
        TGraph *vGraph = new TGraph(48, x_pos, v_chan);
        TF1 *vGauss = new TF1("vGauss","[0] + [1]*TMath::Exp(-((x-[2])*(x-[2]))/(2*[3]*[3]))");
        vGauss->SetParameter(0,-1.*v_sum/48.);
        vGauss->SetParameter(1,-1.*v_min);
        vGauss->SetParameter(2,0.);
        vGauss->SetParameter(3,2.);
        vGauss->SetParLimits(3, 0., 10.);

        vGraph->Fit(vGauss,"Q");

        std::pair<double,double> ret = std::make_pair<double,double>(hGauss->GetParameter(3),vGauss->GetParameter(3));

        hGraph->~TGraph();
        hGauss->~TF1();
        vGraph->~TGraph();
        vGauss->~TF1();

        return ret;
    }

    const SpillCut kNuMISpillBeamWidthCut ( [](const caf::SRSpillProxy *sr)
    {
        std::pair<double,double> beamWidths = BeamWidthVal(sr->hdr.spillnumiinfo);

        return beamWidths.first > 0.57 && beamWidths.first <= 1.88 &&
            beamWidths.second > 0.57 && beamWidths.second <= 1.88;
    });

    // beam quality cuts
    const SpillCut kNuMISpillQualityCut = kNuMISpillHornCurrentCut && kNuMISpillPOTCut && kNuMISpillBeamPosCut && kNuMISpillBeamWidthCut;

    // automatic selection
    const Cut kAutomaticNuMuSelection = kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kMuon_LengthCut && kMuon_Containment &&
                                    kNSelectedProtons && kNoOtherParticle;

    const Cut kAutomaticNuMuSelection_NoTrigger = kNotClearCosmic && kVertexInFV && kFlashMatch && kMuon_LengthCut && kMuon_Containment &&
                                              kNSelectedProtons && kNoOtherParticle;

    const Cut kPreNuMuSelection = kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kMuon_LengthCut && kMuon_Containment;

    const Cut kPreNuMuSelection_NoTrigger = kNotClearCosmic && kVertexInFV && kFlashMatch && kMuon_LengthCut && kMuon_Containment;

    // selections
    struct SelDef {
        std::string suffix = "";
        std::string label = "";
        Cut cut = kNoCut;
        int color = kBlack;
    };

    std::vector<SelDef> SelectionSteps = {
        {"presel", "Presel.",               kNotClearCosmic && kVertexInFV,     kBlack},
        {"flash",  "FM",                    kNotClearCosmic && kVertexInFV && kTrigFlashMatch,     kRed+2},
        {"mulen", "Muon ID",                kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kMuon_LengthCut,     kRed-7},
        {"mucont", "Cont'd",                kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kMuon_LengthCut && kMuon_Containment,     kOrange-3},
        {"proton", "Proton ID",             kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kMuon_LengthCut && kMuon_Containment && kNSelectedProtons,   kGreen-2},
        {"nothingelse", "0#pi",             kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kMuon_LengthCut && kMuon_Containment && kNSelectedProtons && kNoOtherParticle,   kGreen-7},
    };

    std::vector<SelDef> SelectionSteps_NoTrigger = {
        {"presel", "Presel.",               kNotClearCosmic && kVertexInFV,     kBlack},
        {"flash",  "FM",                    kNotClearCosmic && kVertexInFV && kFlashMatch,     kRed+2},
        {"mulen", "Muon ID",                kNotClearCosmic && kVertexInFV && kFlashMatch && kMuon_LengthCut,     kRed-7},
        {"mucont", "Cont'd",                kNotClearCosmic && kVertexInFV && kFlashMatch && kMuon_LengthCut && kMuon_Containment,     kOrange-3},
        {"proton", "Proton ID",             kNotClearCosmic && kVertexInFV && kFlashMatch && kMuon_LengthCut && kMuon_Containment && kNSelectedProtons,   kGreen-2},
        {"nothingelse", "0#pi",             kNotClearCosmic && kVertexInFV && kFlashMatch && kMuon_LengthCut && kMuon_Containment && kNSelectedProtons && kNoOtherParticle,   kGreen-7},
    };
} 