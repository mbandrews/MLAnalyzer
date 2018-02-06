// -*- C++ -*-
//
// Package:    MLAnalyzer/RecHitAnalyzer
// Class:      RecHitAnalyzer
// 
//
// Original Author:  Michael Andrews
//         Created:  Sat, 14 Jan 2017 17:45:54 GMT
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DQM/HcalCommon/interface/Constants.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h" // reco::PhotonCollection defined here
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class RecHitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit RecHitAnalyzer(const edm::ParameterSet&);
    ~RecHitAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    // Tokens
    edm::EDGetTokenT<EcalRecHitCollection> EBRecHitCollectionT_;
    edm::EDGetTokenT<EBDigiCollection>     EBDigiCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> EERecHitCollectionT_;
    edm::EDGetTokenT<HBHERecHitCollection> HBHERecHitCollectionT_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollectionT_;
    edm::EDGetTokenT<reco::PhotonCollection> photonCollectionT_;
    edm::EDGetTokenT<reco::PFJetCollection> jetCollectionT_;
    edm::EDGetTokenT<reco::GenJetCollection> genJetCollectionT_;
    //edm::InputTag trackTags_; //used to select what tracks to read from configuration file

    // Declare some constants
    static const int EE_IZ_MAX = 2;
    static const int HBHE_IETA_MAX_EB = hcaldqm::constants::IETA_MAX_HB + 1;//17
    static const int HBHE_IETA_MAX_iEta20 = 20;
    static const int EE_NC_PER_ZSIDE = EEDetId::IX_MAX*EEDetId::IY_MAX; // 100*100

    // Initialize eta and phi edges
    double eta_bins_HBHE[2*(hcaldqm::constants::IETA_MAX_HE-1)+1] = 
                      {-3.000, -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305,
                       -1.218, -1.131, -1.044, -0.957, -0.870, -0.783, -0.695, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0.000, 
                        0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,  0.695,  0.783,  0.870,  0.957,  1.044,  1.131,  1.218,
                        1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,  1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  3.000}; // 57
    double eta_bins_EEm[5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB)+1] = 
                      {-3.    , -2.93  , -2.86  , -2.79  , -2.72  , -2.65  , -2.62  ,
                       -2.59  , -2.56  , -2.53  , -2.5   , -2.4644, -2.4288, -2.3932,
                       -2.3576, -2.322 , -2.292 , -2.262 , -2.232 , -2.202 , -2.172 ,
                       -2.1462, -2.1204, -2.0946, -2.0688, -2.043 , -2.0204, -1.9978,
                       -1.9752, -1.9526, -1.93  , -1.91  , -1.89  , -1.87  , -1.85  ,
                       -1.83  , -1.812 , -1.794 , -1.776 , -1.758 , -1.74  , -1.7226,
                       -1.7052, -1.6878, -1.6704, -1.653 , -1.6356, -1.6182, -1.6008,
                       -1.5834, -1.566 , -1.5486, -1.5312, -1.5138, -1.4964, -1.479 }; // 56
    double eta_bins_EEp[5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB)+1] = 
                       {1.479 ,  1.4964,  1.5138,  1.5312,  1.5486,  1.566 ,  1.5834,
                        1.6008,  1.6182,  1.6356,  1.653 ,  1.6704,  1.6878,  1.7052,
                        1.7226,  1.74  ,  1.758 ,  1.776 ,  1.794 ,  1.812 ,  1.83  ,
                        1.85  ,  1.87  ,  1.89  ,  1.91  ,  1.93  ,  1.9526,  1.9752,
                        1.9978,  2.0204,  2.043 ,  2.0688,  2.0946,  2.1204,  2.1462,
                        2.172 ,  2.202 ,  2.232 ,  2.262 ,  2.292 ,  2.322 ,  2.3576,
                        2.3932,  2.4288,  2.4644,  2.5   ,  2.53  ,  2.56  ,  2.59  ,
                        2.62  ,  2.65  ,  2.72  ,  2.79  ,  2.86  ,  2.93  ,  3.    }; // 56
    /*
    double phi_bins_EE[EBDetId::MAX_IPHI+1] = 
                       {-3.14128,  -3.12411,  -3.10677,  -3.0896 ,  -3.07225,  -3.05509,  -3.03774,  -3.02057,  -3.00323,  -2.98606,  -2.96478,  -2.94761,
                        -2.93027,  -2.9131 ,  -2.89575,  -2.87858,  -2.86124,  -2.84407,  -2.82672,  -2.80956,  -2.79221,  -2.77505,  -2.7577 ,  -2.74053,
                        -2.72319,  -2.70602,  -2.68867,  -2.67151,  -2.65416,  -2.63699,  -2.61572,  -2.59855,  -2.5812 ,  -2.56403,  -2.54668,  -2.52952,
                        -2.51217,  -2.495  ,  -2.47766,  -2.46049,  -2.44315,  -2.42598,  -2.40863,  -2.39147,  -2.37412,  -2.35696,  -2.33961,  -2.32244,
                        -2.30509,  -2.28793,  -2.26665,  -2.24948,  -2.23213,  -2.21496,  -2.19762,  -2.18045,  -2.16311,  -2.14594,  -2.12859,  -2.11143,
                        -2.09408,  -2.07691,  -2.05957,  -2.0424 ,  -2.02506,  -2.00789,  -1.99054,  -1.97338,  -1.95603,  -1.93886,  -1.91758,  -1.90041,
                        -1.88307,  -1.8659 ,  -1.84855,  -1.83138,  -1.81404,  -1.79687,  -1.77953,  -1.76236,  -1.74502,  -1.72785,  -1.7105 ,  -1.69334,
                        -1.67599,  -1.65882,  -1.64148,  -1.62431,  -1.60696,  -1.5898 ,  -1.56852,  -1.55135,  -1.534  ,  -1.51683,  -1.49949,  -1.48232,
                        -1.46497,  -1.44781,  -1.43046,  -1.41329,  -1.39595,  -1.37878,  -1.36144,  -1.34427,  -1.32692,  -1.30976,  -1.29241,  -1.27525,
                        -1.2579 ,  -1.24073,  -1.21945,  -1.20228,  -1.18494,  -1.16777,  -1.15042,  -1.13325,  -1.11591,  -1.09874,  -1.0814 ,  -1.06423,
                        -1.04688,  -1.02972,  -1.01237,  -0.99521,  -0.97786,  -0.96069,  -0.94335,  -0.92618,  -0.90883,  -0.89166,  -0.87039,  -0.85322,
                        -0.83587,  -0.8187 ,  -0.80136,  -0.78419,  -0.76684,  -0.74967,  -0.73233,  -0.71516,  -0.69782,  -0.68065,  -0.66331,  -0.64614,
                        -0.62879,  -0.61162,  -0.59428,  -0.57711,  -0.55977,  -0.54260,  -0.52132,  -0.50415,  -0.48680,  -0.46964,  -0.45229,  -0.43512,
                        -0.41778,  -0.40061,  -0.38327,  -0.36610,  -0.34875,  -0.33159,  -0.31424,  -0.29707,  -0.27973,  -0.26256,  -0.24521,  -0.22805,
                        -0.21070,  -0.19353,  -0.17226,  -0.15508,  -0.13774,  -0.12057,  -0.10322,  -0.08606,  -0.06871,  -0.05154,  -0.03420,  -0.01703,
                         0.00031,   0.01748,   0.03483,   0.05199,   0.06934,   0.08650,   0.10385,   0.12102,   0.13837,   0.15553,   0.17681,   0.19398,
                         0.21133,   0.22850,   0.24584,   0.26301,   0.28036,   0.29752,   0.31487,   0.33204,   0.34938,   0.36655,   0.38389,   0.40106,
                         0.41840,   0.43557,   0.45292,   0.47008,   0.48743,   0.5046 ,   0.52588,   0.54305,   0.56039,   0.57756,   0.59491,   0.61208,
                         0.62942,   0.64659,   0.66393,   0.68110,   0.69845,   0.71561,   0.73296,   0.75012,   0.76747,   0.78464,   0.80198,   0.81915,
                         0.83650,   0.85367,   0.87494,   0.89211,   0.90946,   0.92663,   0.94397,   0.96114,   0.97849,   0.99566,   1.013  ,   1.03017,
                         1.04751,   1.06468,   1.08202,   1.09919,   1.11654,   1.1337 ,   1.15105,   1.16822,   1.18556,   1.20273,   1.22401,   1.24118,
                         1.25853,   1.27569,   1.29304,   1.31021,   1.32755,   1.34472,   1.36207,   1.37923,   1.39658,   1.41374,   1.43109,   1.44826,
                         1.4656 ,   1.48277,   1.50012,   1.51728,   1.53463,   1.5518 ,   1.57307,   1.59025,   1.60759,   1.62476,   1.64211,   1.65927,
                         1.67662,   1.69379,   1.71113,   1.7283 ,   1.74564,   1.76281,   1.78016,   1.79732,   1.81467,   1.83183,   1.84918,   1.86635,
                         1.8837 ,   1.90086,   1.92214,   1.93931,   1.95666,   1.97383,   1.99117,   2.00834,   2.02568,   2.04285,   2.0602 ,   2.07736,
                         2.09471,   2.11188,   2.12922,   2.14639,   2.16373,   2.1809 ,   2.19825,   2.21541,   2.23276,   2.24993,   2.27121,   2.28838,
                         2.30572,   2.32289,   2.34024,   2.35741,   2.37475,   2.39192,   2.40926,   2.42643,   2.44378,   2.46094,   2.47829,   2.49545,
                         2.5128 ,   2.52997,   2.54731,   2.56448,   2.58183,   2.59899,   2.62027,   2.63744,   2.65479,   2.67196,   2.6893 ,   2.70647,
                         2.72382,   2.74098,   2.75833,   2.7755 ,   2.79284,   2.81001,   2.82735,   2.84452,   2.86187,   2.87903,   2.89638,   2.91355,
                         2.93089,   2.94806,   2.96934,   2.98651,   3.00385,   3.02102,   3.03837,   3.05554,   3.07288,   3.09005,   3.10739,   3.12456,
                         3.14191};
    */

    // Initialize Calorimeter Geometry
    const CaloGeometry* caloGeom;
    //const CaloSubdetectorGeometry* towerGeom; 
    const CaloCellGeometry* cell;

    // Initializer global cell position
    GlobalPoint pos;

    // Diagnostic histograms
    // Cumulative:
    TH2D * hEB_energy; 
    TH2D * hEB_time; 
    TH2D * hEB_adc[EcalDataFrame::MAXSAMPLES]; 
    TH2D * hEE_energy[EE_IZ_MAX]; 
    TH2D * hEE_time[EE_IZ_MAX]; 
    TH2D * hHBHE_energy; 
    TH2D * hHBHE_energy_EB; 
    //TH1D * hHBHE_depth; 
    TH1D * h_pT; 
    TH1D * h_E; 
    TH1D * h_eta; 
    TH1D * h_m0; 
    TH1D * h_leadJetPt; 
    // Single-event:
    TH2D * hEvt_HBHE_EMenergy;
    TH2D * hEvt_EE_energy[EE_IZ_MAX];
    TH2D * hEvt_HBHE_energy;

    // Main TTree
    TTree* RHTree;

    // Objects used to fill RHTree branches
    float eventId_;
    float m0_;
    float diPhoE_;
    float diPhoPt_;
    std::vector<float> vECAL_energy_;
    std::vector<float> vEB_energy_;
    std::vector<float> vEB_time_;
    std::vector<float> vEB_adc_[EcalDataFrame::MAXSAMPLES];
    std::vector<float> vEE_energy_[EE_IZ_MAX];
    std::vector<float> vEE_time_[EE_IZ_MAX];
    //std::vector<float> vHBHE_energy[hcaldqm::constants::DEPTH_NUM];
    std::vector<float> vHBHE_energy_EB_;
    std::vector<float> vHBHE_energy_;
    std::vector<float> vHBHE_EMenergy_;
    //std::vector<float> vFC_inputs_;
    math::PtEtaPhiELorentzVectorD vPho_[2];
    //float vFC_inputs_[5];
  
    // Selection and filling functions
    bool runSelections( const edm::Event&, const edm::EventSetup& );
    bool runSelections_H2GG( const edm::Event&, const edm::EventSetup& );
    bool runSelections_H24G( const edm::Event&, const edm::EventSetup& );
    void fillEBrechits( const edm::Event&, const edm::EventSetup& );
    void fillEErechits( const edm::Event&, const edm::EventSetup& );
    void fillHBHErechits( const edm::Event&, const edm::EventSetup& );
    void fillECALatHCAL();
    void fillECALstitched();
    void fillEBdigis( const edm::Event&, const edm::EventSetup& );
    //void fillFC( const edm::Event&, const edm::EventSetup& );

}; // class RecHitAnalyzer

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RecHitAnalyzer::~RecHitAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}
