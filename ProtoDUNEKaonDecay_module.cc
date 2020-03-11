////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNEKaonDecay
// File:        ProtoDUNEKaonDecay_module.cc
//
// Copy of the ProtoDUNEAnalTree from Georgios Christodoulou,
// but I've cut all the stuff I didn't care about, changed all the
// arrays to use vectors, and changed the name slightly.
//
// Updated carious includes to reflect the new structure of dunetpc
// and the new protoduneana ups modules
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/MCBase/MCHitCollection.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/ArtDataHelper/MVAReader.h"

#include "dune/DuneObj/ProtoDUNEBeamEvent.h"
//#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEbeamsim.h"
//#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEBeamInstrument.h"

//#include "dune/Protodune/singlephase/DataUtils/ProtoDUNETrackUtils.h"
//#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEShowerUtils.h"
//#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
//#include "dune/Protodune/singlephase/DataUtils/ProtoDUNETruthUtils.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
//#include "protoduneana/Utilities/ProtoDUNEDataUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TGraph2D.h"
#include <Math/Vector3D.h>
#include "TLorentzVector.h"
//#include "TMath"

// C++ Includes
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#ifdef __MAKECINT__
#pragma link C++ class std::vector<double>+;
#endif

namespace protoana {
  class ProtoDUNEKaonDecay;
}

class protoana::ProtoDUNEKaonDecay : public art::EDAnalyzer
{
  public:

    explicit ProtoDUNEKaonDecay(fhicl::ParameterSet const & p);

    ProtoDUNEKaonDecay(ProtoDUNEKaonDecay const &) = delete;
    ProtoDUNEKaonDecay(ProtoDUNEKaonDecay &&) = delete;
    ProtoDUNEKaonDecay & operator = (ProtoDUNEKaonDecay const &) = delete;
    ProtoDUNEKaonDecay & operator = (ProtoDUNEKaonDecay &&) = delete;

    virtual void beginJob() override;
    virtual void endJob() override;

    // Required functions.
    void analyze(art::Event const & evt) override;


  private:

    // Typdefs for shorter initialisation
    // Use stsd:numeric_limits as it is more obvious than -999 as being incorrect.
    std::numeric_limits<int>    iLim;
    std::numeric_limits<double> dLim;

    // Backtracker stuff
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;


    // Utility functions
    protoana::ProtoDUNEDataUtils       dataUtil;
    protoana::ProtoDUNEPFParticleUtils pfpUtil;
    protoana::ProtoDUNETrackUtils      trackUtil;
    protoana::ProtoDUNEShowerUtils     showerUtil;
    protoana::ProtoDUNETruthUtils      truthUtil;
    protoana::ProtoDUNEBeamlineUtils   beamUtil;

    // Track momentum algorithm. Calculates momentum based on track range
    trkf::TrackMomentumCalculator trmom;

    // Geometry
    const spacecharge::SpaceCharge* sce = lar::providerFrom<spacecharge::SpaceChargeService>();
    geo::GeometryCore const * fGeometry = &*(art::ServiceHandle<geo::Geometry>());
    const detinfo::DetectorProperties* detprop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->provider();
    art::ServiceHandle< geo::Geometry > geom;

    // Initialise tree variables
    void InitialiseBeamInfo();
    void InitialiseRecoInfo();
    void InitialiseMCInfo();
    void InitialiseMatchedMcInfo();
    void InitialiseAnaTree();

    bool EventSelector( art::Event const & evt,
                        std::vector< const recob::PFParticle* > & PrimaryPFParticles );

    void FillPFParticle( art::Event const & evt,
                         art::ValidHandle< std::vector< recob::PFParticle > > & recoParticles,
                         const recob::PFParticle* particle,
                         int flowLevel );

    void FillShower(  art::Event const & evt,
                      const recob::Shower* thisShower,
                      const recob::PFParticle* particle );

    void FillBeamData( art::Event const & evt );

    void FillConfigTree();


//    void FillSPVecs(  art::Event const & evt,
//                      std::vector< art::Ptr<recob::Hit> > hits );

    const std::vector<const recob::Hit*> GetRecoShowerHitsFromPlane(  const recob::Shower &shower,
                                                                      art::Event const &evt,
                                                                      unsigned int planeID ) const;

    std::vector<double> EstimateEnergyPerHitModBox( const std::vector<const recob::Hit*> &hits ) ;

    void showerEnergyEstimationMethods( art::Event const & evt,
                                        const recob::Shower* thisShower );

    void missingEnergy(  art::Event const & evt,
                         const std::vector< art::Ptr<recob::Hit> > thisShower );

    double dEdx5cm( art::Event const & evt,
                  const recob::Shower* thisShower,
                  std::vector<double> xVec,
                  std::vector<double> yVec,
                  std::vector<double> zVec,
                  std::vector<double> dEdxVec,
                  std::vector<double> pitchVec );

    double MydEdx(  std::vector<double> xVec,
                    std::vector<double> yVec,
                    std::vector<double> zVec,
                    std::vector<double> eVec );

    void oldProjLengthFinder( std::vector<double> xVec, 
                              std::vector<double> yVec, 
                              std::vector<double> zVec );

    double transToPCA(  TVector3 point,
                        TVector3 cent,
                        TVector3 axis );

    double GetYPos( double z,
                  std::vector<TVector3> spointVec );

    std::pair<TVector3, double> xyzSpatialAndPitch( art::Ptr<recob::Hit> hit,
                                                    const recob::Shower* thisShower,
                                                    TVector3 spoint );

    std::pair<double,double> dQcorrections( TVector3 pos,
                                            double dQ );


    void sandbox( art::Event const & evt,
                  std::vector<art::Ptr<recob::Hit>> showerHits,
                  std::vector<art::Ptr<recob::Hit>> mcHits,
                  std::vector<art::Ptr<recob::Hit>> allHits,
                  const simb::MCParticle* mcParticle );

    void showerStartFinder(   std::vector<double> xVec,
                              std::vector<double> yVec,
                              std::vector<double> zVec,
                              std::vector<double> eVec );

    double startFinder( std::vector< double > xyBinStds, double threshold );

    std::vector<double> xyBinStdCalc( std::vector< std::vector<double> > binStds );

    std::vector<double> makeBinEdges( double binWidth, double maxZ );

    std::vector< std::vector<double> > binStdCalc( std::vector< std::vector<TVector3> > binnedPoints );

    std::pair< std::vector< std::vector<TVector3> >, std::vector< std::vector<double> > > binShower(  std::vector<double> binEdges,
                                                                                                      std::vector<TVector3> xyzVec,
                                                                                                      std::vector<double> eVec );

    std::vector<double> dEdtCalc( std::vector< std::vector<double> > binnedEnergies, float binWidth );

      // Make the same typedefs as pandora except pandora::CartesianVector = TVector3
    typedef TVector3 EigenValues;
    typedef std::vector<TVector3> EigenVectors;
    typedef std::pair< const TVector3, double > WeightedPoint;
    typedef std::vector<WeightedPoint> WeightedPointVector;

    std::tuple< TVector3, TVector3, TVector3, TVector3, TVector3 > runPCA( std::vector<TVector3> pointVector, const recob::Shower* thisShower );
    std::tuple< TVector3, TVector3, TVector3, TVector3, TVector3 > recursivePCA( std::vector<TVector3> pointVector, const recob::Shower* thisShower, std::vector<double> weightVector, double prevAveWeight, int recursionLevel );



    // Constants from DUNE docDB 15974 by A Paudel.
    const double rho   = 1.383;   // g/cm^3 (LAr density at 18 psi)
    const double betap = 0.212;   // kV/cm * g/cm / MeV
    const double E0    = 0.4867;  // kV/cm (nominal electric field)
    const double alpha = 0.93;    // ArgoNeuT-determined parameter at E0 kV/cm
    const double Wion  = 23.6e-6; // ArgoNeuT-determined parameter at E0 kV/cm
    //const double recomb = 0.6417; // Recombination from Aaron
    const double recomb = 0.715; // Modal recombination from LArG4
    
    // Get the variable Efield using data driven maps.
    TFile* ef = new TFile("/data/neutrino/dune/phrzzf/DUNE_v08_40_00/prod/corrections/SCE_DataDriven_180kV_v3.root");
    TH3F* xneg = (TH3F*)ef->Get("Reco_ElecField_X_Neg");
    TH3F* yneg = (TH3F*)ef->Get("Reco_ElecField_Y_Neg");
    TH3F* zneg = (TH3F*)ef->Get("Reco_ElecField_Z_Neg");
    TH3F* xpos = (TH3F*)ef->Get("Reco_ElecField_X_Pos");
    TH3F* ypos = (TH3F*)ef->Get("Reco_ElecField_Y_Pos");
    TH3F* zpos = (TH3F*)ef->Get("Reco_ElecField_Z_Pos");

      // Get dQ/dx YZ and X correction factors.
    TFile* yzf = new TFile("/data/neutrino/dune/phrzzf/DUNE_v08_40_00/prod/corrections/YZcalo_sce.root");
    TH2F* yzCorrPosHist =(TH2F*)yzf->Get("correction_dqdx_ZvsY_positiveX_hist_2");
    TH2F* yzCorrNegHist =(TH2F*)yzf->Get("correction_dqdx_ZvsY_negativeX_hist_2");

    TFile* xf = new TFile("/data/neutrino/dune/phrzzf/DUNE_v08_40_00/prod/corrections/Xcalo_sce.root");
    TH1F* xCorrHist = (TH1F*)xf->Get("dqdx_X_correction_hist_2");


    // fcl parameters
    bool fVerbose;
    bool fFillSPVecs;
    int fPDGSelection;
    double fCalibConst;
    double fNormConst;
    const art::InputTag fBeamLabel;
    std::string fCaloLabel;
    std::string fParticleIDLabel;
    std::string fTrackLabel;
    std::string fShowerLabel;
    std::string fShowerCaloLabel;
    std::string fPFParticleLabel;
    std::string fHitLabel;
    std::string fGeneratorLabel;
    std::string fSimulationLabel;
    calo::CalorimetryAlg fCaloAlg;

    // List for all MC particles.

    // Get beam trigger info
    bool beamTriggerEvent;

    // Initial mc particle
    const simb::MCParticle* fInitialParticle;

    // Config tree
    TTree *fConfigTree;

    int fNAPAs;
    int fNChansPerAPA;
    int fNCryostats;
    int fNTPCs;
    int fNChannels;
    int fNPlanes;
    double fTPCBoundsX0;
    double fTPCBoundsX1;
    double fTPCBoundsY0;
    double fTPCBoundsY1;
    double fTPCBoundsZ0;
    double fTPCBoundsZ1;


    // Output tree
    TTree *fBeamTree;
    TTree *fAnaTree;

    // Tree variables

    // ==========================
    // ======== Run Info ========
    // ==========================

    int    fRun;
    int    fSubRun;
    int    fevent;
    int    fNactiveFEMBs0;
    int    fNactiveFEMBs1;
    int    fNactiveFEMBs2;
    int    fNactiveFEMBs3;
    int    fNactiveFEMBs4;
    int    fNactiveFEMBs5;
    double fTimeStamp;


    // ===========================
    // ======== Beam Info ========
    // ===========================

    int    fbeamTrigger;
    int    fcerenkovStatus0;
    int    fcerenkovStatus1;
    double ftof;
    double fbeamTime;
    double fbeamPosX;
    double fbeamPosY;
    double fbeamPosZ;
    double fbeamDirX;
    double fbeamDirY;
    double fbeamDirZ;
    double fbeamMomentum;
    double fcerenkovTime0;
    double fcerenkovTime1;
    double fcerenkovPressure0;
    double fcerenkovPressure1;

    // ==========================
    // ======== Ana Tree ========
    // ==========================

      // Pass Checks
    bool NEWPCAPass;
    bool NEWDEDXPass;
    int anaTree_flow;
    int anaTree_isShower;

      // Beam info
    double anaTree_beamMomentum;


      // PF Particle Info
    int anaTree_nHits;
    int anaTree_nDaughters;
    
      // dEdx Info
    double anaTree_dEdx;
    double anaTree_dEdxAaron;
    std::vector<double> anaTree_correcteddEdx;
    std::vector<double> anaTree_pitch;

      // Energy estimates
    double anaTree_energy;
    double anaTree_energyCorrected;
    double anaTree_energyAaron;

      // Energy Corrections
    double anaTree_depositionCorrection;
    double anaTree_missedHitsCorrection;
    double anaTree_noHitsCorrection;
    double anaTree_mcIDEdiscrep;
    double anaTree_contaminationCorrection;
    double anaTree_pureIDEcorrection;

      // MC Info
    double anaTree_mcInitEnergy;
    double anaTree_mcDepEnergy;

      // Energy Comparison
    double anaTree_recoMinusTrueOverTrue;

      // Misc
    double anaTree_trueEnergyOfShowerHits;

      // Charge Info
    double anaTree_totalCharge;

      // Purity and Completeness
    double anaTree_hitPurity;
    double anaTree_hitCompleteness;
    double anaTree_energyPurity;

      // Angle Info
    double anaTree_pandoraAngleToMC; 
    double anaTree_correctedAngleToMC;
    double anaTree_recursiveAngleToMC;

      // Start info
    double anaTree_showerStart;

      // PCAs
    TVector3 anaTree_pandoraEvals;
    TVector3 anaTree_pandoraAvePos;
    TVector3 anaTree_pandoraPvec;
    TVector3 anaTree_pandoraSvec;
    TVector3 anaTree_pandoraTvec;
    
    TVector3 anaTree_correctedEvals;
    TVector3 anaTree_correctedAvePos;
    TVector3 anaTree_correctedPvec;
    TVector3 anaTree_correctedSvec;
    TVector3 anaTree_correctedTvec;
   
    TVector3 anaTree_recursiveEvals;
    TVector3 anaTree_recursiveAvePos;
    TVector3 anaTree_recursivePvec;
    TVector3 anaTree_recursiveSvec;
    TVector3 anaTree_recursiveTvec;

      // Lengths
        // Eigen value lengths
    double anaTree_pandoraPriEigenValLength;
    double anaTree_pandoraSecEigenValLength;
    double anaTree_pandoraTerEigenValLength;

    double anaTree_correctedPriEigenValLength;
    double anaTree_correctedSecEigenValLength;
    double anaTree_correctedTerEigenValLength;
    
        // Projection lengths
    double anaTree_pandoraPriProjectionLength;
    double anaTree_pandoraSecProjectionLength;
    double anaTree_pandoraTerProjectionLength;

    double anaTree_correctedPriProjectionLength;
    double anaTree_correctedSecProjectionLength;
    double anaTree_correctedTerProjectionLength;

      // 3D Position Info
    std::vector<double> anaTree_nonCorrectedHit3DX;
    std::vector<double> anaTree_nonCorrectedHit3DY;
    std::vector<double> anaTree_nonCorrectedHit3DZ;

    std::vector<double> anaTree_correctedHit3DX;
    std::vector<double> anaTree_correctedHit3DY;
    std::vector<double> anaTree_correctedHit3DZ;

      // Per hit vectors
    std::vector< std::pair< double, std::vector<const sim::IDE*> > > anaTree_hitIDEmap;
    std::vector<double> anaTree_dEdt;
    std::vector<double> anaTree_hitCharge;
    std::vector<double> anaTree_hitTrueEnergy;

      // CNN Shower Score
    double anaTree_cnnShowerScore;


};

// =========================================================================================

protoana::ProtoDUNEKaonDecay::ProtoDUNEKaonDecay(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  dataUtil         ( p.get<fhicl::ParameterSet> ( "DataUtils"       ) ),
  beamUtil         ( p.get<fhicl::ParameterSet> ( "BeamUtils"       ) ),
  fVerbose         ( p.get<bool>                ( "Verbose"         ) ),
  fFillSPVecs      ( p.get<bool>                ( "FillSPVecs"      ) ),
  fPDGSelection    ( p.get<int>                 ( "PDGSelection"    ) ),
  fCalibConst      ( p.get<double>              ( "CalibConst"      ) ),
  fNormConst       ( p.get<double>              ( "NormConst"       ) ),
  fBeamLabel       ( p.get<art::InputTag>       ( "BeamLabel"       ) ),
  fCaloLabel       ( p.get<std::string>         ( "CaloLabel"       ) ),
  fParticleIDLabel ( p.get<std::string>         ( "ParticleIDLabel" ) ),
  fTrackLabel      ( p.get<std::string>         ( "TrackLabel"      ) ),
  fShowerLabel     ( p.get<std::string>         ( "ShowerLabel"     ) ),
  fShowerCaloLabel ( p.get<std::string>         ( "ShowerCaloLabel" ) ),
  fPFParticleLabel ( p.get<std::string>         ( "PFParticleLabel" ) ),
  fHitLabel        ( p.get<std::string>         ( "HitLabel"        ) ),
  fGeneratorLabel  ( p.get<std::string>         ( "GeneratorLabel"  ) ),
  fSimulationLabel ( p.get<std::string>         ( "SimulationLabel" ) ),
  fCaloAlg         ( p.get<fhicl::ParameterSet> ( "CalorimetryAlg"  ) )
{

}

// =========================================================================================

void protoana::ProtoDUNEKaonDecay::beginJob()
{

    // Init the root trees
  art::ServiceHandle<art::TFileService> tfs;

  // Config Tree
  fConfigTree = tfs->make<TTree>("Config", "Configuration Tree");

  fConfigTree->Branch("NAPAs",             &fNAPAs        );
  fConfigTree->Branch("NChansPerAPA",      &fNChansPerAPA );
  fConfigTree->Branch("NCryostats",        &fNCryostats   );
  fConfigTree->Branch("NTPCs",             &fNTPCs        );
  fConfigTree->Branch("NChannels",         &fNChannels    );
  fConfigTree->Branch("NPlanes",           &fNPlanes      );
  fConfigTree->Branch("ActiveTPCBoundsX0", &fTPCBoundsX0  );
  fConfigTree->Branch("ActiveTPCBoundsX1", &fTPCBoundsX1  );
  fConfigTree->Branch("ActiveTPCBoundsY0", &fTPCBoundsY0  );
  fConfigTree->Branch("ActiveTPCBoundsY1", &fTPCBoundsY1  );
  fConfigTree->Branch("ActiveTPCBoundsZ0", &fTPCBoundsZ0  );
  fConfigTree->Branch("ActiveTPCBoundsZ1", &fTPCBoundsZ1  );
  FillConfigTree();

  // Beam info tree
  fBeamTree = tfs->make<TTree>("BeamTree", "Beam information");

  fBeamTree->Branch( "run",               &fRun               );
  fBeamTree->Branch( "subrun",            &fSubRun            );
  fBeamTree->Branch( "event",             &fevent             );
  fBeamTree->Branch( "timestamp",         &fTimeStamp         );
  fBeamTree->Branch( "NactiveFEMBs0",     &fNactiveFEMBs0     );
  fBeamTree->Branch( "NactiveFEMBs1",     &fNactiveFEMBs1     );
  fBeamTree->Branch( "NactiveFEMBs2",     &fNactiveFEMBs2     );
  fBeamTree->Branch( "NactiveFEMBs3",     &fNactiveFEMBs3     );
  fBeamTree->Branch( "NactiveFEMBs4",     &fNactiveFEMBs4     );
  fBeamTree->Branch( "NactiveFEMBs5",     &fNactiveFEMBs5     );
  fBeamTree->Branch( "beamTrigger",       &fbeamTrigger       );
  fBeamTree->Branch( "tof",               &ftof               );
  fBeamTree->Branch( "beamMomentum",      &fbeamMomentum      );
  fBeamTree->Branch( "beamTime",          &fbeamTime          );
  fBeamTree->Branch( "cerenkovStatus0",   &fcerenkovStatus0   );
  fBeamTree->Branch( "cerenkovStatus1",   &fcerenkovStatus1   );
  fBeamTree->Branch( "cerenkovTime0",     &fcerenkovTime0     );
  fBeamTree->Branch( "cerenkovTime1",     &fcerenkovTime1     );
  fBeamTree->Branch( "cerenkovPressure0", &fcerenkovPressure0 );
  fBeamTree->Branch( "cerenkovPressure1", &fcerenkovPressure1 );
  fBeamTree->Branch( "beamPosX",          &fbeamPosX          );
  fBeamTree->Branch( "beamPosY",          &fbeamPosY          );
  fBeamTree->Branch( "beamPosZ",          &fbeamPosZ          );
  fBeamTree->Branch( "beamDirX",          &fbeamDirX          );
  fBeamTree->Branch( "beamDirY",          &fbeamDirY          );
  fBeamTree->Branch( "beamDirZ",          &fbeamDirZ          );

  // Ana info tree
  fAnaTree = tfs->make<TTree>("AnaTree", "Analysis");

  fAnaTree->Branch( "nHits",                        &anaTree_nHits                        );
  fAnaTree->Branch( "nDaughters",                   &anaTree_nDaughters                   );
  fAnaTree->Branch( "dEdx",                         &anaTree_dEdx                         );
  fAnaTree->Branch( "dEdxAaron",                    &anaTree_dEdxAaron                    );
  fAnaTree->Branch( "energy",                       &anaTree_energy                       );
  fAnaTree->Branch( "energyCorrected",              &anaTree_energyCorrected              );
  fAnaTree->Branch( "energyAaron",                  &anaTree_energyAaron                  );
  fAnaTree->Branch( "depositionCorrection",         &anaTree_depositionCorrection         );
  fAnaTree->Branch( "missedHitsCorrection",         &anaTree_missedHitsCorrection         );
  fAnaTree->Branch( "noHitsCorrection",             &anaTree_noHitsCorrection             );
  fAnaTree->Branch( "mcIDEdiscrep",                 &anaTree_mcIDEdiscrep                 );
  fAnaTree->Branch( "contaminationCorrection",      &anaTree_contaminationCorrection      );
  fAnaTree->Branch( "pureIDEcorrection",            &anaTree_pureIDEcorrection            );
  fAnaTree->Branch( "mcInitEnergy",                 &anaTree_mcInitEnergy                 );
  fAnaTree->Branch( "beamMomentum",                 &anaTree_beamMomentum                 );
  fAnaTree->Branch( "mcDepEnergy",                  &anaTree_mcDepEnergy                  );
  fAnaTree->Branch( "recoMinusTrueOverTrue",        &anaTree_recoMinusTrueOverTrue        );
  fAnaTree->Branch( "trueEnergyOfShowerHits",       &anaTree_trueEnergyOfShowerHits       );
  fAnaTree->Branch( "totalCharge",                  &anaTree_totalCharge                  );
  fAnaTree->Branch( "hitPurity",                    &anaTree_hitPurity                    );
  fAnaTree->Branch( "hitCompleteness",              &anaTree_hitCompleteness              );
  fAnaTree->Branch( "energyPurity",                 &anaTree_energyPurity                 );
  fAnaTree->Branch( "pandoraAngleToMC",             &anaTree_pandoraAngleToMC             );
  fAnaTree->Branch( "correctedAngleToMC",           &anaTree_correctedAngleToMC           );
  fAnaTree->Branch( "recursiveAngleToMC",           &anaTree_recursiveAngleToMC           );
  fAnaTree->Branch( "showerStart",                  &anaTree_showerStart                  );
  fAnaTree->Branch( "pandoraEvals",                 &anaTree_pandoraEvals                 );
  fAnaTree->Branch( "pandoraAvePos",                &anaTree_pandoraAvePos                );
  fAnaTree->Branch( "pandoraPvec",                  &anaTree_pandoraPvec                  );
  fAnaTree->Branch( "pandoraSvec",                  &anaTree_pandoraSvec                  );
  fAnaTree->Branch( "pandoraTvec",                  &anaTree_pandoraTvec                  );
  fAnaTree->Branch( "correctedEvals",               &anaTree_correctedEvals               );
  fAnaTree->Branch( "correctedAvePos",              &anaTree_correctedAvePos              );
  fAnaTree->Branch( "correctedPvec",                &anaTree_correctedPvec                );
  fAnaTree->Branch( "correctedSvec",                &anaTree_correctedSvec                );
  fAnaTree->Branch( "correctedTvec",                &anaTree_correctedTvec                );
  fAnaTree->Branch( "recursiveEvals",               &anaTree_recursiveEvals               );
  fAnaTree->Branch( "recursiveAvePos",              &anaTree_recursiveAvePos              );
  fAnaTree->Branch( "recursivePvec",                &anaTree_recursivePvec                );
  fAnaTree->Branch( "recursiveSvec",                &anaTree_recursiveSvec                );
  fAnaTree->Branch( "recursiveTvec",                &anaTree_recursiveTvec                );
  fAnaTree->Branch( "pandoraPriEigenValLength",     &anaTree_pandoraPriEigenValLength     );
  fAnaTree->Branch( "pandoraSecEigenValLength",     &anaTree_pandoraSecEigenValLength     );
  fAnaTree->Branch( "pandoraTerEigenValLength",     &anaTree_pandoraTerEigenValLength     );
  fAnaTree->Branch( "correctedPriEigenValLength",   &anaTree_correctedPriEigenValLength   );
  fAnaTree->Branch( "correctedSecEigenValLength",   &anaTree_correctedSecEigenValLength   );
  fAnaTree->Branch( "correctedTerEigenValLength",   &anaTree_correctedTerEigenValLength   );
  fAnaTree->Branch( "pandoraPriProjectionLength",   &anaTree_pandoraPriProjectionLength   );
  fAnaTree->Branch( "pandoraSecProjectionLength",   &anaTree_pandoraSecProjectionLength   );
  fAnaTree->Branch( "pandoraTerProjectionLength",   &anaTree_pandoraTerProjectionLength   );
  fAnaTree->Branch( "correctedPriProjectionLength", &anaTree_correctedPriProjectionLength );
  fAnaTree->Branch( "correctedSecProjectionLength", &anaTree_correctedSecProjectionLength );
  fAnaTree->Branch( "correctedTerProjectionLength", &anaTree_correctedTerProjectionLength );
  fAnaTree->Branch( "nonCorrectedHit3DX",           &anaTree_nonCorrectedHit3DX           );
  fAnaTree->Branch( "nonCorrectedHit3DY",           &anaTree_nonCorrectedHit3DY           );
  fAnaTree->Branch( "nonCorrectedHit3DZ",           &anaTree_nonCorrectedHit3DZ           );
  fAnaTree->Branch( "correctedHit3DX",              &anaTree_correctedHit3DX              );
  fAnaTree->Branch( "correctedHit3DY",              &anaTree_correctedHit3DY              );
  fAnaTree->Branch( "correctedHit3DZ",              &anaTree_correctedHit3DZ              );
  fAnaTree->Branch( "dEdt",                         &anaTree_dEdt                         );
  fAnaTree->Branch( "hitCharge",                    &anaTree_hitCharge                    );
  fAnaTree->Branch( "hitTrueEnergy",                &anaTree_hitTrueEnergy                );
  fAnaTree->Branch( "cnnShowerScore",               &anaTree_cnnShowerScore               );

}

// =========================================================================================

void protoana::ProtoDUNEKaonDecay::analyze( art::Event const & evt )
{

    // Initialise Ana tree
  InitialiseAnaTree();

    // Lets get all the PFParticles.
  auto recoParticles = evt.getValidHandle< std::vector< recob::PFParticle > >( fPFParticleLabel );
    // Get all the primary PFParticles from the beam slice
  std::vector< const recob::PFParticle* > PrimaryPFParticles = pfpUtil.GetPFParticlesFromBeamSlice( evt, fPFParticleLabel );

  bool doEvent = EventSelector( evt, PrimaryPFParticles );

  if ( doEvent ) {
    // Initialise the tree variables
    InitialiseBeamInfo();

    // Find the relevent run and event info
    fRun    = evt.run();
    fSubRun = evt.subRun();
    fevent  = evt.id().event();
    art::Timestamp ts = evt.time();
    if (ts.timeHigh() == 0) {
      TTimeStamp ts2( ts.timeLow() );
      fTimeStamp = ts2.AsDouble();
    }
    else {
      TTimeStamp ts2( ts.timeHigh(), ts.timeLow() );
      fTimeStamp = ts2.AsDouble();
    }

    // Number of active FEMBs - All in MC, maybe not all in data
    if ( !evt.isRealData() ) { // MC
      fNactiveFEMBs0 = 20;
      fNactiveFEMBs1 = 20;
      fNactiveFEMBs2 = 20;
      fNactiveFEMBs3 = 20;
      fNactiveFEMBs4 = 20;
      fNactiveFEMBs5 = 20;
    }
    else { // Data
      std::vector< double > tempFEMBs;
      for ( size_t i{0}; i < 6; i++ ) {
        tempFEMBs.push_back( dataUtil.GetNActiveFembsForAPA(evt, i) );
      }
      fNactiveFEMBs0 = tempFEMBs[0];
      fNactiveFEMBs1 = tempFEMBs[1];
      fNactiveFEMBs2 = tempFEMBs[2];
      fNactiveFEMBs3 = tempFEMBs[3];
      fNactiveFEMBs4 = tempFEMBs[4];
      fNactiveFEMBs5 = tempFEMBs[5];
    }

    if ( !evt.isRealData() ) { // MC Info

        // Firstly we need to get the list of MCTruth objects from the generator. 
        // The standard protoDUNE simulation has fGeneratorLabel = "generator"
      auto mcTruths = evt.getValidHandle< std::vector< simb::MCTruth > >(fGeneratorLabel);

        // mcTruths is basically a pointer to an std::vector of simb::MCTruth objects. 
        // There should only be one of these, so we pass the first element into the function to get the good particle
      fInitialParticle = truthUtil.GetGeantGoodParticle( (*mcTruths)[0], evt );

        // Fill the beam trigger despite this being MC, because I like to live recklessly
      fbeamTrigger  = 12;
        // Fill the beam tree
      fBeamTree->Fill();

    }

    else FillBeamData( evt ); // Data info

    for ( const recob::PFParticle* particle : PrimaryPFParticles ) { // Reco info
      std::cout << "Fill PFParticle" << std::endl;
      FillPFParticle( evt, recoParticles ,particle, 0 );
    }
  }

  if ( (doEvent) && (!evt.isRealData()) ) {
    
      // My python ana here
    std::cout << anaTree_nHits << " | " << anaTree_flow << " | " << anaTree_dEdx << " | " << fInitialParticle->E() << " | " << anaTree_isShower << std::endl;
    if ( (NEWDEDXPass) && (NEWPCAPass) && (anaTree_nHits > 0) && (anaTree_flow == 0) && (anaTree_dEdx > -1) && (fInitialParticle->E() < 10) && (fInitialParticle->E() > 0) && (anaTree_isShower == 1) ) {
     
        // MC Energy info 
      anaTree_mcInitEnergy = fInitialParticle->E();
      anaTree_mcDepEnergy  = truthUtil.GetDepEnergyMC( evt, fGeometry, fInitialParticle->TrackId(), 2);
      anaTree_depositionCorrection = anaTree_mcInitEnergy - anaTree_mcDepEnergy/1000;

      double sumIDEVector = 0;
      double contam       = 0;
      for ( auto pair : anaTree_hitIDEmap ) {
        double hitTrueEnergy_temp = 0;
        for ( auto ide : pair.second ) {
          anaTree_trueEnergyOfShowerHits += ide->energy;
          hitTrueEnergy_temp += ide->energy;
          sumIDEVector += ( std::abs(ide->trackID) == fInitialParticle->TrackId() ) ? ide->energy : 0;
          contam       += ( std::abs(ide->trackID) != fInitialParticle->TrackId() ) ? ide->energy : 0;
        }
        ////anaTree_hitTrueEnergy.push_back(hitTrueEnergy_temp);
      }


      anaTree_energyPurity = sumIDEVector/(sumIDEVector+contam);


      anaTree_pureIDEcorrection = (sumIDEVector - anaTree_pureIDEcorrection)/1000;
      anaTree_contaminationCorrection = anaTree_contaminationCorrection/1000;
      anaTree_energyCorrected         = anaTree_energy  + anaTree_pureIDEcorrection + anaTree_missedHitsCorrection + anaTree_noHitsCorrection + anaTree_depositionCorrection - anaTree_contaminationCorrection + anaTree_mcIDEdiscrep;
      anaTree_recoMinusTrueOverTrue   = (anaTree_energyCorrected - anaTree_mcInitEnergy)/anaTree_mcInitEnergy;

        // Calculate the angle between reconstructed and mc particles
        // First tracjectory point inside TPC
      int startTPC = truthUtil.GetFirstTrajectoryPointInTPCActiveVolume(*fInitialParticle, fTPCBoundsX0, fTPCBoundsX1, fTPCBoundsY0, fTPCBoundsY1, fTPCBoundsZ0, fTPCBoundsZ1);

        // MC direction
      TVector3 mcDirTVect(fInitialParticle->Px(startTPC),fInitialParticle->Py(startTPC),fInitialParticle->Pz(startTPC));
      mcDirTVect = mcDirTVect.Unit();

        // Normalise all the PCA primary axes
      TVector3 pandoraPvecUnit = anaTree_pandoraPvec.Unit();
      TVector3 correctedPvecUnit = anaTree_correctedPvec.Unit();
      TVector3 recursivePvecUnit = anaTree_recursivePvec.Unit();

        // Find the angles
      anaTree_pandoraAngleToMC   = pandoraPvecUnit.Angle( mcDirTVect );
      anaTree_correctedAngleToMC = correctedPvecUnit.Angle( mcDirTVect );
      anaTree_recursiveAngleToMC = recursivePvecUnit.Angle( mcDirTVect );


        // Lengths
      anaTree_pandoraPriEigenValLength = 6*std::sqrt(anaTree_pandoraEvals.X());
      anaTree_pandoraSecEigenValLength = 6*std::sqrt(anaTree_pandoraEvals.Y());
      anaTree_pandoraTerEigenValLength = 6*std::sqrt(anaTree_pandoraEvals.Z());
      
      anaTree_correctedPriEigenValLength = 6*std::sqrt(anaTree_correctedEvals.X());
      anaTree_correctedSecEigenValLength = 6*std::sqrt(anaTree_correctedEvals.Y());
      anaTree_correctedTerEigenValLength = 6*std::sqrt(anaTree_correctedEvals.Z());

      oldProjLengthFinder( anaTree_nonCorrectedHit3DX, anaTree_nonCorrectedHit3DY, anaTree_nonCorrectedHit3DZ );

      fAnaTree->Fill();
    }
  }
  else if ( (doEvent) && (evt.isRealData()) ) {
    
    std::cout << anaTree_nHits << " | " << anaTree_flow << " | " << anaTree_dEdx << " | " << anaTree_beamMomentum << " | " << anaTree_isShower << std::endl;
    if ( (NEWDEDXPass) && (NEWPCAPass) && (anaTree_nHits > 0) && (anaTree_flow == 0) && (anaTree_dEdx > -1) && (anaTree_beamMomentum < 10) && (anaTree_beamMomentum > 0) && (anaTree_isShower == 1) ) {
      
        // Lengths
      anaTree_pandoraPriEigenValLength = 6*std::sqrt(anaTree_pandoraEvals.X());
      anaTree_pandoraSecEigenValLength = 6*std::sqrt(anaTree_pandoraEvals.Y());
      anaTree_pandoraTerEigenValLength = 6*std::sqrt(anaTree_pandoraEvals.Z());
      
      anaTree_correctedPriEigenValLength = 6*std::sqrt(anaTree_correctedEvals.X());
      anaTree_correctedSecEigenValLength = 6*std::sqrt(anaTree_correctedEvals.Y());
      anaTree_correctedTerEigenValLength = 6*std::sqrt(anaTree_correctedEvals.Z());

      oldProjLengthFinder( anaTree_nonCorrectedHit3DX, anaTree_nonCorrectedHit3DY, anaTree_nonCorrectedHit3DZ );
    
      fAnaTree->Fill();
    } 

  }



}
// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

bool protoana::ProtoDUNEKaonDecay::EventSelector( art::Event const & evt, std::vector< const recob::PFParticle* > & PrimaryPFParticles )
{ // Let's have some logic to select our events

  std::cout << "EventSelector" << std::endl;

  // Set a default of skipping the event
  bool doEvent = false;
  auto recoParticles = evt.getValidHandle< std::vector< recob::PFParticle > >( fPFParticleLabel );

  // Check first if we have any PF particles
  if ( PrimaryPFParticles.size() > 0 ) {
    doEvent = true;
  }

  std::cout << "PrimaryPFParticles.size(): " << PrimaryPFParticles.size() << std::endl;

  if ( doEvent ) {

    // If MC - Check if there is a good geant particle
    //       - Check if pdg code is 321
    if ( !evt.isRealData() ) { // MC
      // Firstly we need to get the list of MCTruth objects from the generator.
      // The standard protoDUNE simulation has fGeneratorLabel = "generator"
      auto mcTruths = evt.getValidHandle< std::vector< simb::MCTruth > >(fGeneratorLabel);
      // mcTruths is basically a pointer to an std::vector of simb::MCTruth objects.
      // There should only be one of these, so we pass the first element into the function to get the good particle
      const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle( (*mcTruths)[0], evt );

      // Check for good geant particle
      if ( geantGoodParticle != nullptr ) {
        // Check for pdg code == 11
        int pdgCode = geantGoodParticle->PdgCode();
        std::cout << "geantGoodParticle->PdgCode(): " << pdgCode << std::endl;
        std::cout << "fPDGSelection: " << fPDGSelection << std::endl;
          if ( std::abs(pdgCode) == fPDGSelection ) doEvent = true;
        else doEvent = false;
      }
      else {
        doEvent = false;
      }

    }

      // If data - Check if it is a beam trigger
    else { // Data
      
      if( !beamUtil.IsGoodBeamlineTrigger( evt ) ){
        std::cout << "Failed quality check" << std::endl;
        return false;
      }

      std::cout << "Passed quality check!" << std::endl << std::endl;
      
        //Access the Beam Event
      auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("beamevent");

      std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
      if( beamHandle.isValid()){
        art::fill_ptr_vector(beamVec, beamHandle);
      }

      const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0)); //Should just have one
        //Check the quality of the event
      std::cout << "Timing Trigger: " << beamEvent.GetTimingTrigger() << std::endl;
      std::cout << "Is Matched: "     << beamEvent.CheckIsMatched() << std::endl << std::endl;

      std::vector<double> momenta;
      int nMomenta;
      //int c0;
      //int c1;
      //double tof;
      //int chan;

        //Access momentum
      const std::vector< double > & the_momenta = beamEvent.GetRecoBeamMomenta();
      std::cout << "Number of reconstructed momenta: " << the_momenta.size() << std::endl;

      if( the_momenta.size() > 0 )
        std::cout << "Measured Momentum: " << the_momenta[0] << std::endl;

      if( the_momenta.size()  == 1)
        anaTree_mcInitEnergy = the_momenta[0];

      momenta.insert( momenta.end(), the_momenta.begin(), the_momenta.end() );
      /*
         for( size_t i = 0; i < the_momenta.size(); ++i ){
         momenta.push_back( the_momenta[i] );
         }
         */
      nMomenta = momenta.size();
      std::cout << "nMomenta: " << nMomenta << std::endl;
      if ( nMomenta == 0 ){
        std::cout << "Skipping event" << std::endl;
        return false;
      }
      /////////////////////////////////////////////////////////////


      std::cout << "Current: " << beamEvent.GetMagnetCurrent() << std::endl;

        //Access time of flight
      const std::vector< double > & the_tofs  = beamEvent.GetTOFs();
      const std::vector< int    > & the_chans = beamEvent.GetTOFChans();

      std::cout << "Number of measured TOF: " << the_tofs.size() << std::endl;
      std::cout << "First TOF: "              << beamEvent.GetTOF()         << std::endl;
      std::cout << "First TOF Channel: "      << beamEvent.GetTOFChan()     << std::endl << std::endl;

      std::cout << "All (TOF, Channels): " << std::endl;
      for( size_t i = 0; i < the_tofs.size(); ++i ){
        std::cout << "\t(" << the_tofs[i] << ", " << the_chans[i] << ")" << std::endl;
      }
      std::cout << std::endl;

      //if( the_tofs.size() > 0){
      //  tof = the_tofs[0];
      //  chan = the_chans[0];
      //}

        //Access Cerenkov info
      std::cout << "Cerenkov status, pressure:" << std::endl;
      //c0 = beamEvent.GetCKov0Status();
      //c1 = beamEvent.GetCKov1Status();
      std::cout << "C0: " << beamEvent.GetCKov0Status() << ", " << beamEvent.GetCKov0Pressure() << std::endl;
      std::cout << "C1: " << beamEvent.GetCKov1Status() << ", " << beamEvent.GetCKov1Pressure() << std::endl << std::endl;

        //Access PID
      std::vector< int > pids = beamUtil.GetPID( beamEvent, 1. );

      std::cout << "Possible particles" << std::endl;

      for( size_t i = 0; i < pids.size(); ++i ){
        std::cout << pids[i] << std::endl;
      }
      std::cout << std::endl;

      PossibleParticleCands candidates = beamUtil.GetPIDCandidates( beamEvent, 1. );
      std::cout << std::left << std::setw(10) << "electron " << candidates.electron << std::endl;
      std::cout << std::left << std::setw(10) << "muon "     << candidates.muon     << std::endl;
      std::cout << std::left << std::setw(10) << "pion "     << candidates.pion     << std::endl;
      std::cout << std::left << std::setw(10) << "kaon "     << candidates.kaon     << std::endl;
      std::cout << std::left << std::setw(10) << "proton "   << candidates.proton   << std::endl << std::endl;

      std::string candidates_string = beamUtil.GetPIDCandidates( beamEvent, 1. );
      std::cout << candidates_string << std::endl;

      if (candidates.kaon) {
        std::cout << "Kaon candidate" << std::endl;
        std::cout << "Doing event" << std::endl;

        anaTree_beamMomentum = the_momenta[0];
        return true;
      }
    }
  }

  // Let's see what Pandora has reconstructed our primary PFParticles as
  // I want to select primary showers, and then primary tracks with a single primary daughter shower.
  if ( doEvent ) {
    auto primary = PrimaryPFParticles[0];
    if ( primary->NumDaughters() > 0 ) { // Only want kaons that have decayed
      return true;
    }
  }

  return doEvent;
}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

void protoana::ProtoDUNEKaonDecay::FillPFParticle( art::Event const & evt,
                                                 art::ValidHandle< std::vector< recob::PFParticle > > & recoParticles,
                                                 const recob::PFParticle* particle,
                                                 int flowLevel )
{

  // Is this particle a track or a shower?
  const recob::Track*  thisTrack  = pfpUtil.GetPFParticleTrack(  *particle, evt, fPFParticleLabel, fTrackLabel );
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower( *particle, evt, fPFParticleLabel, fShowerLabel  );

  // Let's record the point in the particle flow
  anaTree_flow = flowLevel;

  // Get information from recob::Tracks/Showers
  if ( thisTrack != nullptr ) { // Tracks
    std::cout << "PFParticle is track like" << std::endl;
    anaTree_isShower = 0;
  }
  else if ( thisShower != nullptr ) { // Showers
    std::cout << "PFParticle is shower like" << std::endl;
    anaTree_isShower = 1;
    FillShower( evt, thisShower, particle );
  }
  else { // Neither
    std::cout << "PFParticle is neither like" << std::endl;
    anaTree_isShower = 0;
  }
}


// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

void protoana::ProtoDUNEKaonDecay::FillShower( art::Event const & evt, const recob::Shower* thisShower, const recob::PFParticle* particle )
{

  //if ( !evt.isRealData() ) { // MC
  //  fInitialParticle = truthUtil.GetMCParticleFromRecoShower( *thisShower, evt, fShowerLabel );
  //}

  anaTree_nHits = GetRecoShowerHitsFromPlane( *thisShower, evt, 2 ).size();

    // Get the showers' PCA info
  auto showerPCAxis = showerUtil.GetRecoShowerPCAxis( *thisShower, evt, fShowerLabel);
  if ( showerPCAxis.size() > 0 ) {

      // Eigen Values
    anaTree_pandoraEvals.SetX(  ((showerPCAxis[0])->getEigenValues())[0] );
    anaTree_pandoraEvals.SetY(  ((showerPCAxis[0])->getEigenValues())[1] );
    anaTree_pandoraEvals.SetZ(  ((showerPCAxis[0])->getEigenValues())[2] );

      // Centroid
    anaTree_pandoraAvePos.SetX( ((showerPCAxis[0])->getAvePosition())[0] );
    anaTree_pandoraAvePos.SetY( ((showerPCAxis[0])->getAvePosition())[1] );
    anaTree_pandoraAvePos.SetZ( ((showerPCAxis[0])->getAvePosition())[2] );

    std::vector< std::vector< double > > eVecs = (showerPCAxis[0])->getEigenVectors();

      // Primary Axis
    anaTree_pandoraPvec.SetX( eVecs[0][0] );
    anaTree_pandoraPvec.SetY( eVecs[0][1] );
    anaTree_pandoraPvec.SetZ( eVecs[0][2] );

      // Secondary Axis
    anaTree_pandoraSvec.SetX( eVecs[1][0] );
    anaTree_pandoraSvec.SetY( eVecs[1][1] );
    anaTree_pandoraSvec.SetZ( eVecs[1][2] );
      
      // Tertiary Axis
    anaTree_pandoraTvec.SetX( eVecs[2][0] );
    anaTree_pandoraTvec.SetY( eVecs[2][1] );
    anaTree_pandoraTvec.SetZ( eVecs[2][2] );
  }

  showerEnergyEstimationMethods( evt, thisShower );


}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

double protoana::ProtoDUNEKaonDecay::dEdx5cm( art::Event const & evt, const recob::Shower* thisShower, std::vector<double> xVec, std::vector<double> yVec, std::vector<double> zVec, std::vector<double> dEdxVec, std::vector<double> pitchVec )
{


    // Get the projections of the space points onto the primary axis of the shower
  std::vector<double> projectionsVec;
  for ( size_t i = 0 ; i < xVec.size() ; i++ ) {
    TVector3 thisPoint(xVec[i], yVec[i], zVec[i]);
    projectionsVec.push_back( thisPoint.Dot(anaTree_recursiveEvals) );
  }

    // Find space point that is at the start of the shower
  //int minProjIndex = std::min_element(projectionsVec.begin(),projectionsVec.end()) - projectionsVec.begin();
  double minProj   = *std::min_element(projectionsVec.begin(), projectionsVec.end());

  double dEdx5cm = 0;
  double meandEdx = 0;
  int count = 0;
    // Shift the projections and select the dEdx of the ones within 5cm of the first space point
  for ( size_t i = 0 ; i < projectionsVec.size() ; i++ ) {
    double newProj = projectionsVec[i] - minProj;
    if ( newProj <= 5 ) {
      dEdx5cm += (dEdxVec[i] * pitchVec[i] );
      meandEdx += dEdxVec[i];
      count++;
    }
  }

  return meandEdx/count;
}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

double protoana::ProtoDUNEKaonDecay::MydEdx( std::vector<double> xVec, std::vector<double> yVec, std::vector<double> zVec, std::vector<double> eVec )
{

  std::vector<double> xVecTrans;
  std::vector<double> yVecTrans;
  std::vector<double> zVecTrans;

    // Transform the spacepoints onto the axis of the shower - z axis = primary axis | y axis = secondary axis | x axis = tertiary axis
  for ( size_t i = 0 ; i < xVec.size() ; i++ ) {
    TVector3 point( xVec[i], yVec[i], zVec[i] );
    xVecTrans.push_back( transToPCA(point, anaTree_recursiveAvePos, anaTree_recursiveTvec) );
    yVecTrans.push_back( transToPCA(point, anaTree_recursiveAvePos, anaTree_recursiveSvec) );
    zVecTrans.push_back( transToPCA(point, anaTree_recursiveAvePos, anaTree_recursivePvec) );
  }

    // Shift everything so first point is z=0
  double shift = *std::min_element( zVecTrans.begin(), zVecTrans.end() );
  std::vector<TVector3> xyzVec;
  for ( size_t i = 0 ; i < xVecTrans.size() ; i ++ ) {
    TVector3 point( xVecTrans[i], yVecTrans[i], zVecTrans[i] - shift );
    xyzVec.push_back( point );
  }
  
    // Get the length from the projections
  double maxPrimary   = -5000.0;
  double maxSecondary = -5000.0;
  double maxTertiary  = -5000.0;
  for ( auto point : xyzVec ) {
    maxPrimary   = point.Z() > maxPrimary   ? point.Z() : maxPrimary;  
    maxSecondary = std::abs(point.Y()) > maxSecondary ? std::abs(point.Y()) : maxSecondary;  
    maxTertiary  = std::abs(point.X()) > maxTertiary  ? std::abs(point.X()) : maxTertiary;  
  }
  anaTree_correctedPriProjectionLength = maxPrimary;
  anaTree_correctedSecProjectionLength = maxSecondary;
  anaTree_correctedTerProjectionLength = maxTertiary;

    // fa
  double zExtent5cm = -1;
  double energy     = 0;
  int numberInRange = 0;
  for ( size_t i = 0 ; i < xyzVec.size() ; i++ ) {
    if ( (xyzVec[i].Z() <= 5.0) && ( std::sqrt( std::pow(xyzVec[i].X(),2) + std::pow(xyzVec[i].Y(),2) ) < 5.0 ) ) {
      numberInRange++;
      energy += eVec[i];
      zExtent5cm = xyzVec[i].Z() > zExtent5cm ? xyzVec[i].Z() : zExtent5cm;
    }
  }

    // Make sure that the number of 3D points in the dEdx region is at least 5
  if (numberInRange > 4){
    NEWDEDXPass = true;
  }
  else {
    NEWDEDXPass = false;
  }

  //return energy/zExtent5cm;
  return energy/5.0;

}

void protoana::ProtoDUNEKaonDecay::oldProjLengthFinder( std::vector<double> xVec, std::vector<double> yVec, std::vector<double> zVec ) 
{ // Get max projection lengths of the old spacepoints on the old axes
  std::vector<double> xVecTrans;
  std::vector<double> yVecTrans;
  std::vector<double> zVecTrans;

    // Transform the spacepoints onto the axis of the shower - z axis = primary axis | y axis = secondary axis | x axis = tertiary axis
  for ( size_t i = 0 ; i < xVec.size() ; i++ ) {
    TVector3 point( xVec[i], yVec[i], zVec[i] );
    xVecTrans.push_back( transToPCA(point, anaTree_pandoraAvePos, anaTree_pandoraTvec) );
    yVecTrans.push_back( transToPCA(point, anaTree_pandoraAvePos, anaTree_pandoraSvec) );
    zVecTrans.push_back( transToPCA(point, anaTree_pandoraAvePos, anaTree_pandoraPvec) );
  }

  double maxPrimary   = -5000.0;
  double maxSecondary = -5000.0;
  double maxTertiary  = -5000.0;
    // Shift everything so first point is z=0
  double shift = *std::min_element( zVecTrans.begin(), zVecTrans.end() );
  std::vector<TVector3> xyzVec;
  for ( size_t i = 0 ; i < xVecTrans.size() ; i ++ ) {
    maxPrimary   = zVecTrans[i] - shift   > maxPrimary   ? zVecTrans[i] - shift   : maxPrimary;
    maxSecondary = std::abs(yVecTrans[i]) > maxSecondary ? std::abs(yVecTrans[i]) : maxSecondary;
    maxTertiary  = std::abs(xVecTrans[i]) > maxTertiary  ? std::abs(xVecTrans[i]) : maxTertiary;
    TVector3 point( xVecTrans[i], yVecTrans[i], zVecTrans[i] - shift );
    xyzVec.push_back( point );
  }

  anaTree_pandoraPriProjectionLength = maxPrimary;
  anaTree_pandoraSecProjectionLength = maxSecondary;
  anaTree_pandoraTerProjectionLength = maxTertiary;

}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

double protoana::ProtoDUNEKaonDecay::transToPCA( TVector3 point, TVector3 cent, TVector3 axis )
{
    // Vector for cent -> point
  TVector3 pVec = point - cent;

    // Projection along axis
  axis = axis.Unit();
  double proj = pVec.Dot(axis);

  return proj;
}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

void protoana::ProtoDUNEKaonDecay::showerStartFinder( std::vector<double> xVec, std::vector<double> yVec, std::vector<double> zVec, std::vector<double> eVec )
{

  double X0 = 14;//cm - radiation length in LAr
  double Rm = 10;//cm - moliere radius in LAr

  std::vector<double> xVecTrans;
  std::vector<double> yVecTrans;
  std::vector<double> zVecTrans;

    // Transform the spacepoints onto the axis of the shower - z axis = primary axis | y axis = secondary axis | x axis = tertiary axis
  for ( size_t i = 0 ; i < xVec.size() ; i++ ) {
    TVector3 point( xVec[i], yVec[i], zVec[i] );
    xVecTrans.push_back( transToPCA(point, anaTree_recursiveAvePos, anaTree_recursiveTvec) );
    yVecTrans.push_back( transToPCA(point, anaTree_recursiveAvePos, anaTree_recursiveSvec) );
    zVecTrans.push_back( transToPCA(point, anaTree_recursiveAvePos, anaTree_recursivePvec) );
  }

    // Shift everything so first point is z=0
  double shift = *std::min_element( zVecTrans.begin(), zVecTrans.end() );
  std::vector<TVector3> xyzVec;
  for ( size_t i = 0 ; i < xVecTrans.size() ; i ++ ) {
    TVector3 point( xVecTrans[i]/Rm, yVecTrans[i]/Rm, (zVecTrans[i] - shift)/X0 );
    xyzVec.push_back( point );
  }

    // Define the bins for the shower to be binned into
  float binWidth = 1/X0;
  double maxZ = *std::max_element( zVecTrans.begin(), zVecTrans.end() );
  std::vector<double> binEdges = makeBinEdges( binWidth, maxZ );

  if ( binEdges.size() == 0 ) { // Shower is smaller than bin width - so we probably don't care about it really.
    anaTree_showerStart = 0;
    return;
  }

    // Bin the shower
  auto [binnedPoints, binnedEnergies] = binShower( binEdges, xyzVec, eVec );

    // Bin stds
  std::vector< std::vector<double> > binStds = binStdCalc( binnedPoints );

    // Get the combined secondary/tertiary bin stds
  std::vector<double> xyBinStds = xyBinStdCalc( binStds );

    // Find the start of the shower
  double threshold = 0.005; // Should make a fcl param
  int startBinIndex = startFinder( xyBinStds, threshold );

    // Now we get the mean z position of the bin that is the start of the shower
  anaTree_showerStart = binStds[startBinIndex][3];

  std::vector<double> dEdt = dEdtCalc( binnedEnergies, binWidth );
  anaTree_dEdt = dEdt;

}

std::vector<double> protoana::ProtoDUNEKaonDecay::makeBinEdges( double binWidth, double maxZ )
{   // Function to create bins of width binWidth up to the bin including maxZ

  std::vector<double> binEdges;            // Create empty bin vector
  double BE = 0;                           // First edge is always 0
  while ( BE < maxZ ) {                    // Loop up till maxZ
    binEdges.push_back(BE);                // Push the bin edge to the vector
    BE += binWidth;                        // Increase the bin edge value by the bin width
    if (BE > maxZ) binEdges.push_back(BE); // Upper edge of final bin
  }
  return binEdges;

}

std::pair< std::vector< std::vector<TVector3> >, std::vector< std::vector<double> > > protoana::ProtoDUNEKaonDecay::binShower( std::vector<double> binEdges, std::vector<TVector3> xyzVec, std::vector<double> eVec )
{   // Function to bin the coordinates of the shower into predefined bins.

  std::vector< std::vector<TVector3> > bins(  binEdges.size() - 1 ); // Create the number of bins we have
  std::vector< std::vector<double>   > eBins( binEdges.size() - 1 );

  int index = 0;
  for ( auto point : xyzVec ) {                                 // Go through the points
    double z  = point.Z();                                      // Get the Z position
    for ( size_t i = 0 ; i < (binEdges.size() - 1) ; i++ ) {    // Go through the bins to see which one point belongs in
      if ( (binEdges[i] <= z) && (z < binEdges[i+1]) ) {        // If z is in the range of these bin edges
        bins[i].push_back(point);                               // Then add this point to that bin
        eBins[i].push_back(eVec[index]);
        break;                                                  // No point checking any more bin edges
      }
    }
    index++;
  }

  std::pair< std::vector< std::vector<TVector3> >, std::vector< std::vector<double> > > poop(bins,eBins);
  return poop;
}

std::vector< std::vector<double> > protoana::ProtoDUNEKaonDecay::binStdCalc( std::vector< std::vector<TVector3> > binnedPoints )
{   // Function that finds the std in the x,z,y for each bin. Also store the mean of z.

  std::vector< std::vector<double> > stds(binnedPoints.size());

  int counter = 0;
  for ( std::vector<TVector3> bin : binnedPoints ) {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    for ( TVector3 point : bin ) {
      x.push_back(point.X());
      y.push_back(point.Y());
      z.push_back(point.Z());
    }
    stds[counter].push_back(TMath::StdDev(x.begin(),x.end()));
    stds[counter].push_back(TMath::StdDev(y.begin(),y.end()));
    stds[counter].push_back(TMath::StdDev(z.begin(),z.end()));
    stds[counter].push_back(TMath::Mean(z.begin(),z.end()));
    counter++;
  }
  return stds;
}

std::vector<double> protoana::ProtoDUNEKaonDecay::xyBinStdCalc( std::vector< std::vector<double> > binStds )
{   // Function that finds the cross-sectional area of a bin in the secondary-tertiary plane.

  std::vector< double > xyBinStds;

  for ( std::vector<double> bin : binStds ) {
    xyBinStds.push_back(bin[0]*bin[1]);
  }
  return xyBinStds;
}

double protoana::ProtoDUNEKaonDecay::startFinder( std::vector< double > xyBinStds, double threshold )
{   // Function that finds the bin number that the shower starts in

    // If the difference in the cross-sectional area of one bin to the next is bigger than the threshold value
    // call that bin the start of the shower.
  for ( size_t i = 1 ; i < xyBinStds.size() ; i++ ) {
    if ( (xyBinStds[i] - xyBinStds[i-1]) > threshold ) return i;
  }
  std::cout << "Shower start threshold not met. Return zeroth bin" << std::endl;
  return 0; // In the case where the threshold wasn't met
}

std::vector<double> protoana::ProtoDUNEKaonDecay::dEdtCalc( std::vector< std::vector<double> > binnedEnergies, float binWidth )
{

  std::vector<double> bindEdt;
  for ( auto bin : binnedEnergies ) {
    double eSum = 0;
    for ( double elem : bin ) {
      eSum += elem;
    }
    bindEdt.push_back( eSum/binWidth );
  }

  return bindEdt;

}

std::tuple< TVector3, TVector3, TVector3, TVector3, TVector3 > protoana::ProtoDUNEKaonDecay::runPCA( std::vector<TVector3> pointVector, const recob::Shower* thisShower )
{ // Pretty much copied from Pandora

    // Make a weighted vector in case we want weights
  WeightedPointVector weightedPointVector;

  for ( const auto &point : pointVector ) {
    weightedPointVector.push_back( std::make_pair( point, 1. ) );
    //weightedPointVector.push_back( std::make_pair( point, 1./point.Z() ) );
  }

    // Get the centroid position of the points
  double meanPosition[3] = {0., 0., 0.};
  double sumWeight(0.);

  for ( const WeightedPoint &weightedPoint : weightedPointVector ) {

    const TVector3 &point(weightedPoint.first);
    const double weight(weightedPoint.second);

    meanPosition[0] += static_cast<double>(point.X()) * weight;
    meanPosition[1] += static_cast<double>(point.Y()) * weight;
    meanPosition[2] += static_cast<double>(point.Z()) * weight;
    sumWeight += weight;

  }

  meanPosition[0] /= sumWeight;
  meanPosition[1] /= sumWeight;
  meanPosition[2] /= sumWeight;
  TVector3 centroid( meanPosition[0], meanPosition[1], meanPosition[2] );

    // Define covariance matrix
  double xi2(0.);
  double xiyi(0.);
  double xizi(0.);
  double yi2(0.);
  double yizi(0.);
  double zi2(0.);

  for ( const WeightedPoint &weightedPoint : weightedPointVector ) {

    const TVector3 &point(weightedPoint.first);
    const double weight(weightedPoint.second);
    const double x(static_cast<double>((point.X()) - meanPosition[0]));
    const double y(static_cast<double>((point.Y()) - meanPosition[1]));
    const double z(static_cast<double>((point.Z()) - meanPosition[2]));

    xi2  += x * x * weight;
    xiyi += x * y * weight;
    xizi += x * z * weight;
    yi2  += y * y * weight;
    yizi += y * z * weight;
    zi2  += z * z * weight;
  }

    // Use Eigen
  Eigen::Matrix3f sig;

  sig <<  xi2, xiyi, xizi,
          xiyi, yi2, yizi,
          xizi, yizi, zi2;

  sig *= 1. / sumWeight;

  Eigen::SelfAdjointEigenSolver< Eigen::Matrix3f> eigenMat(sig);

  if ( eigenMat.info() != Eigen::ComputationInfo::Success ) {
    std::cout << "Something done gone wrong with the PCA gurl!" << std::endl;
    std::cout << "There appears to be a problem with the decomposition." << std::endl;
    std::cout << "Will throw away this event" << std::endl;
    NEWPCAPass = false;
    TVector3 nullVec(0,0,0);
    return std::tuple< TVector3, TVector3, TVector3, TVector3, TVector3 >{ nullVec, nullVec, nullVec, nullVec, nullVec };
  }

  typedef std::pair<float, size_t> EigenValColPair;
  typedef std::vector<EigenValColPair> EigenValColVector;

  EigenValColVector eigenValColVector;
  const auto &resultEigenMat( eigenMat.eigenvalues() );
  eigenValColVector.emplace_back( resultEigenMat(0), 0 );
  eigenValColVector.emplace_back( resultEigenMat(1), 1 );
  eigenValColVector.emplace_back( resultEigenMat(2), 2 );

  std::sort(eigenValColVector.begin(), eigenValColVector.end(), [](const EigenValColPair &left, const EigenValColPair &right){return left.first > right.first;});

    // Get the eigen values
  EigenValues eigenVals(eigenValColVector.at(0).first, eigenValColVector.at(1).first, eigenValColVector.at(2).first);

    // Principal axes
  const Eigen::Matrix3f &eigenVecs(eigenMat.eigenvectors());

  EigenVectors eigenVectors;
  for (const EigenValColPair &pair : eigenValColVector) {
    eigenVectors.emplace_back(eigenVecs(0, pair.second), eigenVecs(1, pair.second), eigenVecs(2, pair.second));
  }

    // Ensure always pointing downstream
  const float testProjection( eigenVectors.at(0).Dot( thisShower->ShowerStart() - centroid ) );
  const float directionScaleFactor((testProjection > std::numeric_limits<float>::epsilon()) ? -1.f : 1.f);

  const TVector3 primaryAxis(eigenVectors.at(0) * directionScaleFactor);
  const TVector3 secondaryAxis(eigenVectors.at(1) * directionScaleFactor);
  const TVector3 tertiaryAxis(eigenVectors.at(2) * directionScaleFactor);

  NEWPCAPass = true;
  return std::tuple< TVector3, TVector3, TVector3, TVector3, TVector3 >{ eigenVals, primaryAxis, secondaryAxis, tertiaryAxis, centroid };

}

std::tuple< TVector3, TVector3, TVector3, TVector3, TVector3 > protoana::ProtoDUNEKaonDecay::recursivePCA( std::vector<TVector3> pointVector, const recob::Shower* thisShower, std::vector<double> weightVector, double prevAveWeight, int recursionLevel )
{ // Pretty much copied from Pandora

    // Make a weighted vector in case we want weights
  WeightedPointVector weightedPointVector;

  for ( size_t i = 0 ; i < pointVector.size() ; i++ ) {
    weightedPointVector.push_back( std::make_pair( pointVector[i], weightVector[i] ) );
  }

    // Get the centroid position of the points
  double meanPosition[3] = {0., 0., 0.};
  double sumWeight(0.);

  for ( const WeightedPoint &weightedPoint : weightedPointVector ) {

    const TVector3 &point(weightedPoint.first);
    const double weight(weightedPoint.second);

    meanPosition[0] += static_cast<double>(point.X()) * weight;
    meanPosition[1] += static_cast<double>(point.Y()) * weight;
    meanPosition[2] += static_cast<double>(point.Z()) * weight;
    sumWeight += weight;

  }

  meanPosition[0] /= sumWeight;
  meanPosition[1] /= sumWeight;
  meanPosition[2] /= sumWeight;
  TVector3 centroid( meanPosition[0], meanPosition[1], meanPosition[2] );

    // Define covariance matrix
  double xi2(0.);
  double xiyi(0.);
  double xizi(0.);
  double yi2(0.);
  double yizi(0.);
  double zi2(0.);

  for ( const WeightedPoint &weightedPoint : weightedPointVector ) {

    const TVector3 &point(weightedPoint.first);
    const double weight(weightedPoint.second);
    const double x(static_cast<double>((point.X()) - meanPosition[0]));
    const double y(static_cast<double>((point.Y()) - meanPosition[1]));
    const double z(static_cast<double>((point.Z()) - meanPosition[2]));

    xi2  += x * x * weight;
    xiyi += x * y * weight;
    xizi += x * z * weight;
    yi2  += y * y * weight;
    yizi += y * z * weight;
    zi2  += z * z * weight;
  }

    // Use Eigen
  Eigen::Matrix3f sig;

  sig <<  xi2, xiyi, xizi,
          xiyi, yi2, yizi,
          xizi, yizi, zi2;

  sig *= 1. / sumWeight;

  Eigen::SelfAdjointEigenSolver< Eigen::Matrix3f> eigenMat(sig);

  if ( eigenMat.info() != Eigen::ComputationInfo::Success ) {
    std::cout << "Something done gone wrong with the PCA gurl!" << std::endl;
    std::cout << "There appears to be a problem with the decomposition." << std::endl;
    std::cout << "Will throw away this event" << std::endl;
    NEWPCAPass = false;
    TVector3 nullVec(0,0,0);
    return std::tuple< TVector3, TVector3, TVector3, TVector3, TVector3 >{ nullVec, nullVec, nullVec, nullVec, nullVec };
  }

  typedef std::pair<float, size_t> EigenValColPair;
  typedef std::vector<EigenValColPair> EigenValColVector;

  EigenValColVector eigenValColVector;
  const auto &resultEigenMat( eigenMat.eigenvalues() );
  eigenValColVector.emplace_back( resultEigenMat(0), 0 );
  eigenValColVector.emplace_back( resultEigenMat(1), 1 );
  eigenValColVector.emplace_back( resultEigenMat(2), 2 );

  std::sort(eigenValColVector.begin(), eigenValColVector.end(), [](const EigenValColPair &left, const EigenValColPair &right){return left.first > right.first;});

    // Get the eigen values
  EigenValues eigenVals(eigenValColVector.at(0).first, eigenValColVector.at(1).first, eigenValColVector.at(2).first);

    // Principal axes
  const Eigen::Matrix3f &eigenVecs(eigenMat.eigenvectors());

  EigenVectors eigenVectors;
  for (const EigenValColPair &pair : eigenValColVector) {
    eigenVectors.emplace_back(eigenVecs(0, pair.second), eigenVecs(1, pair.second), eigenVecs(2, pair.second));
  }

    // Ensure always pointing downstream
  const float testProjection( eigenVectors.at(0).Dot( thisShower->ShowerStart() - centroid ) );
  const float directionScaleFactor((testProjection > std::numeric_limits<float>::epsilon()) ? 1.f : -1.f);

  const TVector3 primaryAxis(eigenVectors.at(0) * directionScaleFactor);
  const TVector3 secondaryAxis(eigenVectors.at(1) * directionScaleFactor);
  const TVector3 tertiaryAxis(eigenVectors.at(2) * directionScaleFactor);


  weightVector.clear();
  TVector3 start = centroid - (3*std::sqrt(eigenVals.X()))*primaryAxis;
  TVector3 end   = centroid + (3*std::sqrt(eigenVals.X()))*primaryAxis;
  for ( const auto point : pointVector ) {
    TVector3 pointToStart = point - start;
    TVector3 pointToEnd   = point - end;
    TVector3 startToEnd   = end - start;
    TVector3 crossProd    = pointToStart.Cross(pointToEnd);
    double distToLine     = crossProd.Mag() / startToEnd.Mag(); 
    weightVector.push_back( 1./(distToLine*point.Z()) );
  }
  double maxWeight = *std::max_element(weightVector.begin(), weightVector.end());
  double aveWeight = 0;
  for ( size_t i = 0 ; i < weightVector.size() ; i++ ) {
    weightVector[i] /= maxWeight;
    aveWeight += weightVector[i];
  }
  aveWeight /= weightVector.size();
  
  if ( (std::abs(aveWeight - prevAveWeight) < 0.00001) || (recursionLevel == 50) ) {
    NEWPCAPass = true;
    return std::tuple< TVector3, TVector3, TVector3, TVector3, TVector3 >{ eigenVals, primaryAxis, secondaryAxis, tertiaryAxis, centroid };
  }
  else {
    recursionLevel++;
    std::cout << "Recursion number " << recursionLevel << std::endl;
    return recursivePCA( pointVector, thisShower, weightVector, aveWeight, recursionLevel );
  }



}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

void protoana::ProtoDUNEKaonDecay::showerEnergyEstimationMethods( art::Event const & evt, const recob::Shower* thisShower )
{

  int showerID = showerUtil.GetShowerIndex( *thisShower, evt, fShowerLabel );

    // CNN scores
  anab::MVAReader<recob::Hit,4> hitResults(evt, "emtrkmichelid:emtrkmichel" );
  int cnnCount = 0;

    // Let's get the hits from our shower
  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(fShowerLabel);
  art::FindManyP<recob::Hit> findHits(recoShowers,evt,fShowerLabel);
  const std::vector<art::Ptr<recob::Hit>> hits = findHits.at(showerID);

    // Let's make some variables for storing the energy estimations
  double aaronEstimate = 0;
  double jamesEstimate = 0;

    // We'll need a vector of the space points so we can find average Y positions of hits that don't have spacepoints
  std::vector<TVector3> spointVec;
  std::vector<double>   eVec;

    // Go through the hits.
    // We'll have to do this twice, as we need to know the Y-positons of hits with spacepoints
    // so we can then apply corrections to hits without spacepoints.
  for ( art::Ptr<recob::Hit> hit : hits ) {

      // Skip hits not on the collection plane
    if ( hit->WireID().Plane != 2 ) continue;

    std::array<float,4> cnn_out = hitResults.getOutput( hit );
    //double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ];
    //double cnn_score = cnn_out[ hitResults.getIndex("em") ]/p_trk_or_sh;
    anaTree_cnnShowerScore += cnn_out[ hitResults.getIndex("em") ];
    cnnCount++;

      // Get the spacepoints for this hit
    auto recoHits = evt.getValidHandle< std::vector< recob::Hit > >(fHitLabel);
    art::FindManyP< recob::SpacePoint > HitToSpoints( recoHits, evt, "pandora" );
    std::vector< art::Ptr<recob::SpacePoint> > SpacePoints = HitToSpoints.at(hit.key());

      // Skip hits without spacepoints for now
    if ( SpacePoints.size() == 0 ) continue;

      // Loop through the spacepoints. Ideally there should only be one.
    for ( art::Ptr<recob::SpacePoint> spoint : SpacePoints ) {

      double dQ = hit->Integral();

      double x = spoint->XYZ()[0];
      double y = spoint->XYZ()[1];
      double z = spoint->XYZ()[2];

        // Record the non corrected position
      anaTree_nonCorrectedHit3DX.push_back( x );
      anaTree_nonCorrectedHit3DY.push_back( y );
      anaTree_nonCorrectedHit3DZ.push_back( z );

        // Apply spacial corrections to x,y,z positions
      TVector3 tempVec(x,y,z);
      auto [thisSpoint, pitch] = xyzSpatialAndPitch( hit, thisShower, tempVec );

        // Put the correct space point in the vector for later
      spointVec.push_back(thisSpoint);

        // Record the corrected position
      anaTree_correctedHit3DX.push_back( thisSpoint.X() );
      anaTree_correctedHit3DY.push_back( thisSpoint.Y() );
      anaTree_correctedHit3DZ.push_back( thisSpoint.Z() );

        // Apply corrections to the dQ of the hit
      auto [ dQCorr, Ef ] = dQcorrections( thisSpoint, dQ );

        // Collect all the charge
      anaTree_totalCharge += dQCorr;
      anaTree_hitCharge.push_back(dQCorr);

        // Find the dQ/dx
      double dQdxCorr       = dQCorr/pitch;

        // Scale the dQ/dx for mod box
      double dQdxCorrScaled = dQdxCorr/fCalibConst;

        // Calculate the dE/dx from the mod box equation
      double dEdxCorr = ( exp( dQdxCorrScaled*( betap / (rho * Ef) * Wion) ) - alpha ) / ( betap / (rho * Ef) );

        // Pusht the dE/dx and the pitch into vectors for use with finding the dE/dx at the begining of the shower
        // via Aaron's method
      anaTree_correcteddEdx.push_back(  dEdxCorr   );
      anaTree_pitch.push_back( pitch );

        // Estimate the energy by mine and Aaron's method
        // Physically, there is no diffference. Mine just goes straight to GeV.
      aaronEstimate += (dQCorr * Wion) / (fCalibConst * recomb);    // See https://indico.fnal.gov/event/22724/contribution/1/material/slides/0.pdf
      jamesEstimate += (dQCorr) / (4.237e7 * recomb * fCalibConst); // 4.237e7 is a scaling factor to go from number of electrons -> GeV
        // Push the hit energy estimation into a vector
      double eEstimation = 1000 * (dQCorr) / (4.237e7 * recomb * fCalibConst);
      eVec.push_back(eEstimation);
    
      
      if ( !evt.isRealData() ) { // MC
          // Get the sim IDEs that created this hit.
          // Push the energy of it into a vector so we can compare our individual hit estimation to the raw ide
        std::vector< const sim::IDE* > hitIDEs = bt_serv->HitToSimIDEs_Ps( hit );
        double thisHitTrueE = 0;
        double thisHitPureE = 0;
        double thisHitContE = 0;
        for ( auto ide : hitIDEs ) {
          if ( std::abs(ide->trackID) == fInitialParticle->TrackId() ) {
            // This ide is from the mc particle
            thisHitPureE += ide->energy > -1 ? ide->energy : 0; 
            thisHitTrueE += ide->energy > -1 ? ide->energy : 0; 
          }
          else{
            // This ide is from a different mc particle
            thisHitContE += ide->energy > -1 ? ide->energy : 0; 
            thisHitTrueE += ide->energy > -1 ? ide->energy : 0; 
          }
        }
        double thisHitPurePcent = thisHitTrueE > 0 ? thisHitPureE / thisHitTrueE : 0;
        double thisHitContPcent = thisHitTrueE > 0 ? thisHitContE / thisHitTrueE : 0;
        anaTree_hitTrueEnergy.push_back(thisHitTrueE);

        anaTree_contaminationCorrection += eEstimation * thisHitContPcent;
        anaTree_pureIDEcorrection       += eEstimation * thisHitPurePcent;

        anaTree_hitIDEmap.push_back( std::pair<double,std::vector<const sim::IDE*> >(eEstimation,hitIDEs) );

        double hitTrueE = 0;
        for ( auto ide : hitIDEs ) {
          hitTrueE += ide->energy;
        }
      }
    }
  }

  double aveY = 0;
  for ( auto point : spointVec ) {
    aveY += point.Y();
  }
  aveY /= spointVec.size();

    // Now we have the positions of the hits with spacepoints we can go through the hits that don't have them
  for ( art::Ptr<recob::Hit> hit : hits ) {

      // Skip hits not on the collection plane
    if ( hit->WireID().Plane != 2 ) continue;

      // Get the spacepoints for this hit
    auto recoHits = evt.getValidHandle< std::vector< recob::Hit > >(fHitLabel);
    art::FindManyP< recob::SpacePoint > HitToSpoints( recoHits, evt, "pandora" );
    std::vector< art::Ptr<recob::SpacePoint> > SpacePoints = HitToSpoints.at(hit.key());

      // Skip hits WITH spacepoints now
    if ( SpacePoints.size() != 0 ) continue;

      // Get the dQ
    double dQ = hit->Integral();

      // Now we change up how we get x,y,z
      // x: We convert the time into an xposition
      // z: We take the wire number
      // y: We take the average of the spacepoints either side of this hit in z
    double x = detprop->ConvertTicksToX(hit->PeakTime(),hit->WireID());
    double z = fGeometry->Wire(hit->WireID()).GetCenter().Z() ;
    //double y = GetYPos(z,spointVec);
    double y = aveY;
      
      // Apply spacial corrections to x,y,z positions
    TVector3 tempVec(x,y,z);
    auto [thisSpoint, pitch] = xyzSpatialAndPitch( hit, thisShower, tempVec );

      // Apply corrections to the dQ of the hit
    auto [ dQCorr, Ef ] = dQcorrections( thisSpoint, dQ );

      // Collect all the charge
    anaTree_totalCharge += dQCorr;
    //anaTree_hitCharge.push_back(dQCorr);

      // Find the dQ/dx
    double dQdxCorr       = dQCorr/pitch;

      // Scale the dQ/dx for mod box
    double dQdxCorrScaled = dQdxCorr/fCalibConst;

      // Calculate the dE/dx from the mod box equation
    double dEdxCorr = ( exp( dQdxCorrScaled*( betap / (rho * Ef) * Wion) ) - alpha ) / ( betap / (rho * Ef) );

      // Pusht the dE/dx and the pitch into vectors for use with finding the dE/dx at the begining of the shower
      // via Aaron's method
    anaTree_correcteddEdx.push_back(  dEdxCorr   );
    anaTree_pitch.push_back( pitch );

      // Estimate the energy by mine and Aaron's method
      // Physically, there is no diffference. Mine just goes straight to GeV.
    aaronEstimate += (dQCorr * Wion) / (fCalibConst * recomb);    // See https://indico.fnal.gov/event/22724/contribution/1/material/slides/0.pdf
    jamesEstimate += (dQCorr) / (4.237e7 * recomb * fCalibConst); // 4.237e7 is a scaling factor to go from number of electrons -> GeV

      // Push the hit energy estimation into a vector
    double eEstimation = 1000 * (dQCorr) / (4.237e7 * recomb * fCalibConst);
    //eVec.push_back(eEstimation);


    if ( !evt.isRealData() ) { // MC
        // Get the sim IDE that created this hit.
        // Push the energy of it into a vector so we can compare our individual hit estimation to the raw ide
      std::vector< const sim::IDE* > hitIDEs = bt_serv->HitToSimIDEs_Ps( hit );
      double thisHitTrueE = 0;
      double thisHitPureE = 0;
      double thisHitContE = 0;
      for ( auto ide : hitIDEs ) {
        if ( std::abs(ide->trackID) == fInitialParticle->TrackId() ) {
          // This ide is from the mc particle
          thisHitPureE += ide->energy > -1 ? ide->energy : 0; 
          thisHitTrueE += ide->energy > -1 ? ide->energy : 0; 
        }
        else{
          // This ide is from a different mc particle
          thisHitContE += ide->energy > -1 ? ide->energy : 0; 
          thisHitTrueE += ide->energy > -1 ? ide->energy : 0; 
        }
      }
      double thisHitPurePcent = thisHitTrueE > 0 ? thisHitPureE / thisHitTrueE : 0;
      double thisHitContPcent = thisHitTrueE > 0 ? thisHitContE / thisHitTrueE : 0;

      anaTree_contaminationCorrection += eEstimation * thisHitContPcent;
      anaTree_pureIDEcorrection       += eEstimation * thisHitPurePcent;
      anaTree_hitIDEmap.push_back( std::pair<double,std::vector<const sim::IDE*> >(eEstimation,hitIDEs) );

      double hitTrueE = 0;
      for ( auto ide : hitIDEs ) {
        hitTrueE += ide->energy;
      }
    }
  }

    // Rerun PCA on corrected 3D positions
    // Make a vector of tvector3s. Yes I know I could have done some things here more efficiently - but I can't be bothered to change things now
  std::vector<TVector3> pointVector;
  for ( size_t i = 0 ; i < anaTree_correctedHit3DX.size() ; i++ ) {
    TVector3 vec( anaTree_correctedHit3DX[i], anaTree_correctedHit3DY[i], anaTree_correctedHit3DZ[i] );
    pointVector.push_back(vec);
  }
  auto [ eigenVals, primaryAxis, secondaryAxis, tertiaryAxis, centroid ] = runPCA( pointVector, thisShower );

  if ( NEWPCAPass ) {
      // Eigen Values
    anaTree_correctedEvals.SetXYZ( eigenVals.X(), eigenVals.Y(), eigenVals.Z() );
      // Centroid
    anaTree_correctedAvePos.SetXYZ( centroid.X(), centroid.Y(), centroid.Z() );
      // Primary Axis
    anaTree_correctedPvec.SetXYZ( primaryAxis.X(), primaryAxis.Y(), primaryAxis.Z() );
      // Secondary Axis
    anaTree_correctedSvec.SetXYZ( secondaryAxis.X(), secondaryAxis.Y(), secondaryAxis.Z() );
      // Tertiary Axis
    anaTree_correctedTvec.SetXYZ( tertiaryAxis.X(), tertiaryAxis.Y(), tertiaryAxis.Z() );

      // Calculate PCA weights
    std::vector<double> weightsVector;
    TVector3 start = anaTree_correctedAvePos - (3*std::sqrt(anaTree_correctedEvals.X()))*anaTree_correctedPvec;
    TVector3 end   = anaTree_correctedAvePos + (3*std::sqrt(anaTree_correctedEvals.X()))*anaTree_correctedPvec;
    for ( const auto point : pointVector ) {
      TVector3 pointToStart = point - start;
      TVector3 pointToEnd   = point - end;
      TVector3 startToEnd   = end - start;
      TVector3 crossProd    = pointToStart.Cross(pointToEnd);
      double distToLine     = crossProd.Mag() / startToEnd.Mag(); 
      weightsVector.push_back( 1./(distToLine*point.Z()) );
    }
    double maxWeight = *std::max_element(weightsVector.begin(), weightsVector.end());
    double aveWeight = 0;
    for ( size_t i = 0 ; i < weightsVector.size() ; i++ ) {
      weightsVector[i] /= maxWeight;
      aveWeight += weightsVector[i];
    }
    aveWeight /= weightsVector.size();

      // Rerun PCA using weights and recursion
    auto [ eigenValsRecursive, primaryAxisRecursive, secondaryAxisRecursive, tertiaryAxisRecursive, centroidRecursive ] = recursivePCA( pointVector, thisShower, weightsVector, aveWeight, 0 );

      // Eigen Values
    anaTree_recursiveEvals.SetXYZ( eigenValsRecursive.X(), eigenValsRecursive.Y(), eigenValsRecursive.Z() );
      // Centroid
    anaTree_recursiveAvePos.SetXYZ( centroidRecursive.X(), centroidRecursive.Y(), centroidRecursive.Z() );
      // Primary Axis
    anaTree_recursivePvec.SetXYZ( primaryAxisRecursive.X(), primaryAxisRecursive.Y(), primaryAxisRecursive.Z() );
      // Secondary Axis
    anaTree_recursiveSvec.SetXYZ( secondaryAxisRecursive.X(), secondaryAxisRecursive.Y(), secondaryAxisRecursive.Z() );
      // Tertiary Axis
    anaTree_recursiveTvec.SetXYZ( tertiaryAxisRecursive.X(), tertiaryAxisRecursive.Y(), tertiaryAxisRecursive.Z() );

    // =================== Calculate the dE/dx for the first 5cm of the shower ===================
      // This is Aaron's method
    anaTree_dEdxAaron = dEdx5cm( evt, thisShower, anaTree_correctedHit3DX, anaTree_correctedHit3DY, anaTree_correctedHit3DZ, anaTree_correcteddEdx, anaTree_pitch );

      // This is my geometrical method
    if ( anaTree_correctedHit3DX.size() > 5 ) {
      anaTree_dEdx = MydEdx( anaTree_correctedHit3DX, anaTree_correctedHit3DY, anaTree_correctedHit3DZ, eVec );
    }
    else {
      anaTree_dEdx = 0;
    }
    if ( anaTree_correctedHit3DX.size() > 5 ) { // Only do this if there are more than 5 spacepoints
      showerStartFinder( anaTree_correctedHit3DX, anaTree_correctedHit3DY, anaTree_correctedHit3DZ, eVec );
    }
    else {
      anaTree_showerStart = 0;
    }

      // Calculate the average shower score from the cnn
    anaTree_cnnShowerScore /= cnnCount; 
    
      // Now we put the energy estimations into the root trees
    anaTree_energyAaron = aaronEstimate/1000;
    anaTree_energy      = jamesEstimate;
    
    if ( !evt.isRealData() ) { // MC
        // Now we want to estimate how much energy we're missing.
      missingEnergy( evt, hits );
    }

  }
}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

std::pair<TVector3, double> protoana::ProtoDUNEKaonDecay::xyzSpatialAndPitch( art::Ptr<recob::Hit> hit, const recob::Shower* thisShower, TVector3 spoint )
{
  // This all comes from Aaron's shower calorimetry module in larreco. I just do it here so I don't have to loop through the calo vector
  // and then the hit vector and all that. So much easier this way.

    // Calculate the direction corrected pitch
  float wire_pitch = fGeometry->WirePitch( hit->View() );
  float pitch = wire_pitch;
  float angleToVert = fGeometry->WireAngleToVertical( hit->View(), hit->WireID() );
  angleToVert -= 0.5*::util::pi<>();
  float cosgamma = std::abs( sin(angleToVert) * thisShower->Direction().Y() + cos(angleToVert) * thisShower->Direction().Z() );
  if ( cosgamma > 0 )       pitch /= cosgamma;
  if ( pitch < wire_pitch ) pitch = wire_pitch;

  // Correct the position of the spacepoints
  geo::Vector_t posOffsets = {0., 0., 0.};
  geo::Vector_t dirOffsets = {0., 0., 0.};

  posOffsets = sce->GetCalPosOffsets(geo::Point_t(spoint),hit->WireID().TPC);

  dirOffsets = sce->GetCalPosOffsets( geo::Point_t{ spoint.X() + pitch*thisShower->Direction().X(),
                                                    spoint.Y() + pitch*thisShower->Direction().Y(),
                                                    spoint.Z() + pitch*thisShower->Direction().Z() },
                                      hit->WireID().TPC );

  TVector3 dir_corr = { pitch*thisShower->Direction().X() - dirOffsets.X() + posOffsets.X(),
                        pitch*thisShower->Direction().Y() + dirOffsets.Y() - posOffsets.Y(),
                        pitch*thisShower->Direction().Z() + dirOffsets.Z() - posOffsets.Z() };

  pitch = dir_corr.Mag();

  spoint.SetX(spoint.X() - posOffsets.X());
  spoint.SetY(spoint.Y() + posOffsets.Y());
  spoint.SetZ(spoint.Z() + posOffsets.Z());

  std::pair<TVector3,double> myPair(spoint, pitch);
  return myPair;
}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

std::pair<double,double> protoana::ProtoDUNEKaonDecay::dQcorrections( TVector3 pos, double dQ )
{

    // Unpack the TVector3
  double x = pos.X();
  double y = pos.Y();
  double z = pos.Z();

    // Get the E-field and Y-Z correction
    // These come from the root files included at the top of the file
  double Ef;
  double yzCorr;
  if ( x >= 0 ) {
    double Ex = E0 + (E0 * xpos->GetBinContent( xpos->FindBin(x,y,z) ));
    double Ey = E0 * ypos->GetBinContent( ypos->FindBin(x,y,z) );
    double Ez = E0 * zpos->GetBinContent( zpos->FindBin(x,y,z) );
    Ef     = sqrt( (Ex * Ex) + (Ey * Ey) + (Ez * Ez) );
    yzCorr = yzCorrPosHist->GetBinContent( yzCorrPosHist->FindBin(z,y) );
  }
  else {
    double Ex = E0 + (E0 * xneg->GetBinContent( xneg->FindBin(x,y,z) ));
    double Ey = E0 * yneg->GetBinContent( yneg->FindBin(x,y,z) );
    double Ez = E0 * zneg->GetBinContent( zneg->FindBin(x,y,z) );
    Ef     = sqrt( (Ex * Ex) + (Ey * Ey) + (Ez * Ez) );
    yzCorr = yzCorrNegHist->GetBinContent( yzCorrNegHist->FindBin(z,y) );
  }

    // Get the X correction
  double xCorr  = xCorrHist->GetBinContent( xCorrHist->FindBin(x) );

    // Calculate the corrected dQ
  double dQCorr = (dQ * fNormConst * xCorr * yzCorr);

    // Package up the corrected dQ and the E-field correction
  std::pair<double,double> myPair( dQCorr, Ef );
  return myPair;
}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================


void protoana::ProtoDUNEKaonDecay::missingEnergy( art::Event const & evt, const std::vector<art::Ptr<recob::Hit>> showerHits )
{


    // Lets get the hits for the initial MC particle
    // First get all the hits in the event
    // I'm not using the truth util function to do this because I need art::Ptr<recob::Hit>s so you can compare their .key()
  art::Handle< std::vector<recob::Hit> > hitHandle;
  std::vector< art::Ptr<recob::Hit> > allHits;
  if ( evt.getByLabel( fHitLabel, hitHandle ) ) art::fill_ptr_vector( allHits, hitHandle );

    // Get the initial mc particle
  const simb::MCParticle* mcParticle = fInitialParticle;

    // Match the hits in the event to the ones that backtrack to the initial particle
  std::vector< art::Ptr<recob::Hit> > mcHits;
  for ( auto hit : allHits ) {
    // bt_serv->HitToTrackIds( *(hit.get()) ) gives us a vector of track ids that this hit is linked to
    for ( const int trackId : bt_serv->HitToTrackIds( *(hit.get()) ) ) {

        // If one of the track ids is our initial particle, then we push the hit into our vector
        // We take the absolute value of the track id, as negative ids are actually from the positive id.
        // Ask Leigh for more info if you need. Sorry Leigh!!
      //if ( pi_serv->TrackIdToParticle_P(std::abs(trackId)) == pi_serv->TrackIdToParticle_P(mcParticle->TrackId()) ) {
      if ( std::abs(trackId) == mcParticle->TrackId() ) {
        mcHits.push_back(hit);
        break;
      }
    }
  }

    // Now we want to see which hits of the mc particle DIDN'T get associated to the shower
  std::vector<bool> hitInShower;
  for ( size_t i = 0 ; i < mcHits.size() ; i++ ) {

      // Set the default case
    bool hitInShowerBool = false;

      // Go through the hits in the shower
    for ( size_t j = 0 ; j < showerHits.size() ; j++ ) {

        // If the mc hit matchs a hit in the shower, then the mc hit is in the shower (obviously....)
      if (mcHits[i].key() == showerHits[j].key()) {
        hitInShower.push_back( true );
        hitInShowerBool = true;
        break;
      }
    }
      // If there is no match, then our mc hit didn't get incorporated into the shower.
      // Doesn't matter if it's in another object, we just care it isn't in our shower.
    if (!hitInShowerBool) hitInShower.push_back(false);
  }

    // Create variable for energy from hits not in the shower
  double eHitsNotInShower = 0;

    // Go through the mc hits that aren't in the shower
  for ( size_t i=0 ; i < hitInShower.size() ; i++ ) {
    if ( hitInShower[i] == true ) continue;

      // Make sure that we only look at collection plane
    if ( mcHits[i]->WireID().Plane != 2 ) continue;

      // Get the ide's and find the true energy of the hit
    std::vector< const sim::IDE* > mcIDEs = bt_serv->HitToSimIDEs_Ps(mcHits[i]);
    for ( auto ide : mcIDEs ) {
      eHitsNotInShower += ide->energy;
    }
  }

  sandbox( evt, showerHits, mcHits, allHits, mcParticle );

}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

double protoana::ProtoDUNEKaonDecay::GetYPos( double z, std::vector<TVector3> spointVec )
{

    // Find the points either side of our z position
  TVector3 lessThanZ;
  TVector3 moreThanZ;
  for ( TVector3 point : spointVec ) {
    if ( point.Z() < z ) lessThanZ = point.Z() > lessThanZ.Z() ? point : lessThanZ;
    if ( point.Z() > z ) moreThanZ = point.Z() < moreThanZ.Z() ? point : moreThanZ;
  }


    // If there is no point less than us in Z, then set Y to the point more than us in Z
  if ( ( lessThanZ.X() == lessThanZ.Y() ) && ( lessThanZ.X() == lessThanZ.Z() ) ) return moreThanZ.Y();
    // If there is no point more than us in Z, then set Y to the point less than us in Z
  if ( ( moreThanZ.X() == moreThanZ.Y() ) && ( moreThanZ.X() == moreThanZ.Z() ) ) return lessThanZ.Y();

    // Return the average Y position of these points
  return (lessThanZ.Y() + moreThanZ.Y())/2;
}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================


void protoana::ProtoDUNEKaonDecay::sandbox( art::Event const & evt, std::vector<art::Ptr<recob::Hit>> showerHits, std::vector<art::Ptr<recob::Hit>> mcHits, std::vector<art::Ptr<recob::Hit>> allHits, const simb::MCParticle* mcParticle )
{ // Use this as a testing grounds - Or, as with all things in particle physics, this is now where I calculate the missing energy stuff. *sigh*

    // Get all the plane 2 sim::IDEs for the mc particle
  std::vector< const sim::IDE* > idesFromBT = bt_serv->TrackIdToSimIDEs_Ps(mcParticle->TrackId(), fGeometry->View(2) );
  double idesFromBTenergy = 0;
  for ( auto ide : idesFromBT  ) {
    idesFromBTenergy += ide->energy;
  }

    // Get all the ides from all the hits in the event
  std::vector< const sim::IDE* > idesFromAllHits;
  for ( auto hit : allHits ) {
    auto ideVec = bt_serv->HitToSimIDEs_Ps(hit);
    for ( auto hitIDE : ideVec ) idesFromAllHits.push_back(hitIDE);
  }


    // Compare the ides from the mc particle with the ides from all the hits.
    // Seperate the ones with hits from the ones without hits
  std::vector< const sim::IDE* > idesFromBTmatchHits;
  std::vector< const sim::IDE* > idesFromBTnoMatchHits;
  double matchEnergy   = 0;
  double noMatchEnergy = 0;
  for ( auto BTide : idesFromBT ) {
    bool match = false;
    for ( auto ALLide : idesFromAllHits ) {
      if ( BTide == ALLide ) {
        matchEnergy += BTide->energy;
        idesFromBTmatchHits.push_back(BTide);
        match = true;
        break;
      }
    }
    if (!match) {
      noMatchEnergy += std::abs(BTide->trackID) == mcParticle->TrackId() ? BTide->energy : 0;
      idesFromBTnoMatchHits.push_back(BTide);
    }
  }

    // Go through all the hits from the shower and sum the ide energy from them
  double showerIDEsum = 0;
  double showerIDEsumTrackIDmatch = 0;
  for ( auto hit : showerHits ) {
    if ( hit->WireID().Plane != 2) continue;
    auto ideVec = bt_serv->HitToSimIDEs_Ps(hit);
    for ( auto ide : ideVec ) {
      showerIDEsum += ide->energy;
      if ( std::abs(ide->trackID) == mcParticle->TrackId() ) showerIDEsumTrackIDmatch += ide->energy;
    }
  }

    // Go through the hits of the mc particle and see if the ides match up with what we already have
  int nSharedHits = 0;
  double mcIDEsum = 0;
  double mcShowerMatchIDEsum = 0;
  double mcShowerNotMatchIDEsum = 0;
  for ( auto hit : mcHits ) {
    if ( hit->WireID().Plane != 2 ) continue;
    bool match = false;
    for ( auto shit : showerHits ) {
      if ( shit.key() == hit.key() ) {
        nSharedHits++;
        match = true;
      }
    }
    auto ideVec = bt_serv->HitToSimIDEs_Ps(hit);
    for ( auto ide : ideVec ) {
      if ( (std::abs(ide->trackID) != mcParticle->TrackId()) ) continue;
      mcIDEsum += ide->energy;
      if (match) mcShowerMatchIDEsum += ide->energy;
      if (!match) mcShowerNotMatchIDEsum += ide->energy;
    }
  }


  int nShowerHitsPlane2 = 0;
  int nMCHitsPlane2     = 0;
  for (auto hit : showerHits) {
    if (hit->WireID().Plane == 2) nShowerHitsPlane2++;
  }
  for (auto hit : mcHits) {
    if (hit->WireID().Plane == 2) nMCHitsPlane2++;
  }

  anaTree_hitPurity       = (double)nSharedHits / (double)nShowerHitsPlane2;
  anaTree_hitCompleteness = (double)nSharedHits / (double)nMCHitsPlane2;
  
  anaTree_missedHitsCorrection = mcShowerNotMatchIDEsum/1000;
  anaTree_noHitsCorrection     = noMatchEnergy/1000;
  anaTree_mcIDEdiscrep         = (matchEnergy - mcIDEsum)/1000;

}


// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

void protoana::ProtoDUNEKaonDecay::FillBeamData( art::Event const & evt )
{
  art::Handle< std::vector< beam::ProtoDUNEBeamEvent > > pdbeamHandle;
  std::vector< art::Ptr< beam::ProtoDUNEBeamEvent > > beamInfo;
  if ( evt.getByLabel( fBeamLabel, pdbeamHandle ) ) {
    art::fill_ptr_vector(beamInfo, pdbeamHandle);
  }

  for ( size_t i{0} ; i < beamInfo.size() ; i++ ) { // beamInfo loop

    auto & tracks = beamInfo[i]->GetBeamTracks();

    if ( !tracks.empty() ) {
      // Beam particles may have more than one corresponding particle.
      for ( size_t j{0}; j < tracks.size(); j++ ) {

        fbeamTrigger = beamInfo[i]->GetTimingTrigger();
        fbeamTime   = (double)beamInfo[i]->GetRDTimestamp();

        // TOF 0-3 means good match corresponding to the different pair-wise combinations of the upstream and downstream channels
        if ( beamInfo[i]->GetTOFChan() >= 0 ) {
          ftof = beamInfo[i]->GetTOF();
        }

        // Cerenkov info
        if ( beamInfo[i]->GetBITrigger() == 1 ) {
          fcerenkovPressure0 = beamInfo[i]->GetCKov0Pressure();
          fcerenkovPressure1 = beamInfo[i]->GetCKov1Pressure();
          fcerenkovStatus0   = beamInfo[i]->GetCKov0Status();
          fcerenkovStatus1   = beamInfo[i]->GetCKov1Status();
          fcerenkovTime0     = beamInfo[i]->GetCKov0Time();
          fcerenkovTime1     = beamInfo[i]->GetCKov1Time();
        }

        // Position info
        fbeamPosX = tracks[j].End().X();
        fbeamPosY = tracks[j].End().Y();
        fbeamPosZ = tracks[j].End().Z();
        fbeamDirX = tracks[j].StartDirection().X();
        fbeamDirY = tracks[j].StartDirection().Y();
        fbeamDirZ = tracks[j].StartDirection().Z();

        // Momenta info
        const std::vector< double > & beamMom = beamInfo[i]->GetRecoBeamMomenta();
        if ( !beamMom.empty() ) {
          fbeamMomentum = beamMom[0];
        }

        // Fill the tree here
        // This means that each entry in the tree will be a particle (if there is more than one)
        fBeamTree->Fill();

      } // track loop

    } // if !tracks

  } // End beamInfo loop

}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

void protoana::ProtoDUNEKaonDecay::FillConfigTree()
{
  fNCryostats     = fGeometry->Ncryostats();
  fNTPCs          = fGeometry->NTPC();
  fNChannels      = fGeometry->Nchannels();
  fNPlanes        = fGeometry->Nplanes();
  fNAPAs          = fGeometry->NTPC()*fGeometry->Ncryostats()/2;
  if(fNAPAs > 0) {
    fNChansPerAPA = fGeometry->Nchannels()/fNAPAs;
  }

  fTPCBoundsX0 = fTPCBoundsY0 = fTPCBoundsZ0 = 1000000.0;
  fTPCBoundsX1 = fTPCBoundsY1 = fTPCBoundsZ1 = -1000000.0;

  if(fVerbose > 0){
    std::cout << "Detector Name: " << fGeometry->DetectorName() << std::endl;
    std::cout << "GDML file: " << fGeometry->GDMLFile() << std::endl;
  }

  for (geo::TPCGeo const& TPC: fGeometry->IterateTPCs()) {
    // get center in world coordinates
    double origin[3] = {0.};
    double center[3] = {0.};
    TPC.LocalToWorld(origin, center);
    double tpc[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() };
    if(center[0] - tpc[0] < fTPCBoundsX0) fTPCBoundsX0 = center[0] - tpc[0];
    if(center[0] + tpc[0] > fTPCBoundsX1) fTPCBoundsX1 = center[0] + tpc[0];
    if(center[1] - tpc[1] < fTPCBoundsY0) fTPCBoundsY0 = center[1] - tpc[1];
    if(center[1] + tpc[1] > fTPCBoundsY1) fTPCBoundsY1 = center[1] + tpc[1];
    if(center[2] - tpc[2] < fTPCBoundsZ0) fTPCBoundsZ0 = center[2] - tpc[2];
    if(center[2] + tpc[2] > fTPCBoundsZ1) fTPCBoundsZ1 = center[2] + tpc[2];
  }

//  fConfigTree->Fill();

}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

const std::vector<const recob::Hit*> protoana::ProtoDUNEKaonDecay::GetRecoShowerHitsFromPlane(const recob::Shower &shower, art::Event const &evt, unsigned int planeID ) const
{

  std::vector<const recob::Hit*> showerHits;
  if( planeID > 2 ){
    std::cout << "Please input plane 0, 1, or 2" << std::endl;
    return showerHits;
  }

  // Shower.ID is sometimes at a default value - make sure we get the correct one
  int actualIndex = showerUtil.GetShowerIndex(shower,evt,fShowerLabel);

  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(fShowerLabel);
  art::FindManyP<recob::Hit> findHits(recoShowers,evt,fShowerLabel);
  std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(actualIndex);

  for(const art::Ptr<recob::Hit> hit : inputHits){
    unsigned int thePlane = hit.get()->WireID().asPlaneID().Plane;
    if( thePlane != planeID ) continue;

    showerHits.push_back(hit.get());

  }

  return showerHits;

}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

void protoana::ProtoDUNEKaonDecay::endJob(){

}

// ========================================================================================================================================
// ========================================================================================================================================
// ========================================================================================================================================

void protoana::ProtoDUNEKaonDecay::InitialiseBeamInfo()
{

  beamTriggerEvent = false;

  // Run and event info
  fRun           = iLim.lowest();
  fSubRun        = iLim.lowest();
  fevent         = iLim.lowest();
  fTimeStamp     = dLim.lowest();
  fNactiveFEMBs0 = iLim.lowest();
  fNactiveFEMBs1 = iLim.lowest();
  fNactiveFEMBs2 = iLim.lowest();
  fNactiveFEMBs3 = iLim.lowest();
  fNactiveFEMBs4 = iLim.lowest();
  fNactiveFEMBs5 = iLim.lowest();

  // Beam info
  fbeamTrigger       = iLim.lowest();
  fcerenkovStatus0   = iLim.lowest();
  fcerenkovStatus1   = iLim.lowest();
  ftof               = dLim.lowest();
  fbeamTime          = dLim.lowest();
  fbeamPosX          = dLim.lowest();
  fbeamPosY          = dLim.lowest();
  fbeamPosZ          = dLim.lowest();
  fbeamDirX          = dLim.lowest();
  fbeamDirY          = dLim.lowest();
  fbeamDirZ          = dLim.lowest();
  fbeamMomentum      = dLim.lowest();
  fcerenkovTime0     = dLim.lowest();
  fcerenkovTime1     = dLim.lowest();
  fcerenkovPressure0 = dLim.lowest();
  fcerenkovPressure1 = dLim.lowest();

}

// ========================================================================================================================================


void protoana::ProtoDUNEKaonDecay::InitialiseAnaTree()
{
    // Pass Checks
  NEWPCAPass  = false;
  NEWDEDXPass = false;
  anaTree_flow     = -1;
  anaTree_isShower = 0;

    // Beam info
  anaTree_beamMomentum = dLim.lowest();

    // PF Particle Info
  anaTree_nHits      = iLim.lowest();
  anaTree_nDaughters = 0;
  
    // dEdx Info
  anaTree_dEdx      = dLim.lowest();
  anaTree_dEdxAaron = dLim.lowest();
  anaTree_correcteddEdx.clear();
  anaTree_pitch.clear();

    // Energy estimates
  anaTree_energy          = dLim.lowest();
  anaTree_energyCorrected = dLim.lowest();
  anaTree_energyAaron     = dLim.lowest();

    // Energy Corrections
  anaTree_depositionCorrection    = dLim.lowest();
  anaTree_missedHitsCorrection    = dLim.lowest();
  anaTree_noHitsCorrection        = dLim.lowest();
  anaTree_mcIDEdiscrep            = dLim.lowest();
  anaTree_contaminationCorrection = 0;
  anaTree_pureIDEcorrection       = 0;

    // MC Info
  anaTree_mcInitEnergy = dLim.lowest();
  anaTree_mcDepEnergy  = dLim.lowest();

    // Energy Comparison
  anaTree_recoMinusTrueOverTrue = dLim.lowest();

    // Misc
  anaTree_trueEnergyOfShowerHits = 0;

    // Charge Info
  anaTree_totalCharge = 0;

    // Purity and Completeness
  anaTree_hitPurity       = dLim.lowest();
  anaTree_hitCompleteness = dLim.lowest();
  anaTree_energyPurity    = dLim.lowest();

    // Angle Info
  anaTree_pandoraAngleToMC   = dLim.lowest();
  anaTree_correctedAngleToMC = dLim.lowest();
  anaTree_recursiveAngleToMC = dLim.lowest();
      
    // Start info
  anaTree_showerStart = dLim.lowest();

    // PCAs
  anaTree_pandoraEvals.Clear();
  anaTree_pandoraAvePos.Clear();
  anaTree_pandoraPvec.Clear();
  anaTree_pandoraSvec.Clear();
  anaTree_pandoraTvec.Clear();
  
  anaTree_correctedEvals.Clear();
  anaTree_correctedAvePos.Clear();
  anaTree_correctedPvec.Clear();
  anaTree_correctedSvec.Clear();
  anaTree_correctedTvec.Clear();
 
  anaTree_recursiveEvals.Clear();
  anaTree_recursiveAvePos.Clear();
  anaTree_recursivePvec.Clear();
  anaTree_recursiveSvec.Clear();
  anaTree_recursiveTvec.Clear();

    // Lengths
      // Eigen value lengths
  anaTree_pandoraPriEigenValLength = dLim.lowest();
  anaTree_pandoraSecEigenValLength = dLim.lowest();
  anaTree_pandoraTerEigenValLength = dLim.lowest();

  anaTree_correctedPriEigenValLength = dLim.lowest();
  anaTree_correctedSecEigenValLength = dLim.lowest();
  anaTree_correctedTerEigenValLength = dLim.lowest();
  
      // Projection lengths
  anaTree_pandoraPriProjectionLength = dLim.lowest();
  anaTree_pandoraSecProjectionLength = dLim.lowest();
  anaTree_pandoraTerProjectionLength = dLim.lowest();

  anaTree_correctedPriProjectionLength = dLim.lowest();
  anaTree_correctedSecProjectionLength = dLim.lowest();
  anaTree_correctedTerProjectionLength = dLim.lowest();

    // 3D Position Info
  anaTree_nonCorrectedHit3DX.clear();
  anaTree_nonCorrectedHit3DY.clear();
  anaTree_nonCorrectedHit3DZ.clear();

  anaTree_correctedHit3DX.clear();
  anaTree_correctedHit3DY.clear();
  anaTree_correctedHit3DZ.clear();

    // Per hit vectors
  anaTree_hitIDEmap.clear();
  anaTree_dEdt.clear();
  anaTree_hitCharge.clear();
  anaTree_hitTrueEnergy.clear();

    // CNN Shower Score
  anaTree_cnnShowerScore = 0;
}
DEFINE_ART_MODULE(protoana::ProtoDUNEKaonDecay
   )
