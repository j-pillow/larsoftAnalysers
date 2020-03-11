#ifndef PROTODUNE_SHOWER_UTILS_H
#define PROTODUNE_SHOWER_UTILS_H

///////////////////////////////////////////////////////////////
// ProtoDUNEShowerUtils
//  - Class to help analysers access useful shower information
// 
///////////////////////////////////////////////////////////////

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "art/Framework/Principal/Event.h"

namespace protoana {

  class ProtoDUNEShowerUtils {

  public:

    ProtoDUNEShowerUtils();
    ~ProtoDUNEShowerUtils();

    /// Get the hits from a given reco shower
    const std::vector<const recob::Hit*> GetRecoShowerHits(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const;
    const std::vector< art::Ptr< recob::Hit > >  GetRecoShowerHitsArtPtr(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const;
   
    /// Get the hits from a given reco shower from a specific plane
    const std::vector<const recob::Hit*> GetRecoShowerHitsFromPlane(const recob::Shower &shower, art::Event const &evt, const std::string showerModule, unsigned int planeID) const;
   
    /// Get the number of hits from a given reco shower
    unsigned int GetNumberRecoShowerHits(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const;
   
    /// Get the associated PCAxis object (from a principal component analysis)
    std::vector<const recob::PCAxis*> GetRecoShowerPCAxis(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const;
   
    /// Estimate the energy of a shower from its hits
    std::vector<double> EstimateEnergyFromHitCharge( const std::vector<const recob::Hit*> &hits, calo::CalorimetryAlg caloAlg);
    std::vector<double> EstimateEnergyFromHitArea(   const std::vector<const recob::Hit*> &hits, calo::CalorimetryAlg caloAlg, double T0, const TVector3 showerDir );
    std::vector<double> EstimateEnergyFromHitAmp(    const std::vector<const recob::Hit*> &hits, calo::CalorimetryAlg caloAlg, double T0, const TVector3 showerDir );
   
    /// If the shower.ID() isn't filled we must find the actual shower index ourselves
    int GetShowerIndex(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const;

    /// Get shower calo info
    std::vector<anab::Calorimetry> GetRecoShowerCalorimetry(const recob::Shower &shower, art::Event const &evt, const std::string showerModule, const std::string caloModule) const;
   
    /// Calculate dEdx for showers 
    std::vector<double> GetdEdx5cm( const recob::Shower & shower, art::Event const & evt, double T0, calo::CalorimetryAlg caloAlg, const std::string showerModule, const std::string hitModule );
    
  private:


  };

}

#endif

