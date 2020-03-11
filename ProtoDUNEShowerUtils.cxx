#include "dune/Protodune/Analysis/ProtoDUNEShowerUtils.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"

#include<TMatrixD.h>

protoana::ProtoDUNEShowerUtils::ProtoDUNEShowerUtils(){

}

protoana::ProtoDUNEShowerUtils::~ProtoDUNEShowerUtils(){

}

//=================================================================================================================================================

// Get the hits from a given reco shower
const std::vector<const recob::Hit*> protoana::ProtoDUNEShowerUtils::GetRecoShowerHits(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const{

  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(showerModule);
  art::FindManyP<recob::Hit> findHits(recoShowers,evt,showerModule);

  // Shower.ID is sometimes at a default value - make sure we get the correct one
  int actualIndex = GetShowerIndex(shower,evt,showerModule);

  std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(actualIndex);
  std::vector<const recob::Hit*> showerHits;

  for(const art::Ptr<recob::Hit> hit : inputHits){

    showerHits.push_back(hit.get());

  }

  return showerHits;  

}

//=================================================================================================================================================

// Get the hits from a given reco shower
unsigned int protoana::ProtoDUNEShowerUtils::GetNumberRecoShowerHits(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const{

  return GetRecoShowerHits(shower,evt,showerModule).size();

}

//=================================================================================================================================================

// Get the PCAxis object from the reco shower
std::vector<const recob::PCAxis*> protoana::ProtoDUNEShowerUtils::GetRecoShowerPCAxis(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const{

  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(showerModule);
  art::FindManyP<recob::PCAxis> findPCA(recoShowers,evt,showerModule);

  // Shower.ID is sometimes at a default value - make sure we get the correct one
  int actualIndex = GetShowerIndex(shower,evt,showerModule);

  std::vector<const recob::PCAxis*> pcaVec;
  for(auto const pca : findPCA.at(actualIndex)){
    pcaVec.push_back(pca.get());
  }

  return pcaVec;
}

//=================================================================================================================================================

std::vector<double> protoana::ProtoDUNEShowerUtils::EstimateEnergyFromHitCharge(const std::vector<const recob::Hit*> &hits, calo::CalorimetryAlg caloAlg) 
{
  double kGeVtoElectrons { 4.237e7 }; // obtained from utils class.. Copied for now, should use class (although this is a physical constant, so hopefully doesn't change).
  double recombination   { 1/0.63 };

  std::vector<double> showerEnergy = {0,0,0};
  
  // Find the total charge on each plane
  for ( size_t h{0} ; h < hits.size() ; h++ ) {
    const recob::Hit* hit = hits[h];
    const int plane = hit->WireID().Plane;
    showerEnergy[ plane ] += ( caloAlg.ElectronsFromADCArea( hit->Integral(), plane) * caloAlg.LifetimeCorrection(hit->PeakTime()) ) / kGeVtoElectrons;
  }
  
  showerEnergy[0] *= recombination;
  showerEnergy[1] *= recombination;
  showerEnergy[2] *= recombination;

  // caloAlg.ElectronsFromADCArea( hit->Integral(), plane) -> Does hit->Integral()/AreaConstants(plane)
  // AreaConstants(plane) is defined in calorimetry_pdune.fcl. Although that fcl file has a typo.
  // These probably need tuning for protodune data.

  return showerEnergy;
}

//=================================================================================================================================================

std::vector<double> protoana::ProtoDUNEShowerUtils::EstimateEnergyFromHitArea(const std::vector<const recob::Hit*> &hits, calo::CalorimetryAlg caloAlg, double T0, const TVector3 showerDir ) 
{
  
  art::ServiceHandle< geo::Geometry > geo;
  std::vector<double> planeEnergies = { 0, 0, 0 };

  for ( size_t h{0} ; h < hits.size() ; h++ ) {
    const recob::Hit* hit = hits[h];
    
    // Correct the pitch for the hit direction
    float wire_pitch = geo->WirePitch( hit->View() );
    float this_pitch = wire_pitch;
    float angleToVert = geo->WireAngleToVertical( hit->View(), hit->WireID() );
    angleToVert -= 0.5*::util::pi<>();
    float cosgamma = std::abs( sin(angleToVert)*showerDir.Y() + cos(angleToVert)*showerDir.Z() );
    if ( cosgamma > 0 ) this_pitch /= cosgamma;
    if ( this_pitch < wire_pitch ) this_pitch = wire_pitch;

    //std::cout << "LifeT: " << caloAlg.LifetimeCorrection(hit->PeakTime()) << std::endl;
    double dE_hit = caloAlg.dEdx_AREA( *hit, this_pitch, T0 ) * this_pitch;
    //std::cout << dE_hit << std::endl;
    if ( dE_hit < 1000 ) planeEnergies[hit->WireID().Plane] += dE_hit;
  }
   
  return planeEnergies; 
}

//=================================================================================================================================================

std::vector<double> protoana::ProtoDUNEShowerUtils::EstimateEnergyFromHitAmp(const std::vector<const recob::Hit*> &hits, calo::CalorimetryAlg caloAlg, double T0, const TVector3 showerDir ) 
{
  
  art::ServiceHandle< geo::Geometry > geo;
  std::vector<double> planeEnergies = { 0, 0, 0 };

  for ( size_t h{0} ; h < hits.size() ; h++ ) {
    const recob::Hit* hit = hits[h];
    
    // Correct the pitch for the hit direction
    float wire_pitch = geo->WirePitch( hit->View() );
    float this_pitch = wire_pitch;
    float angleToVert = geo->WireAngleToVertical( hit->View(), hit->WireID() );
    angleToVert -= 0.5*::util::pi<>();
    float cosgamma = std::abs( sin(angleToVert)*showerDir.Y() + cos(angleToVert)*showerDir.Z() );
    if ( cosgamma > 0 ) this_pitch /= cosgamma;
    if ( this_pitch < wire_pitch ) this_pitch = wire_pitch;

    double dE_hit = caloAlg.dEdx_AMP( *hit, this_pitch, T0 ) * this_pitch;
    if ( dE_hit < 1000 ) planeEnergies[hit->WireID().Plane] += dE_hit;
  }
   
  return planeEnergies; 

}


//=================================================================================================================================================

// If the shower.ID() isn't filled we must find the actual shower index ourselves
int protoana::ProtoDUNEShowerUtils::GetShowerIndex(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const{

  if(shower.ID() != -999) return shower.ID();

  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(showerModule);

  // Iterate through all showers to find the matching one to our shower
  int actualIndex = shower.ID();
  if(shower.ID() < 0){
    for(unsigned int s = 0; s < recoShowers->size(); ++s){
      const recob::Shower thisShower = (*recoShowers)[s];
      // Can't compare actual objects so look at a property
      if(fabs(thisShower.Length() - shower.Length()) < 1e-5){
        actualIndex = s;
        continue;
      }
    }
  }

  return actualIndex;

}

//=================================================================================================================================================

std::vector<anab::Calorimetry> protoana::ProtoDUNEShowerUtils::GetRecoShowerCalorimetry(const recob::Shower &shower, art::Event const &evt, const std::string showerModule, const std::string caloModule) const{

  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(showerModule);
  std::vector<anab::Calorimetry> caloInfo;
  
  try{
    const art::FindManyP<anab::Calorimetry> findCalorimetry(recoShowers,evt,caloModule);
    int actualIndex = GetShowerIndex(shower,evt,showerModule);
    std::vector<art::Ptr<anab::Calorimetry>> theseCalos = findCalorimetry.at(actualIndex);

    for( auto calo : theseCalos){
      caloInfo.push_back(*calo);
    }
  }
  catch(...){
    std::cerr << "No calorimetry object found... returning empty vector" << std::endl;
  }

  return caloInfo;


}

//=================================================================================================================================================

// Get the hits from a given reco shower
const std::vector< art::Ptr< recob::Hit > >  protoana::ProtoDUNEShowerUtils::GetRecoShowerHitsArtPtr(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const{

  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(showerModule);
  art::FindManyP<recob::Hit> findHits(recoShowers,evt,showerModule);

  // Shower.ID is sometimes at a default value - make sure we get the correct one
  int actualIndex = GetShowerIndex(shower,evt,showerModule);

  std::vector<art::Ptr<recob::Hit>> hits = findHits.at(actualIndex);

  return hits;

}

//=================================================================================================================================================

const std::vector<const recob::Hit*> protoana::ProtoDUNEShowerUtils::GetRecoShowerHitsFromPlane(const recob::Shower &shower, art::Event const &evt, const std::string showerModule, unsigned int planeID ) const{

  std::vector<const recob::Hit*> showerHits;
  if( planeID > 2 ){
    std::cout << "Please input plane 0, 1, or 2" << std::endl;
    return showerHits;
  }

  // Shower.ID is sometimes at a default value - make sure we get the correct one
  int actualIndex = GetShowerIndex(shower,evt,showerModule);

  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(showerModule);
  art::FindManyP<recob::Hit> findHits(recoShowers,evt,showerModule);
  std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(actualIndex);

  for(const art::Ptr<recob::Hit> hit : inputHits){
    unsigned int thePlane = hit.get()->WireID().asPlaneID().Plane;
    if( thePlane != planeID ) continue;

    showerHits.push_back(hit.get());

  }

  return showerHits;

}

//=================================================================================================================================================

std::vector<double> protoana::ProtoDUNEShowerUtils::GetdEdx5cm( const recob::Shower & shower, art::Event const & evt, double T0, calo::CalorimetryAlg caloAlg, const std::string showerModule, const std::string hitModule )
{
  // Get the geometry service
  art::ServiceHandle<geo::Geometry const> geo;

  // Get the shower's hits
  const std::vector< art::Ptr< recob::Hit > > hits = GetRecoShowerHitsArtPtr( shower, evt, showerModule);

  // Get the shower's start, direction, and length
  TVector3 showerStart = shower.ShowerStart();
  TVector3 showerDir   = shower.Direction().Unit();
  double   length      = shower.Length();
  double   pcent5      = length * 0.10;

  double kGeVtoElectrons { 4.237e7 };
  double recombination   { 1/0.63 };

  // Get the space points
  auto recoHits = evt.getValidHandle< std::vector< recob::Hit > >( hitModule );
  art::FindManyP< recob::SpacePoint > HitToSpoints( recoHits, evt, "pandora" );

  // Totals for the dE/dx and the number of hits
  //double totaldEdx = 0;
  std::vector<double> dEdxVector = { 0., 0., 0. };

  //std::cout << "There are " << hits.size() << std::endl;

  //std::map< double, art::Ptr< recob::Hit >> hitSpacePointMap;
  //int hitCounter    = 0;
  //int spointCounter = 0;
  for ( size_t h = 0 ; h < hits.size() ; h++ ) { // Loop over the hits in the track
    // Get the hit
    art::Ptr< recob::Hit > hit = hits[h];

    // Check to see if hit has a space point
    std::vector< art::Ptr<recob::SpacePoint> > sPoints = HitToSpoints.at( hit.key() );
    if ( sPoints.size() > 0 ) {
      // Get the associated space point
      art::Ptr< recob::SpacePoint > spoint = sPoints[0];

      // Vector of space point
      TVector3 XYZ = spoint->XYZ();

      // Vector from space point to start of shower
      TVector3 startToSpoint = XYZ - showerStart;

      // Projection of spoint along shower vector
      double projection = startToSpoint.Dot( showerDir );

      // If the space point's projection isn't in the range of interest, skip it.
      if ( projection > pcent5 ) {
        continue;
      }
    
      else {
        const int plane = hit->WireID().Plane;
        double hitEnergy = (( caloAlg.ElectronsFromADCArea( hit->Integral(), plane) * caloAlg.LifetimeCorrection(hit->PeakTime()) ) / kGeVtoElectrons) * recombination;
        // Convert to MeV
        hitEnergy *= 1000;

        // Find pitch JAT way
        float wire_pitch = geo->WirePitch( hit->View() );
        float this_pitch = wire_pitch;
        float angleToVert = geo->WireAngleToVertical( hit->View(), hit->WireID() );
        angleToVert -= 0.5*::util::pi<>();
        float cosgamma = std::abs( sin(angleToVert)*showerDir.Y() + cos(angleToVert)*showerDir.Z() );
        if ( cosgamma > 0 ) this_pitch /= cosgamma;
        if ( this_pitch < wire_pitch ) this_pitch = wire_pitch;

        double dEdx = hitEnergy/this_pitch;
        if ( dEdx < 5000.0 ) {
          dEdxVector[plane] += dEdx;
        }
      }
    }
  }
  return dEdxVector;
}

