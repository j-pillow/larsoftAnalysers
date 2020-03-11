#include "dune/Protodune/Analysis/ProtoDUNETruthUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEShowerUtils.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "art/Framework/Principal/Event.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TVector3.h"

protoana::ProtoDUNETruthUtils::ProtoDUNETruthUtils(){

}

protoana::ProtoDUNETruthUtils::~ProtoDUNETruthUtils(){

}

// Match Reco to True by Number of Hits contributed as opposed to Energy contributed
//
// List
template <typename T>
std::vector< std::pair< const simb::MCParticle *, size_t > > 
  protoana::ProtoDUNETruthUtils::GetMCParticleListByHits( const T &recobj, const art::Event &evt, std::string recoModule, std::string hitModule) const {
  
  //Get the particle list. This is sorted by energy contributed, but we'll turn this into nHits
  const std::vector< std::pair< const simb::MCParticle*, double > > mcparts = GetMCParticleListFromReco( recobj, evt, recoModule ); 

  using weightedMCPair = std::pair<const simb::MCParticle*, size_t>;
  std::vector< weightedMCPair > results;

  for( const std::pair< const simb::MCParticle*, double > part : mcparts ){
    const std::vector<const recob::Hit*> hits = GetSharedHits(*(part.first), recobj, evt, recoModule);
    results.push_back( std::make_pair( part.first, hits.size() ) );
  }

  std::sort(results.begin(), results.end(),
        [](weightedMCPair a, weightedMCPair b){ return a.second > b.second;});

  return results;
}

template std::vector< std::pair< const simb::MCParticle *, size_t > > protoana::ProtoDUNETruthUtils::GetMCParticleListByHits< recob::PFParticle >
  (const recob::PFParticle&, const art::Event&, std::string, std::string) const;
template std::vector< std::pair< const simb::MCParticle *, size_t > > protoana::ProtoDUNETruthUtils::GetMCParticleListByHits< recob::Track >
  (const recob::Track&, const art::Event&, std::string, std::string) const;
template std::vector< std::pair< const simb::MCParticle *, size_t > > protoana::ProtoDUNETruthUtils::GetMCParticleListByHits< recob::Shower >
  (const recob::Shower&, const art::Event&, std::string, std::string) const;
/////////////////

// Best Match
template <typename T>
const simb::MCParticle * protoana::ProtoDUNETruthUtils::GetMCParticleByHits ( const T &recobj, const art::Event &evt, std::string recoModule, std::string hitModule ) const {

  auto list = GetMCParticleListByHits( recobj, evt, recoModule, hitModule );

  if( !list.size() ) return 0x0;

  return list[0].first;
}

template const simb::MCParticle * protoana::ProtoDUNETruthUtils::GetMCParticleByHits< recob::PFParticle >
  (const recob::PFParticle&, const art::Event&, std::string, std::string) const;
template const simb::MCParticle * protoana::ProtoDUNETruthUtils::GetMCParticleByHits< recob::Track >
  (const recob::Track&, const art::Event&, std::string, std::string) const;
template const simb::MCParticle * protoana::ProtoDUNETruthUtils::GetMCParticleByHits< recob::Shower >
  (const recob::Shower&, const art::Event&, std::string, std::string) const;
/////////////////

// Function to find the best matched true particle to a reconstructed particle. In case of problems, returns a null pointer
std::vector< std::pair<const simb::MCParticle*, std::pair<int, double> > > protoana::ProtoDUNETruthUtils::GetMCParticleFromRecoTrackAll(const recob::Track &track, art::Event const & evt, std::string trackModule) const
{

  std::vector< std::pair<const simb::MCParticle*, std::pair<int, double> > > vecToReturn;

  // We must have MC for this module to make sense
  if(evt.isRealData()) {
    return vecToReturn;
  }

  // Get the reconstructed tracks
  auto allRecoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);

  // We need the association between the tracks and the hits
  const art::FindManyP<recob::Hit> findTrackHits(allRecoTracks, evt, trackModule);

  unsigned int trackIndex = track.ID();

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  std::unordered_map<int, double> trkIDE;
  for (auto const & h : findTrackHits.at(trackIndex))
  {
    for (auto const & ide : bt_serv->HitToTrackIDEs(h)) // loop over std::vector<sim::TrackIDE>
    {
        trkIDE[ide.trackID] += ide.energy; // sum energy contribution by each track ID
    }
  }

  for ( auto const & thisPair : trkIDE ) {

    auto thisID = thisPair.first; // What is the id of the first particle
    
    int negID = 0;
    if (thisID < 0) {            // Is the id negative?
      thisID = std::abs(thisID); // If so, abs it.
      negID = 1;                 // Record that it was negative
    }

    auto thisMcParticle = pi_serv->TrackIdToParticle_P(thisID); // Get the simb::MCParticle
    std::pair<int, double> negContrPair(negID, thisPair.second); // Pair up where it was negative, and the energy contribution
    std::pair<const simb::MCParticle*, std::pair<int, double>> particlePairPair(thisMcParticle, negContrPair); // Pair up the simb::MCParticle and the pair of the negative id and the energy contribution

    vecToReturn.push_back( particlePairPair ); // Push it to the vector
  }

  return vecToReturn;
}

// =========================================================================================================

// Function to find the best matched true particle to a reconstructed shower. In case of problems, returns a null pointer
std::vector< std::pair<const simb::MCParticle*, std::pair<int, double> > > protoana::ProtoDUNETruthUtils::GetMCParticleFromRecoShowerAll(const recob::Shower &shower, art::Event const & evt, std::string showerModule) const
{ 

  std::vector< std::pair<const simb::MCParticle*, std::pair<int, double> > > vecToReturn;
  
  if( evt.isRealData() ) { // We must have MC for this module to make sense
    return vecToReturn;
  }

  // Get the reconstructed tracks
  auto allRecoShowers = evt.getValidHandle<std::vector<recob::Shower> >(showerModule);

  // We need the association between the tracks and the hits
  const art::FindManyP<recob::Hit> findShowerHits(allRecoShowers, evt, showerModule);

  protoana::ProtoDUNEShowerUtils showerUtil;

  int showerIndex = showerUtil.GetShowerIndex( shower, evt, showerModule );

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  std::unordered_map<int, double> showerIDE;
  for ( auto const & h : findShowerHits.at(showerIndex) ) {
    
    for ( auto const & ide : bt_serv->HitToTrackIDEs(h) ) { // loop over std::vector<sim::TrackIDE>
        showerIDE[ide.trackID] += ide.energy; // sum energy contribution by each track ID
    }

  }

  for ( auto const & thisPair : showerIDE ) {

    auto thisID = thisPair.first; // What is the id of the first particle
    
    int negID = 0;
    if (thisID < 0) {            // Is the id negative?
      thisID = std::abs(thisID); // If so, abs it.
      negID = 1;                 // Record that it was negative
    }

    auto thisMcParticle = pi_serv->TrackIdToParticle_P(thisID); // Get the simb::MCParticle
    std::pair<int, double> negContrPair(negID, thisPair.second); // Pair up where it was negative, and the energy contribution
    std::pair<const simb::MCParticle*, std::pair<int, double>> particlePairPair(thisMcParticle, negContrPair); // Pair up the simb::MCParticle and the pair of the negative id and the energy contribution

    vecToReturn.push_back( particlePairPair ); // Push it to the vector
  }

  return vecToReturn;
}

// =========================================================================================================

// Function to find the best matched true particle to a reconstructed particle. In case of problems, returns a null pointer
std::vector< std::pair<const simb::MCParticle*, std::pair<int, double> > > protoana::ProtoDUNETruthUtils::GetMCParticleFromHitsAll(const std::vector<const recob::Hit*>& hitVec, art::Event const & evt ) const
{

  std::vector< std::pair<const simb::MCParticle*, std::pair<int, double> > > vecToReturn;

  // We must have MC for this module to make sense
  if(evt.isRealData()) {
    return vecToReturn;
  }


  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  std::unordered_map<int, double> trkIDE;
  for (auto const & h : hitVec)
  {
    for (auto const & ide : bt_serv->HitToTrackIDEs(*h)) // loop over std::vector<sim::TrackIDE>
    {
        trkIDE[ide.trackID] += ide.energy; // sum energy contribution by each track ID
    }
  }

  for ( auto const & thisPair : trkIDE ) {

    auto thisID = thisPair.first; // What is the id of the first particle
    
    int negID = 0;
    if (thisID < 0) {            // Is the id negative?
      thisID = std::abs(thisID); // If so, abs it.
      negID = 1;                 // Record that it was negative
    }

    auto thisMcParticle = pi_serv->TrackIdToParticle_P(thisID); // Get the simb::MCParticle
    std::pair<int, double> negContrPair(negID, thisPair.second); // Pair up where it was negative, and the energy contribution
    std::pair<const simb::MCParticle*, std::pair<int, double>> particlePairPair(thisMcParticle, negContrPair); // Pair up the simb::MCParticle and the pair of the negative id and the energy contribution

    vecToReturn.push_back( particlePairPair ); // Push it to the vector
  }

  return vecToReturn;
}

// =========================================================================================================

// Function that produces a vector of hit pointers that appear in both the given MCParticle and
// PFParticle.
std::vector<const recob::Hit*> protoana::ProtoDUNETruthUtils::GetSharedHits
  (const simb::MCParticle& mcpart, const recob::PFParticle& pfpart, const art::Event& evt, std::string pfparticleModule) const {

  std::vector<const recob::Hit*> outVec;

  // Get the hits associated with this PFParticle
  protoana::ProtoDUNEPFParticleUtils pfpUtils;
  std::vector<const recob::Hit*> pfpHits = pfpUtils.GetPFParticleHits(pfpart, evt, pfparticleModule);

  // Push back all hits that are shared with this MCParticle
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  for(const recob::Hit* hitp : pfpHits) {
    for(const int& trackId : bt_serv->HitToTrackIds(*hitp)) {
      if(trackId == mcpart.TrackId()) {
        outVec.push_back(hitp);
        break;
      }
    }
  }

  return outVec;
}

// Function that produces a vector of hit pointers that appear in both the given MCParticle and
// reconstructed track.
std::vector<const recob::Hit*> protoana::ProtoDUNETruthUtils::GetSharedHits
  (const simb::MCParticle& mcpart, const recob::Track& track, const art::Event& evt, std::string trackModule) const {

  std::vector<const recob::Hit*> outVec;

  // Get the hits associated with this track
  protoana::ProtoDUNETrackUtils trackUtils;
  std::vector<const recob::Hit*> trackHits = trackUtils.GetRecoTrackHits(track, evt, trackModule);

  // Push back all hits that are shared with this MCParticle
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  for(const recob::Hit* hitp : trackHits) {
    for(const int& trackId : bt_serv->HitToTrackIds(*hitp)) {
      if(trackId == mcpart.TrackId()) {
        outVec.push_back(hitp);
        break;
      }
    }
  }

  return outVec;
}

// Function that produces a vector of hit pointers that appear in both the given MCParticle and
// reconstructed shower.
std::vector<const recob::Hit*> protoana::ProtoDUNETruthUtils::GetSharedHits
  (const simb::MCParticle& mcpart, const recob::Shower& shower, const art::Event& evt, std::string showerModule) const {

  std::vector<const recob::Hit*> outVec;

  // Get the hits associated with this shower
  protoana::ProtoDUNEShowerUtils shUtils;
  std::vector<const recob::Hit*> shHits = shUtils.GetRecoShowerHits(shower, evt, showerModule);

  // Push back all hits that are shared with this MCParticle
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  for(const recob::Hit* hitp : shHits) {
    for(const sim::TrackIDE& ide : bt_serv->HitToEveTrackIDEs(*hitp)) {
      if(ide.trackID == mcpart.TrackId()) {
        outVec.push_back(hitp);
        break;
      }
    }
  }

  return outVec;
}

// Function that loops over all hits in an event and returns those that an MCParticle
// contributed to.
std::vector<const recob::Hit*> protoana::ProtoDUNETruthUtils::GetMCParticleHits(
  const simb::MCParticle &mcpart, const art::Event &evt, std::string hitModule) const {

  std::vector<const recob::Hit*> outVec;
  
  art::Handle<std::vector<recob::Hit>> hitHandle;
  if(!evt.getByLabel(hitModule, hitHandle)) {
    std::cout << "protoana::ProtoDUNETruthUtils::GetMCParticleHits: could not find hits in event.\n";
    return outVec;
  }

  // Backtrack all hits to verify whether they belong to the current MCParticle.
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  for(const recob::Hit& hit : *hitHandle) {
    for(const int trackId : bt_serv->HitToTrackIds(hit)) {
      if(pi_serv->TrackIdToParticle_P(trackId) == 
         pi_serv->TrackIdToParticle_P(mcpart.TrackId())) {
        outVec.push_back(&hit);
        break;
      }
    }
  }

  return outVec;
}

// Function to return the completeness of a reconstructed object. By definition:
// Completeness = number of shared hits between MCParticle and reco object /
//                number of hits in MCParticle
template <typename T>
double protoana::ProtoDUNETruthUtils::GetCompleteness(const T &recobj,
  const art::Event &evt, std::string recoModule, std::string hitModule) const {

  // Find the hits of the corresponding MCParticle.
  const simb::MCParticle* mcmatch = GetMCParticleFromReco(recobj, evt, recoModule);
  if(mcmatch == 0x0) { return 0; }
  // std::cout << "MCmatch: " << mcmatch << '\n';
  const std::vector<const recob::Hit*> MCHits = GetMCParticleHits(*mcmatch, evt, hitModule);
  if(MCHits.size() == 0) { return 0; }
  // Get shared hits between reco object and MCParticle.
  const std::vector<const recob::Hit*> sharedHits = GetSharedHits(*mcmatch, recobj, evt, recoModule);

  return (double)sharedHits.size() / MCHits.size();
}
template double protoana::ProtoDUNETruthUtils::GetCompleteness<recob::PFParticle>
  (const recob::PFParticle&, const art::Event&, std::string, std::string) const;
template double protoana::ProtoDUNETruthUtils::GetCompleteness<recob::Track>
  (const recob::Track&, const art::Event&, std::string, std::string) const;
template double protoana::ProtoDUNETruthUtils::GetCompleteness<recob::Shower>
  (const recob::Shower&, const art::Event&, std::string, std::string) const;

// Function to return the purity of a PFParticle. By definition:
// Purity = number of shared hits between MCParticle and reco object /
//          number of hits in reco object
double protoana::ProtoDUNETruthUtils::GetPurity(const recob::PFParticle &pfpart,
  const art::Event &evt, std::string pfparticleModule) const {
  
  // Find the corresponding MCParticle and shared hits.
  const simb::MCParticle* mcmatch = GetMCParticleFromReco(pfpart, evt, pfparticleModule);
  if(mcmatch == 0x0) { return 0; }
  // Get shared hits between PFParticle and MCParticle.
  const std::vector<const recob::Hit*> sharedHits = GetSharedHits(*mcmatch, pfpart, evt, pfparticleModule);

  // Get the hits associated with this PFParticle
  protoana::ProtoDUNEPFParticleUtils pfpUtils;
  std::vector<const recob::Hit*> pfpHits = pfpUtils.GetPFParticleHits(pfpart, evt, pfparticleModule);

  return (double)sharedHits.size() / pfpHits.size();
}

// Function to return the purity of a reco track. By definition:
// Purity = number of shared hits between MCParticle and reco object /
//          number of hits in reco object
double protoana::ProtoDUNETruthUtils::GetPurity(const recob::Track &track,
  const art::Event &evt, std::string trackModule) const {
  
  // Find the corresponding MCParticle and shared hits.
  const simb::MCParticle* mcmatch = GetMCParticleFromReco(track, evt, trackModule);
  if(mcmatch == 0x0) { return 0; }
  // Get shared hits between track and MCParticle.
  const std::vector<const recob::Hit*> sharedHits = GetSharedHits(*mcmatch, track, evt, trackModule);

  // Get the hits associated with this track
  protoana::ProtoDUNETrackUtils trackUtils;
  std::vector<const recob::Hit*> trackHits = trackUtils.GetRecoTrackHits(track, evt, trackModule);

  return (double)sharedHits.size() / trackHits.size();
}

// Function to return the purity of a reco shower. By definition:
// Purity = number of shared hits between MCParticle and reco object /
//          number of hits in reco object
double protoana::ProtoDUNETruthUtils::GetPurity(const recob::Shower &shower,
  const art::Event &evt, std::string showerModule) const {
  
  // Find the corresponding MCParticle and shared hits.
  const simb::MCParticle* mcmatch = GetMCParticleFromReco(shower, evt, showerModule);
  if(mcmatch == 0x0) { return 0; }
  // Get shared hits between shower and MCParticle.
  const std::vector<const recob::Hit*> sharedHits = GetSharedHits(*mcmatch, shower, evt, showerModule);

  // Get the hits associated with this shower
  protoana::ProtoDUNEShowerUtils showerUtils;
  std::vector<const recob::Hit*> showerHits = showerUtils.GetRecoShowerHits(shower, evt, showerModule);

  return (double)sharedHits.size() / showerHits.size();
}

// Function to find a weighted and sorted vector of the true particles matched to a
// vector of reconstructed track hits. In case of trouble, return an empty vector.
std::vector<std::pair<const simb::MCParticle*, double>>
  protoana::ProtoDUNETruthUtils::GetMCParticleListFromTrackHits
  (const std::vector<const recob::Hit*>& hitVec) const {

  using weightedMCPair = std::pair<const simb::MCParticle*, double>;
  std::vector<weightedMCPair> outVec;
  
  // Loop over all hits in the input vector and record the contributing MCParticles.
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  std::unordered_map<const simb::MCParticle*, double> mcEMap;
  double hitTotalE = 0;
  for(const recob::Hit* hitp : hitVec) {
    for(const sim::TrackIDE& ide : bt_serv->HitToTrackIDEs(*hitp)) {
      const simb::MCParticle* curr_part = pi_serv->TrackIdToParticle_P(ide.trackID);
      mcEMap[curr_part] += ide.energy;
      hitTotalE += ide.energy;
    }
  }

  // Fill and sort the output vector
  for (weightedMCPair const& p : mcEMap) {
    outVec.push_back(p);
  }
  std::sort(outVec.begin(), outVec.end(),
        [](weightedMCPair a, weightedMCPair b){ return a.second > b.second;});

  // Normalise the weights by the total track energy.
  if (hitTotalE < 1e-5) { hitTotalE = 1; } // Protect against zero division
  for (weightedMCPair& p : outVec) {
    p.second /= hitTotalE;
  }

  return outVec;
}

// Function to find a weighted and sorted vector of the true particles matched to a
// vector of reconstructed shower hits. In case of trouble, return an empty vector.
std::vector<std::pair<const simb::MCParticle*, double>>
  protoana::ProtoDUNETruthUtils::GetMCParticleListFromShowerHits
  (const std::vector<const recob::Hit*>& hitVec) const {

  using weightedMCPair = std::pair<const simb::MCParticle*, double>;
  std::vector<weightedMCPair> outVec;
  
  // Loop over all hits in the input vector and record the contributing MCParticles.
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  std::unordered_map<const simb::MCParticle*, double> mcEMap;
  double hitTotalE = 0;
  for(const recob::Hit* hitp : hitVec) {
    for(const sim::TrackIDE& ide : bt_serv->HitToEveTrackIDEs(*hitp)) {
      const simb::MCParticle* curr_part = pi_serv->TrackIdToParticle_P(ide.trackID);
      mcEMap[curr_part] += ide.energy;
      hitTotalE += ide.energy;
    }
  }

  // Fill and sort the output vector
  for (weightedMCPair const& p : mcEMap) {
    outVec.push_back(p);
  }
  std::sort(outVec.begin(), outVec.end(),
        [](weightedMCPair a, weightedMCPair b){ return a.second > b.second;});

  // Normalise the weights by the total track energy.
  if (hitTotalE < 1e-5) { hitTotalE = 1; } // Protect against zero division
  for (weightedMCPair& p : outVec) {
    p.second /= hitTotalE;
  }

  return outVec;
}

// Function to find a weighted and sorted vector of the true particles matched to
// a PFParticle. Return an empty vector in case of trouble.
std::vector<std::pair<const simb::MCParticle*, double>>
  protoana::ProtoDUNETruthUtils::GetMCParticleListFromPFParticle
  (const recob::PFParticle &pfpart, art::Event const & evt, std::string pfparticleModule) const {

  std::vector<std::pair<const simb::MCParticle*, double>> outVec;

  // We must have MC for this module to make sense
  if(evt.isRealData()) return outVec;

  // Get the hits associated with this PFParticle
  protoana::ProtoDUNEPFParticleUtils pfpUtils;
  std::vector<const recob::Hit*> pfpHits = pfpUtils.GetPFParticleHits(pfpart, evt, pfparticleModule);

  // Call the appropriate function to determine MC matches.
  if(abs(pfpart.PdgCode()) == 11 || pfpart.PdgCode() == 22) {
    outVec = GetMCParticleListFromShowerHits(pfpHits);
  } else {
    outVec = GetMCParticleListFromTrackHits(pfpHits);
  }

  return outVec;
}

std::vector<std::pair<const recob::PFParticle*, double>>
  protoana::ProtoDUNETruthUtils::GetPFParticleListFromMCParticle
  (const simb::MCParticle &part, art::Event const & evt, std::string pfparticleModule) const {

  using weightedPFPair = std::pair<const recob::PFParticle*, double>;

  std::vector<weightedPFPair> outVec;

  // We must have MC for this module to make sense
  if(evt.isRealData()) return outVec;

  // Get the reconstructed PFParticles
  auto allPFParticles = evt.getValidHandle<std::vector<recob::PFParticle> >(pfparticleModule);

  // Loop over all PFParticles and save the energetic overlap with each MCParticle in a map
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  std::unordered_map<const recob::PFParticle*, double> pfpEMap;
  double pfpTotalE = 0;
  for(const recob::PFParticle& pfpart : *allPFParticles) {
    for(const recob::Hit* hitp : GetSharedHits(part, pfpart, evt, pfparticleModule)) {
      for(const sim::TrackIDE& ide : bt_serv->HitToTrackIDEs(*hitp)) {
        if(ide.trackID == part.TrackId()) {
          pfpEMap[&pfpart] += ide.energy;
        }
        pfpTotalE += ide.energy;
      }
    }
  }

  // Fill and sort the output vector
  for (weightedPFPair const& p : pfpEMap) {
    outVec.push_back(p);
  }
  std::sort(outVec.begin(), outVec.end(),
        [](weightedPFPair a, weightedPFPair b){ return a.second > b.second;});

  // Normalise the weights by the total track energy.
  if (pfpTotalE < 1e-5) { pfpTotalE = 1; } // Protect against zero division
  for (weightedPFPair& p : outVec) {
    p.second /= pfpTotalE;
  }

  return outVec;
}

// Function to find a weighted and sorted vector of the true particles matched to
// a track. Return an empty vector in case of trouble.
std::vector<std::pair<const simb::MCParticle*, double>>
  protoana::ProtoDUNETruthUtils::GetMCParticleListFromRecoTrack
  (const recob::Track &track, art::Event const & evt, std::string trackModule) const {

  // We must have MC for this module to make sense
  if(evt.isRealData()) return {};

  // Get the reconstructed track hits
  protoana::ProtoDUNETrackUtils trackUtils;
  std::vector<const recob::Hit*> trackHits = trackUtils.GetRecoTrackHits(track, evt, trackModule);
  return GetMCParticleListFromTrackHits(trackHits);
}

// Function to find the best matched reconstructed particle tracks to a true
// particle. In case of problems, or if the particle was not reconstruced as a
// track, returns an empty vector.
std::vector<std::pair<const recob::Track*, double>>
  protoana::ProtoDUNETruthUtils::GetRecoTrackListFromMCParticle
  (const simb::MCParticle &part, art::Event const & evt, std::string trackModule) const{

  using weightedTrackPair = std::pair<const recob::Track*, double>;

  std::vector<weightedTrackPair> outVec;

  // We must have MC for this module to make sense
  if(evt.isRealData()) return outVec;

  // Get the reconstructed tracks
  auto allRecoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);

  // We need the association between the tracks and the hits
  const art::FindManyP<recob::Hit> findTrackHits(allRecoTracks, evt, trackModule);

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  std::unordered_map<int, double> recoTrack;

  // Record the energy contribution to the MCParticle of every relevant reco track
  double part_E = 0;
  for (recob::Track const & track : *allRecoTracks)
  {
    // This finds all hits that are shared between a reconstructed track and MCParticle
    for (art::Ptr<recob::Hit> const & hptr :
          bt_serv->TrackIdToHits_Ps(part.TrackId(), findTrackHits.at(track.ID())))
    {
      // Loop over hit IDEs to find our particle's part in it
      for (const sim::IDE* ide : bt_serv->HitToSimIDEs_Ps(hptr))
      {
        if (ide->trackID == part.TrackId())
        {
          recoTrack[track.ID()] += ide->energy;
        }
        part_E += ide->energy;
      } // sim::IDE*
    } // art::Ptr<recob::Hit>
  } // const recob::Track&

  // Fill and sort the output vector
  for (std::pair<int, double> const& p : recoTrack)
  {
    auto const trackIt = std::find_if(allRecoTracks->begin(), allRecoTracks->end(),
                               [&](recob::Track tr){ return tr.ID() == p.first; });
    outVec.push_back(std::make_pair(&*trackIt, p.second));
  }
  std::sort(outVec.begin(), outVec.end(),
    [](weightedTrackPair a, weightedTrackPair b){ return a.second > b.second;});

  // Normalise the vector weights
  if (part_E < 1e-5) { part_E = 1; } // Protect against zero division
  for (weightedTrackPair& p : outVec)
  {
    p.second /= part_E;
  }

  return outVec;
}

// Function to find the best matched true particle to a reconstructed particle
// shower. In case of problems, returns a null pointer
std::vector<std::pair<const simb::MCParticle*, double>>
  protoana::ProtoDUNETruthUtils::GetMCParticleListFromRecoShower
  (const recob::Shower &shower, art::Event const & evt, std::string showerModule) const{
    
  // We must have MC for this module to make sense
  if(evt.isRealData()) return {};

  // Get the reconstructed shower hits
  protoana::ProtoDUNEShowerUtils showerUtils;
  std::vector<const recob::Hit*> showerHits = showerUtils.GetRecoShowerHits(shower, evt, showerModule);
  return GetMCParticleListFromShowerHits(showerHits);
}

// Function to find the best matched reconstructed particle tracks to a true
// particle. In case of problems, or if the particle was not reconstruced as a
// shower, returns an empty vector.
std::vector<std::pair<const recob::Shower*, double>>
  protoana::ProtoDUNETruthUtils::GetRecoShowerListFromMCParticle
  (const simb::MCParticle &part, art::Event const & evt, std::string showerModule) const{

  using weightedShowerPair = std::pair<const recob::Shower*, double>;

  std::vector<weightedShowerPair> outVec;

  // We must have MC for this module to make sense
  if(evt.isRealData()) return outVec;

  // Get the reconstructed showers
  auto allRecoShowers = evt.getValidHandle<std::vector<recob::Shower> >(showerModule);

  // We need the association between the showers and the hits
  const art::FindManyP<recob::Hit> findShowerHits(allRecoShowers, evt, showerModule);

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  ProtoDUNEShowerUtils shUtils;
  std::unordered_map<int, double> recoShower;

  // Record the energy contribution to the MCParticle of every relevant reco shower
  double part_E = 0;
  for (recob::Shower const & shower : *allRecoShowers)
  {
    const unsigned showerIndex = shUtils.GetShowerIndex(shower, evt, showerModule);
    // Since MCParticles do not have shower hits associated with them, loop over all shower hits
    for (art::Ptr<recob::Hit> hit : findShowerHits.at(showerIndex))
    {
      // Loop over hit IDEs to find our particle's part in it
      for (const sim::TrackIDE& ide : bt_serv->HitToEveTrackIDEs(hit))
      {
        if (ide.trackID == part.TrackId())
        {
          recoShower[showerIndex] += ide.energy;
        }
        part_E += ide.energy;
      } // sim::IDE*
    } // art::Ptr<recob::Hit>
  } // const recob::Shower&

  // Fill and sort the output vector
  for (std::pair<int, double> const& p : recoShower)
  {
    auto const showerIt = std::find_if(allRecoShowers->begin(), allRecoShowers->end(),
      [&](recob::Shower sh){ return shUtils.GetShowerIndex(sh, evt, showerModule) == p.first; });
    outVec.push_back(std::make_pair(&*showerIt, p.second));
  }
  std::sort(outVec.begin(), outVec.end(),
    [](weightedShowerPair a, weightedShowerPair b){ return a.second > b.second;});

  // Normalise the vector weights
  if (part_E < 1e-5) { part_E = 1; } // Protect against zero division
  std::for_each(outVec.begin(), outVec.end(),
          [&](weightedShowerPair& p){ p.second /= part_E; });

  return outVec;
}

// Function to find the best matched true particle to a reconstructed PFParticle.
// In case of problems, returns a null pointer
const simb::MCParticle* protoana::ProtoDUNETruthUtils::GetMCParticleFromPFParticle
  (const recob::PFParticle &pfpart, art::Event const & evt, std::string pfparticleModule) const{

  const simb::MCParticle* outPart = 0x0;

  const std::vector<std::pair<const simb::MCParticle*, double>> inVec =
    GetMCParticleListFromPFParticle(pfpart, evt, pfparticleModule);

  if (inVec.size() != 0) outPart = inVec[0].first;

  return outPart;
}

// Function to find the best matched reconstructed PFParticle to a true particle. In
// case of problems, or if the true particle is not a primary contributor to any
// PFParticle, return a null pointer
const recob::PFParticle* protoana::ProtoDUNETruthUtils::GetPFParticleFromMCParticle
  (const simb::MCParticle &part, art::Event const & evt, std::string pfparticleModule) const {

  const recob::PFParticle* outPFPart = 0x0;

  // We must have MC for this module to make sense
  if(evt.isRealData()) return outPFPart;

  // Get the list of contributing PFParticles for the MCParticle and see if any of
  // them have the MCParticle as a primary contributor.
  for (std::pair<const recob::PFParticle*, double> p :
                        GetPFParticleListFromMCParticle(part, evt, pfparticleModule))
  {
    const recob::PFParticle* pfp = p.first;
    if (GetMCParticleFromPFParticle(*pfp, evt, pfparticleModule)->TrackId() == part.TrackId())
    {
      outPFPart = pfp;
      break;
    }
  }

  return outPFPart;
}

// Function to find the best matched true particle to a reconstructed particle
// track. In case of problems, returns a null pointer
const simb::MCParticle* protoana::ProtoDUNETruthUtils::GetMCParticleFromRecoTrack
  (const recob::Track &track, art::Event const & evt, std::string trackModule) const{

  const simb::MCParticle* outPart = 0x0;

  const std::vector<std::pair<const simb::MCParticle*, double>> inVec =
    GetMCParticleListFromRecoTrack(track, evt, trackModule);

  if (inVec.size() != 0) outPart = inVec[0].first;

  return outPart;
}

// Function to find the best matched reconstructed track to a true particle. In
// case of problems, or if the true particle is not a primary contributor to any
// track, return a null pointer
const recob::Track* protoana::ProtoDUNETruthUtils::GetRecoTrackFromMCParticle
  (const simb::MCParticle &part, art::Event const & evt, std::string trackModule) const {

  const recob::Track* outTrack = 0x0;

  // We must have MC for this module to make sense
  if(evt.isRealData()) return outTrack;

  // Get the list of contributing tracks for the MCParticle and see if any of
  // them have the MCParticle as a primary contributor.
  for (std::pair<const recob::Track*, double> p :
                        GetRecoTrackListFromMCParticle(part, evt, trackModule))
  {
    const recob::Track* tr = p.first;
    if (GetMCParticleFromRecoTrack(*tr, evt, trackModule)->TrackId() == part.TrackId())
    {
      outTrack = tr;
      break;
    }
  }

  return outTrack;
}

// Function to find the best matched true particle to a reconstructed particle
// shower. In case of problems, returns a null pointer
const simb::MCParticle* protoana::ProtoDUNETruthUtils::GetMCParticleFromRecoShower
  (const recob::Shower &shower, art::Event const & evt, std::string showerModule) const{

  const simb::MCParticle* outPart = 0x0;

  const std::vector<std::pair<const simb::MCParticle*, double>> inVec =
    GetMCParticleListFromRecoShower(shower, evt, showerModule);

  if (inVec.size() != 0) outPart = inVec[0].first;

  return outPart;
}

// Function to find the best matched reconstructed shower to a true particle. In
// case of problems, or if the true particle is not a primary contributor to any
// shower, returns a null pointer
const recob::Shower* protoana::ProtoDUNETruthUtils::GetRecoShowerFromMCParticle
  (const simb::MCParticle &part, art::Event const & evt, std::string showerModule) const {

  const recob::Shower* outShower = 0x0;

  // We must have MC for this module to make sense
  if(evt.isRealData()) return outShower;

  // Get the list of contributing showers for the MCParticle and see if any of
  // them have the MCParticle as a primary contributor.
  for (std::pair<const recob::Shower*, double> p :
                        GetRecoShowerListFromMCParticle(part, evt, showerModule))
  {
    const recob::Shower* tr = p.first;
    if (GetMCParticleFromRecoShower(*tr, evt, showerModule)->TrackId() == part.TrackId())
    {
      outShower = tr;
      break;
    }
  }

  return outShower;
}

// General overloaded functions to get the best MCParticle match for reconstructed objects
const simb::MCParticle* protoana::ProtoDUNETruthUtils::GetMCParticleFromReco
  (const recob::PFParticle &pfpart, art::Event const &evt, std::string pfparticleModule) const {
  return GetMCParticleFromPFParticle(pfpart, evt, pfparticleModule);
}
const simb::MCParticle* protoana::ProtoDUNETruthUtils::GetMCParticleFromReco
  (const recob::Track &track, art::Event const &evt, std::string trackModule) const {
  return GetMCParticleFromRecoTrack(track, evt, trackModule);
}
const simb::MCParticle* protoana::ProtoDUNETruthUtils::GetMCParticleFromReco
  (const recob::Shower &shower, art::Event const &evt, std::string showerModule) const {
  return GetMCParticleFromRecoShower(shower, evt, showerModule);
}

// General overloaded functions to get the contributing MCParticles to a reconstructed object
std::vector<std::pair<const simb::MCParticle*, double>>
  protoana::ProtoDUNETruthUtils::GetMCParticleListFromReco
  (const recob::PFParticle &pfpart, art::Event const &evt, std::string pfparticleModule) const {
  return GetMCParticleListFromPFParticle(pfpart, evt, pfparticleModule);
}
std::vector<std::pair<const simb::MCParticle*, double>>
  protoana::ProtoDUNETruthUtils::GetMCParticleListFromReco
  (const recob::Track &track, art::Event const &evt, std::string trackModule) const {
  return GetMCParticleListFromRecoTrack(track, evt, trackModule);
}
std::vector<std::pair<const simb::MCParticle*, double>>
  protoana::ProtoDUNETruthUtils::GetMCParticleListFromReco
  (const recob::Shower &shower, art::Event const &evt, std::string showerModule) const {
  return GetMCParticleListFromRecoShower(shower, evt, showerModule);
}


const simb::MCParticle* protoana::ProtoDUNETruthUtils::MatchPduneMCtoG4( const simb::MCParticle & pDunePart, const art::Event & evt )
{  // Function that will match the protoDUNE MC particle to the Geant 4 MC particle, and return the matched particle (or a null pointer).

  // Find the energy of the procided MC Particle
  double pDuneEnergy = pDunePart.E();

  // Get list of the g4 particles. plist should be a std::map< int, simb::MCParticle* >
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  const sim::ParticleList & plist = pi_serv->ParticleList();

  // Check if plist is empty
  if ( !plist.size() ) {
    std::cerr << "\n\n#####################################\n"
              << "\nEvent " << evt.id().event() << "\n"
              << "sim::ParticleList from cheat::ParticleInventoryService is empty\n"
              << "A null pointer will be returned\n"
              << "#####################################\n\n";
    return nullptr;
  }

  // Go through the list of G4 particles
  for ( auto partIt = plist.begin() ; partIt != plist.end() ; partIt++ ) {

    const simb::MCParticle* pPart = partIt->second;
    if (!pPart) {
      std::cerr << "\n\n#####################################\n"
                << "\nEvent " << evt.id().event() << "\n"
                << "GEANT particle #" << partIt->first << " returned a null pointer\n"
                << "This is not necessarily bad. It just means at least one\n"
                << "of the G4 particles returned a null pointer. It may well\n"
                << "have still matched a PD particle and a G4 particle.\n"
                << "#####################################\n\n";
      continue;
    }

    // If the initial energy of the g4 particle is very close to the energy of the protoDUNE particle, call it a day and have a cuppa.
    if ( (pDunePart.PdgCode() == pPart->PdgCode()) && fabs(pPart->E() - pDuneEnergy) < 0.00001 ) {
      return pPart;
    }

  }  // G4 particle list loop end.

  std::cout << "No G4 particle was matched for Event " << evt.id().event() << ". Null pointer returned\n";
  return nullptr;

}  // End MatchPduneMCtoG4

const simb::MCParticle* protoana::ProtoDUNETruthUtils::GetGeantGoodParticle(const simb::MCTruth &genTruth, const art::Event &evt) const{

  // Get the good particle from the MCTruth
  simb::MCParticle goodPart;
  bool found = false;
  for(int t = 0; t < genTruth.NParticles(); ++t){
    simb::MCParticle part = genTruth.GetParticle(t);
    if(part.Process() == "primary"){
      goodPart = part;
      found = true;
      break;
    }
  }

  if(!found){
    std::cerr << "No good particle found, returning null pointer" << std::endl;
    return nullptr;
  }

  // Now loop over geant particles to find the one that matches
  // Get list of the g4 particles. plist should be a std::map< int, simb::MCParticle* >
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  const sim::ParticleList & plist = pi_serv->ParticleList();

  for(auto const part : plist){
    if((goodPart.PdgCode() == part.second->PdgCode()) && fabs(part.second->E() - goodPart.E()) < 1e-5){
      return part.second;
    }
  }

  // If we get here something has gone wrong
  std::cerr << "No G4 version of the good particle was found, returning null pointer" << std::endl;
  return nullptr;

}

// Converting times in LArSoft can be a bit of a minefield. These functions convert true times in ns
// to pandora times in ns
const float protoana::ProtoDUNETruthUtils::ConvertTrueTimeToPandoraTimeNano(const simb::MCParticle &part) const{
  return ConvertTrueTimeToPandoraTimeNano(part.T());
}

const float protoana::ProtoDUNETruthUtils::ConvertTrueTimeToPandoraTimeNano(const float trueTime) const{
  return 1000. * ConvertTrueTimeToPandoraTimeMicro(trueTime);
}

// Microsecond versions
const float protoana::ProtoDUNETruthUtils::ConvertTrueTimeToPandoraTimeMicro(const simb::MCParticle &part) const{
  return ConvertTrueTimeToPandoraTimeMicro(part.T());
}

const float protoana::ProtoDUNETruthUtils::ConvertTrueTimeToPandoraTimeMicro(const float trueTime) const{

  // Use the clocks service to account for the offset between the Geant4 time and the electronics clock
  auto const* detclock = lar::providerFrom<detinfo::DetectorClocksService>();

  return detclock->G4ToElecTime(trueTime);
}

// Get process key.
int protoana::ProtoDUNETruthUtils::GetProcessKey(std::string process){

  if(process.compare("primary") == 0)                    return 0;
  else if(process.compare("hadElastic") == 0)            return 1;
  else if(process.compare("pi-Inelastic") == 0)          return 2;
  else if(process.compare("pi+Inelastic") == 0)          return 3;
  else if(process.compare("kaon-Inelastic") == 0)        return 4;
  else if(process.compare("kaon+Inelastic") == 0)        return 5;
  else if(process.compare("protonInelastic") == 0)       return 6;
  else if(process.compare("neutronInelastic") == 0)      return 7;
  else if(process.compare("kaon0SInelastic") == 0)       return 8;
  else if(process.compare("kaon0LInelastic") == 0)       return 9;
  else if(process.compare("lambdaInelastic") == 0)       return 10;
  else if(process.compare("omega-Inelastic") == 0)       return 11;
  else if(process.compare("sigma+Inelastic") == 0)       return 12;
  else if(process.compare("sigma-Inelastic") == 0)       return 13;
  else if(process.compare("sigma0Inelastic") == 0)       return 14;
  else if(process.compare("xi-Inelastic") == 0)          return 15;
  else if(process.compare("xi0Inelastic") == 0)          return 16;
  else if(process.compare("anti_protonInelastic") == 0)  return 20;
  else if(process.compare("anti_neutronInelastic") == 0) return 21;
  else if(process.compare("anti_lambdaInelastic") == 0)  return 22;
  else if(process.compare("anti_omega-Inelastic") == 0)  return 23;
  else if(process.compare("anti_sigma+Inelastic") == 0)  return 24;
  else if(process.compare("anti_sigma-Inelastic") == 0)  return 25;
  else if(process.compare("anti_xi-Inelastic") == 0)     return 26;
  else if(process.compare("anti_xi0Inelastic") == 0)     return 27;

  else if(process.compare("Decay") == 0)                 return 30;
  else if(process.compare("FastScintillation") == 0)     return 31;
  else if(process.compare("nKiller") == 0)               return 32; // Remove unwanted neutrons: neutron kinetic energy threshold (default 0) or time limit for neutron track
  else if(process.compare("nCapture") == 0)              return 33; // Neutron capture

  else if(process.compare("compt") == 0)                 return 40; // Compton Scattering
  else if(process.compare("rayleigh") == 0)              return 41; // Rayleigh Scattering
  else if(process.compare("phot") == 0)                  return 42; // Photoelectric Effect
  else if(process.compare("conv") == 0)                  return 43; // Pair production
  else if(process.compare("CoupledTransportation") == 0) return 44; //

  else return -1;
}

// Get estimated particle energy deposit. The G4 trackId must be provided
double protoana::ProtoDUNETruthUtils::GetDepEnergyMC(const art::Event &evt, geo::GeometryCore const * fGeom, int trackid, int whichview) const {

  double edep = 0.0;

  art::Handle< std::vector<sim::SimChannel> > simchannelHandle;
  if(evt.getByLabel("largeant", simchannelHandle)){
    // Loop over sim channels
    for(auto const& simchannel : (*simchannelHandle)){
      // Only keep channels in the selected view
      if(fGeom->View(simchannel.Channel()) != whichview) continue;
      // Get all time slices
      auto const& alltimeslices = simchannel.TDCIDEMap();
      // Loop over time slices
      for(auto const& tslice : alltimeslices){
	auto const& simide = tslice.second;
	// Loop over energy deposits
	for(auto const& eDep : simide){
	  if(eDep.trackID == trackid || eDep.trackID == -trackid)
	    edep += eDep.energy;
	}
      }
    }
  }

  return edep;

}

// Get first trajectory point in TPC active volume
int protoana::ProtoDUNETruthUtils::GetFirstTrajectoryPointInTPCActiveVolume(const simb::MCParticle& mcpart, double tpcactiveXlow, double tpcactiveXhigh, double tpcactiveYlow, double tpcactiveYhigh, double tpcactiveZlow, double tpcactiveZhigh){

  int firstpoint = -999;
  for(unsigned int i = 0; i < mcpart.NumberTrajectoryPoints(); ++i) {
    if(mcpart.Vx(i) >= tpcactiveXlow && mcpart.Vx(i) <= tpcactiveXhigh && mcpart.Vy(i) >= tpcactiveYlow && mcpart.Vy(i) <= tpcactiveYhigh && mcpart.Vz(i) >= tpcactiveZlow && mcpart.Vz(i) <= tpcactiveZhigh){

      firstpoint = i;
      break;
    }
  }

  return firstpoint;
}

// Get MC Particle length in TPC active volume
double protoana::ProtoDUNETruthUtils::GetMCParticleLengthInTPCActiveVolume(const simb::MCParticle& mcpart, double tpcactiveXlow, double tpcactiveXhigh, double tpcactiveYlow, double tpcactiveYhigh, double tpcactiveZlow, double tpcactiveZhigh){

  double length = 0.0;
  int firstpoint = GetFirstTrajectoryPointInTPCActiveVolume(mcpart, tpcactiveXlow, tpcactiveXhigh, tpcactiveYlow, tpcactiveYhigh, tpcactiveZlow, tpcactiveZhigh);

  if(firstpoint < 0) return length;

  TVector3 pos =  mcpart.Position(firstpoint).Vect();
  for(unsigned int i = firstpoint+1; i < mcpart.NumberTrajectoryPoints(); ++i) {
    if(mcpart.Vx(i) >= tpcactiveXlow && mcpart.Vx(i) <= tpcactiveXhigh && mcpart.Vy(i) >= tpcactiveYlow && mcpart.Vy(i) <= tpcactiveYhigh && mcpart.Vz(i) >= tpcactiveZlow && mcpart.Vz(i) <= tpcactiveZhigh){

      pos -= mcpart.Position(i).Vect();
      length += pos.Mag();
      pos = mcpart.Position(i).Vect();
    }
  }

  return length;
}

// Estimate last point energy loss
double protoana::ProtoDUNETruthUtils::GetDepEnergyAtLastTrajPoint(const simb::MCParticle& mcpart, double tpcactiveXlow, double tpcactiveXhigh, double tpcactiveYlow, double tpcactiveYhigh, double tpcactiveZlow, double tpcactiveZhigh){

  // First trajectory point inside TPC active volume
  int firstpoint = GetFirstTrajectoryPointInTPCActiveVolume(mcpart, tpcactiveXlow, tpcactiveXhigh, tpcactiveYlow, tpcactiveYhigh, tpcactiveZlow, tpcactiveZhigh);
  if(firstpoint < 0) return 0.0;

  // True track length in TPC active volume
  double length = GetMCParticleLengthInTPCActiveVolume(mcpart, tpcactiveXlow, tpcactiveXhigh, tpcactiveYlow, tpcactiveYhigh, tpcactiveZlow, tpcactiveZhigh);

  // Get the true track length up to the next to the last trajectory point
  int ntrajpoints = mcpart.NumberTrajectoryPoints();
  double plength = 0.0;
  TVector3 pos =  mcpart.Position(firstpoint).Vect();
  for(int i = firstpoint+1; i < ntrajpoints-1; ++i) {
    if(mcpart.Vx(i) >= tpcactiveXlow && mcpart.Vx(i) <= tpcactiveXhigh && mcpart.Vy(i) >= tpcactiveYlow && mcpart.Vy(i) <= tpcactiveYhigh && mcpart.Vz(i) >= tpcactiveZlow && mcpart.Vz(i) <= tpcactiveZhigh){

      pos -= mcpart.Position(i).Vect();
      plength += pos.Mag();
      pos = mcpart.Position(i).Vect();
    }
  }
  if(plength == 0.0) return 0.0;

  double length_laststep = length - plength;
  if(length_laststep <= 0.0) return 0.0;
  
  // Kinetic energy of the track
  double kinene_start = std::sqrt(mcpart.P(firstpoint)*mcpart.P(firstpoint) + mcpart.Mass()*mcpart.Mass()) - mcpart.Mass();
  double kinene_end = std::sqrt(mcpart.P(ntrajpoints-2)*mcpart.P(ntrajpoints-2) + mcpart.Mass()*mcpart.Mass()) - mcpart.Mass();
  double eneloss = (kinene_start - kinene_end);
  if(eneloss < 0.0 || kinene_start <= 0.0) return 0.0;

  // Average energy loss due to last trajectory point
  double aveneloss_laststep = length_laststep*eneloss/plength;

  return aveneloss_laststep;
  
}

// Get particle's kinetic energy at the interaction vertex
double protoana::ProtoDUNETruthUtils::GetKinEnergyAtVertex(const simb::MCParticle& mcpart, double kinene_lastpoint){

  int ntrajpoints = mcpart.NumberTrajectoryPoints();
  double kinene_end = std::sqrt(mcpart.P(ntrajpoints-2)*mcpart.P(ntrajpoints-2) + mcpart.Mass()*mcpart.Mass()) - mcpart.Mass();
  
  return (kinene_end - kinene_lastpoint);
}

// Get the sim::IDEs from the MCParticle, organized by the trajectory points
std::map< size_t, std::vector< const sim::IDE * > > protoana::ProtoDUNETruthUtils::GetSimIDEs( const simb::MCParticle & mcpart ){
  art::ServiceHandle< cheat::BackTrackerService >       bt_serv;
  //art::ServiceHandle< geo::Geometry >                   geom;

  const simb::MCTrajectory & mctraj = mcpart.Trajectory();

  //Get all ides from the MCParticle
  std::vector< const sim::IDE * > ides = bt_serv->TrackIdToSimIDEs_Ps( mcpart.TrackId()/*, geom->View(2)*/ );
  //Sort by z position, assume traveling in beam direction.
  std::sort( ides.begin(), ides.end(), sort_IDEs );

  std::map< size_t, std::vector< const sim::IDE * > > results;
  
  size_t ide_start = 0;

  for( size_t i = 0; i < mctraj.size(); ++i ){
    results[i] = std::vector< const sim::IDE * >();  

    for( size_t j = ide_start; j < ides.size(); ++j ){

      if( ( ides[j]->z > mctraj.Z(i) ) && ( ides[j]->z < mctraj.Z(i+1) ) ){
        results[i].push_back( ides[j] );
      }

      else if( ( ides[j]->z < mctraj.Z(i) ) && ( ides[j]->z < mctraj.Z(i+1) ) ){
        continue; 
      }

      else{
        ide_start = j;
        break;
      }
    }
  }

  return results;
}
