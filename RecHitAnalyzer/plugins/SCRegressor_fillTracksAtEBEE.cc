#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Fill Tracks in EB+EE ////////////////////////////////
// Store tracks in EB+EE projection

// Initialize branches ____________________________________________________________//
void SCRegressor::branchesTracksAtEBEE ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  //tree->Branch("Tracks_EB",   &vTracks_EB_);
  tree->Branch("TracksPt_EB", &vTracksPt_EB_);
  //tree->Branch("TracksQPt_EB", &vTracksQPt_EB_);

  // Histograms for monitoring
  hTracks_EB = fs->make<TH2F>("Tracks_EB", "N(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );
  hTracksPt_EB = fs->make<TH2F>("TracksPt_EB", "pT(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );

  /*
  char hname[50], htitle[50];
  for ( int iz(0); iz < nEE; iz++ ) {
    // Branches for images
    const char *zside = (iz > 0) ? "p" : "m";
    sprintf(hname, "Tracks_EE%s",zside);
    tree->Branch(hname,        &vTracks_EE_[iz]);
    sprintf(hname, "TracksPt_EE%s",zside);
    tree->Branch(hname,        &vTracksPt_EE_[iz]);
    sprintf(hname, "TracksQPt_EE%s",zside);
    tree->Branch(hname,        &vTracksQPt_EE_[iz]);

    // Histograms for monitoring
    sprintf(hname, "Tracks_EE%s",zside);
    sprintf(htitle,"N(ix,iy);ix;iy");
    hTracks_EE[iz] = fs->make<TH2F>(hname, htitle,
        EE_MAX_IX, EE_MIN_IX-1, EE_MAX_IX,
        EE_MAX_IY, EE_MIN_IY-1, EE_MAX_IY );
    sprintf(hname, "TracksPt_EE%s",zside);
    sprintf(htitle,"pT(ix,iy);ix;iy");
    hTracksPt_EE[iz] = fs->make<TH2F>(hname, htitle,
        EE_MAX_IX, EE_MIN_IX-1, EE_MAX_IX,
        EE_MAX_IY, EE_MIN_IY-1, EE_MAX_IY );
  } // iz
  */

} // branchesEB()

// Fill TRK rechits at EB/EE ______________________________________________________________//
void SCRegressor::fillTracksAtEBEE ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  //int ix_, iy_, iz_;
  int iphi_, ieta_, idx_; // rows:ieta, cols:iphi
  float eta, phi;
  GlobalPoint pos;

  vTracks_EB_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  vTracksPt_EB_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  vTracksQPt_EB_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  /*
  for ( int iz(0); iz < nEE; iz++ ) {
    vTracks_EE_[iz].assign( EE_NC_PER_ZSIDE, 0. );
    vTracksPt_EE_[iz].assign( EE_NC_PER_ZSIDE, 0. );
    vTracksQPt_EE_[iz].assign( EE_NC_PER_ZSIDE, 0. );
  }
  */

  //edm::Handle<reco::TrackCollection> tracksH_;
  edm::Handle<pat::IsolatedTrackCollection> tracksH_;
  iEvent.getByToken( trackCollectionT_, tracksH_ );

  // Provides access to global cell position
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  //reco::Track::TrackQuality tkQt_ = reco::Track::qualityByName("highPurity");

  //for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
  for ( pat::IsolatedTrackCollection::const_iterator iTk = tracksH_->begin();
        iTk != tracksH_->end(); ++iTk ) {
    //if ( !(iTk->quality(tkQt_)) ) continue;

    eta = iTk->eta();
    phi = iTk->phi();
    //if ( std::abs(eta) > 3. ) continue;
    if ( std::abs(eta) > 1.5 ) continue;

    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalBarrel ) {
      EBDetId ebId( id );
      iphi_ = ebId.iphi() - 1;
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
      // Fill histograms for monitoring
      hTracks_EB->Fill( iphi_, ieta_ );
      hTracksPt_EB->Fill( iphi_, ieta_, iTk->pt() );
      idx_ = ebId.hashedIndex(); // (ieta_+EB_IETA_MAX)*EB_IPHI_MAX + iphi_
      // Fill vectors for images
      vTracks_EB_[idx_] += 1.;
      vTracksPt_EB_[idx_] += iTk->pt();
      vTracksQPt_EB_[idx_] += (iTk->charge()*iTk->pt());
    /*
    } else if ( id.subdetId() == EcalEndcap ) {
      EEDetId eeId( id );
      ix_ = eeId.ix() - 1;
      iy_ = eeId.iy() - 1;
      iz_ = (eeId.zside() > 0) ? 1 : 0;
      // Fill histograms for monitoring
      hTracks_EE[iz_]->Fill( ix_, iy_ );
      hTracksPt_EE[iz_]->Fill( ix_, iy_, iTk->pt() );
      // Create hashed Index: maps from [iy][ix] -> [idx_]
      idx_ = iy_*EE_MAX_IX + ix_;
      // Fill vectors for images
      vTracks_EE_[iz_][idx_] += 1.;
      vTracksPt_EE_[iz_][idx_] += iTk->pt();
      vTracksQPt_EE_[iz_][idx_] += (iTk->charge()*iTk->pt());
    */
    } 
  } // tracks

} // fillEB()
