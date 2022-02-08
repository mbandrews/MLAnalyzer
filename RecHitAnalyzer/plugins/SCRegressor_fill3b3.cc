#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Initialize branches _____________________________________________________//
void SCRegressor::branches3b3 ( TTree* tree, edm::Service<TFileService> &fs )
{

  tree->Branch("pho_mass3b3", &vPho_mass3b3_);
  tree->Branch("pho_dR3b3",   &vPho_dR3b3_);
  tree->Branch("pho_ieta03b3",   &vPho_ieta03b3_);
  tree->Branch("pho_ieta13b3",   &vPho_ieta13b3_);
  tree->Branch("pho_iphi03b3",   &vPho_iphi03b3_);
  tree->Branch("pho_iphi13b3",   &vPho_iphi13b3_);

  hPho_mass3b3 = fs->make<TH1F>("hPho_mass3b3", "pho_mass3b3;mass;N", 48, 0., 1.2);

} // branches3b3()

// Fill 3b3 rechits _________________________________________________________________//
void SCRegressor::fill3b3 ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  // Adapted from DQMOffline/EGamma/plugins/PiZeroAnalyzer.cc

  using namespace std;
  using namespace edm;

  // Set default mass values to each regressable selected photon object
  // To be overwritten only if pi0 cand found
  float default_mass = -9.;
  std::vector<float> vPho_mass3b3;
  vPho_mass3b3_.assign( vRegressPhoIdxs_.size(), default_mass );
  vPho_dR3b3_.assign( vRegressPhoIdxs_.size(), default_mass );
  //std::cout << ">> N regressable photons:" << vRegressPhoIdxs_.size() << std::endl;
  //std::cout << ">> vPho_mass3b3_.size():" << vPho_mass3b3_.size() << std::endl;
  vPho_ieta03b3_.assign( vRegressPhoIdxs_.size(), default_mass );
  vPho_ieta13b3_.assign( vRegressPhoIdxs_.size(), default_mass );
  vPho_iphi03b3_.assign( vRegressPhoIdxs_.size(), default_mass );
  vPho_iphi13b3_.assign( vRegressPhoIdxs_.size(), default_mass );
  int iphi_, ieta_; // rows:ieta, cols:iphi

  //////////////// Get handles ///////////////////////////

  // Rechits
  edm::Handle<EcalRecHitCollection> rhEB;
  iEvent.getByToken(EBRecHitCollectionT_, rhEB);

  edm::Handle<EcalRecHitCollection> rhEE;
  iEvent.getByToken(EBRecHitCollectionT_, rhEE);

  edm::Handle<PhotonCollection> photonH_;
  iEvent.getByToken(photonCollectionT_, photonH_);

  const EcalRecHitCollection *hitCollection_p = rhEB.product();

  // Geometry
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);

  edm::ESHandle<CaloTopology> theCaloTopology;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopology);

  const CaloSubdetectorGeometry *geometry_p;
  const CaloSubdetectorTopology *topology_p;
  const CaloSubdetectorGeometry *geometryES_p;
  geometry_p = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  geometryES_p = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
  const CaloGeometry* caloGeom = geoHandle.product();

  //////////////// pi0 reconstruction ///////////////////////////

  // Parameters for the position calculation:
  PositionCalc posCalculator_ = PositionCalc(posCalcParameters_);

  std::map<DetId, EcalRecHit> recHitsEB_map;

  std::vector<EcalRecHit> seeds;
  seeds.clear();

  vector<EBDetId> usedXtals;
  usedXtals.clear();

  EcalRecHitCollection::const_iterator itb;

  static const int MAXCLUS = 2000;
  int nClus=0;
  vector<float> eClus;
  vector<float> etClus;
  vector<float> etaClus;
  vector<float> phiClus;
  vector<EBDetId> max_hit;
  vector< vector<EcalRecHit> > RecHitsCluster;
  vector<float> s4s9Clus;

  // __ Find cluster seeds around selected photons ____________________________________//

  bool passedDR;
  float dR;
  float eta, phi;
  GlobalPoint pos;

  // Loop over rechit collection
  for(itb=rhEB->begin(); itb!=rhEB->end(); ++itb){
    EBDetId id(itb->id());
    double energy = itb->energy();
    if (energy > seleXtalMinEnergy_) {
      std::pair<DetId, EcalRecHit> map_entry(itb->id(), *itb);
      recHitsEB_map.insert(map_entry);
    }
    //if (energy > clusSeedThr_) seeds.push_back(*itb);
    if (energy < clusSeedThr_) continue;

    // Get position of seed crystal
    pos = caloGeom->getPosition( id );
    eta = pos.eta();
    phi = pos.phi();
    passedDR = false;
    for ( unsigned int iP = 0; iP < vRegressPhoIdxs_.size(); iP++ ) {
      dR = reco::deltaR( eta,phi, vPho_eta_[iP],vPho_phi_[iP] );
      //if ( dR < 2.*0.04 ) {
      if ( dR < 16*0.0174 ) {
        passedDR = true;
        break;
      }
    } // diphotons
    if ( !passedDR ) continue;
    seeds.push_back(*itb);
    //std::cout << ">> seed (eta,phi,E):" << eta << ","<<phi<< ","<<energy << std::endl;
    //std::cout << ">> seed (ieta,iphi,E):" << id.ieta() << ","<<id.iphi()<< ","<<energy << std::endl;
  } // Eb rechits
  //std::cout << ">> N seeds:" << seeds.size() << std::endl;

  //___ Get rechits associated with seeds to form clusters ____________________________________//

  // Loop over seeds
  // Sort seeds by energy => clusters arranged by energy
  //sort(seeds.begin(), seeds.end(), ecalRecHitLess());
  sort(seeds.begin(), seeds.end(), [](auto& x, auto& y){return (x.energy() > y.energy());});
  for (std::vector<EcalRecHit>::iterator itseed=seeds.begin(); itseed!=seeds.end(); itseed++) {
    EBDetId seed_id = itseed->id();
    std::vector<EBDetId>::const_iterator usedIds;

    // Remove duplicate seeds
    bool seedAlreadyUsed=false;
    for(usedIds=usedXtals.begin(); usedIds!=usedXtals.end(); usedIds++){
      if(*usedIds==seed_id){
        seedAlreadyUsed=true;
        //cout<< " Seed with energy "<<itseed->energy()<<" was used !"<<endl;
        break;
      }
    }
    if(seedAlreadyUsed)continue;
    topology_p = theCaloTopology->getSubdetectorTopology(DetId::Ecal,EcalBarrel);
    std::vector<DetId> clus_v = topology_p->getWindow(seed_id,clusEtaSize_,clusPhiSize_);
    //std::vector<DetId> clus_used;
    std::vector<std::pair<DetId, float> > clus_used;

    vector<EcalRecHit> RecHitsInWindow;

    double simple_energy = 0;

    // Get rechits associated with seed that make the cluster
    //std::cout << "clus_v.size:"<<clus_v.size()<<std::endl;
    for (std::vector<DetId>::iterator det=clus_v.begin(); det!=clus_v.end(); det++) {
      // EBDetId EBdet = *det;
      //      cout<<" det "<< EBdet<<" ieta "<<EBdet.ieta()<<" iphi "<<EBdet.iphi()<<endl;
      bool  HitAlreadyUsed=false;
      for(usedIds=usedXtals.begin(); usedIds!=usedXtals.end(); usedIds++){
        if(*usedIds==*det){
          HitAlreadyUsed=true;
          break;
        }
      }
      if(HitAlreadyUsed)continue;
      if (recHitsEB_map.find(*det) != recHitsEB_map.end()){
        //      cout<<" Used det "<< EBdet<<endl;
        std::map<DetId, EcalRecHit>::iterator aHit;
        aHit = recHitsEB_map.find(*det);
        usedXtals.push_back(*det);
        RecHitsInWindow.push_back(aHit->second);
        clus_used.push_back( std::pair<DetId, float>(*det, 1.) );
        simple_energy = simple_energy + aHit->second.energy();
        //std::cout << " simple E:"<<simple_energy<<std::endl;
      }
    }

    // Store cluster kinematics
    math::XYZPoint clus_pos = posCalculator_.Calculate_Location(clus_used,hitCollection_p,geometry_p,geometryES_p);
    float theta_s = 2. * atan(exp(-clus_pos.eta()));
    float p0x_s = simple_energy * sin(theta_s) * cos(clus_pos.phi());
    float p0y_s = simple_energy * sin(theta_s) * sin(clus_pos.phi());
    //      float p0z_s = simple_energy * cos(theta_s);
    float et_s = sqrt( p0x_s*p0x_s + p0y_s*p0y_s);

    //std::cout << "       Simple Clustering: E,Et,px,py,pz: "<<simple_energy<<" "<<et_s<<" "<<p0x_s<<" "<<p0y_s<<" "<<std::endl;

    eClus.push_back(simple_energy);
    etClus.push_back(et_s);
    etaClus.push_back(clus_pos.eta());
    phiClus.push_back(clus_pos.phi());
    max_hit.push_back(seed_id);
    RecHitsCluster.push_back(RecHitsInWindow);

    //Compute S4/S9 variable
    float s4s9_[4];
    for(int i=0;i<4;i++)s4s9_[i]= itseed->energy();
    //We are not sure to have 9 RecHits so need to check eta and phi:
    for(unsigned int j=0; j<RecHitsInWindow.size();j++){
      //cout << " Simple cluster rh, ieta, iphi : "<<((EBDetId)RecHitsInWindow[j].id()).ieta()<<" "<<((EBDetId)RecHitsInWindow[j].id()).iphi()<<endl;
      if((((EBDetId)RecHitsInWindow[j].id()).ieta() == seed_id.ieta()-1 && seed_id.ieta()!=1 ) || ( seed_id.ieta()==1 && (((EBDetId)RecHitsInWindow[j].id()).ieta() == seed_id.ieta()-2))){
        if(((EBDetId)RecHitsInWindow[j].id()).iphi() == seed_id.iphi()-1 ||((EBDetId)RecHitsInWindow[j].id()).iphi()-360 == seed_id.iphi()-1 ){
          s4s9_[0]+=RecHitsInWindow[j].energy();
        }else{
          if(((EBDetId)RecHitsInWindow[j].id()).iphi() == seed_id.iphi()){
            s4s9_[0]+=RecHitsInWindow[j].energy();
            s4s9_[1]+=RecHitsInWindow[j].energy();
          }else{
            if(((EBDetId)RecHitsInWindow[j].id()).iphi() == seed_id.iphi()+1 ||((EBDetId)RecHitsInWindow[j].id()).iphi()-360 == seed_id.iphi()+1 ){
              s4s9_[1]+=RecHitsInWindow[j].energy();
            }
          }
        }
      }else{
        if(((EBDetId)RecHitsInWindow[j].id()).ieta() == seed_id.ieta()){
          if(((EBDetId)RecHitsInWindow[j].id()).iphi() == seed_id.iphi()-1 ||((EBDetId)RecHitsInWindow[j].id()).iphi()-360 == seed_id.iphi()-1 ){
            s4s9_[0]+=RecHitsInWindow[j].energy();
            s4s9_[3]+=RecHitsInWindow[j].energy();
          }else{
            if(((EBDetId)RecHitsInWindow[j].id()).iphi() == seed_id.iphi()+1 ||((EBDetId)RecHitsInWindow[j].id()).iphi()-360 == seed_id.iphi()+1 ){
              s4s9_[1]+=RecHitsInWindow[j].energy();
              s4s9_[2]+=RecHitsInWindow[j].energy();
            }
          }
        }else{
          if((((EBDetId)RecHitsInWindow[j].id()).ieta() == seed_id.ieta()+1 && seed_id.ieta()!=-1 ) || ( seed_id.ieta()==-1 && (((EBDetId)RecHitsInWindow[j].id()).ieta() == seed_id.ieta()+2))){
            if(((EBDetId)RecHitsInWindow[j].id()).iphi() == seed_id.iphi()-1 ||((EBDetId)RecHitsInWindow[j].id()).iphi()-360 == seed_id.iphi()-1 ){
              s4s9_[3]+=RecHitsInWindow[j].energy();
            }else{
              if(((EBDetId)RecHitsInWindow[j].id()).iphi() == seed_id.iphi()){
                s4s9_[2]+=RecHitsInWindow[j].energy();
                s4s9_[3]+=RecHitsInWindow[j].energy();
              }else{
                if(((EBDetId)RecHitsInWindow[j].id()).iphi() == seed_id.iphi()+1 ||((EBDetId)RecHitsInWindow[j].id()).iphi()-360 == seed_id.iphi()+1 ){
                  s4s9_[2]+=RecHitsInWindow[j].energy();
                }
              }
            }
          }else{
            cout<<" (EBDetId)RecHitsInWindow[j].id()).ieta() "<<((EBDetId)RecHitsInWindow[j].id()).ieta()<<" seed_id.ieta() "<<seed_id.ieta()<<endl;
            cout<<" Problem with S4 calculation "<<endl;return;
          }
        }
      }
    }
    s4s9Clus.push_back(*max_element( s4s9_,s4s9_+4)/simple_energy);
    //    cout<<" s4s9Clus[0] "<<s4s9_[0]/simple_energy<<" s4s9Clus[1] "<<s4s9_[1]/simple_energy<<" s4s9Clus[2] "<<s4s9_[2]/simple_energy<<" s4s9Clus[3] "<<s4s9_[3]/simple_energy<<endl;
    //    cout<<" Max "<<*max_element( s4s9_,s4s9_+4)/simple_energy<<endl;
    nClus++;
    //std::cout << ">> seed (ieta,iphi,E):" << seed_id.ieta() << ","<<seed_id.iphi()<< std::endl;
    if (nClus == MAXCLUS) return;
  }  //  End loop over seeds to form clusters from rechits around seed

  // cout<< " Pi0 clusters: "<<nClus<<endl;
  //std::cout << ">> nClus:" << nClus << std::endl;
  if (nClus < 2) return;

  //___ Make pi0 candidates by pairing clusters ____________________________________//

  float dR_i, dR_j;
  static const int MAXPI0S = 200;
  int npi0_s=0;
  vector<float> vM0;

  vector<EBDetId> scXtals;
  scXtals.clear();

  unsigned int nPi0Cands = 0;
  std::vector<int> vClusPairIdxs;
  // Clusters sorted by energy so most energetic cands will be formed first
  // Loop over cluster candidate i
  for(Int_t i=0 ; i<nClus ; i++){

    if ( std::find(vClusPairIdxs.begin(),vClusPairIdxs.end(),i) != vClusPairIdxs.end() ) continue;
    // dont attempt to construct more pi0 cands than there are photons to regress
    if ( nPi0Cands >= vRegressPhoIdxs_.size() ) break;

    // Loop over cluster candidate j
    int nJ = 0;
    for(Int_t j=i+1 ; j<nClus ; j++){

      // if clus_j was already used to form a succesful pi0 cand with clus_i,
      // dont attempt to construct another pi0 cand with the same clus_i
      // i.e. choose only most energetic succesful clus pariring to assign as pi0 cand
      if ( nJ > 1 ) break;

      // Do not attempt to pair distant clusters
      dR = reco::deltaR( etaClus[i],phiClus[i], etaClus[j],phiClus[j] );
      if ( dR > 16*0.0174 ) continue;

      //std::cout << ">> etClus[i]:"<<etClus[i] <<",etClus[j]:"<<etClus[j]<<std::endl;
      //if( etClus[i]>selePtGammaOne_ && etClus[j]>selePtGammaTwo_ && s4s9Clus[i]>seleS4S9GammaOne_ && s4s9Clus[j]>seleS4S9GammaTwo_){
      if( etClus[i]>selePtGammaOne_ && etClus[j]>selePtGammaTwo_ ){
        float theta_0 = 2. * atan(exp(-etaClus[i]));
        float theta_1 = 2. * atan(exp(-etaClus[j]));

        float p0x = eClus[i] * sin(theta_0) * cos(phiClus[i]);
        float p1x = eClus[j] * sin(theta_1) * cos(phiClus[j]);
        float p0y = eClus[i] * sin(theta_0) * sin(phiClus[i]);
        float p1y = eClus[j] * sin(theta_1) * sin(phiClus[j]);
        float p0z = eClus[i] * cos(theta_0);
        float p1z = eClus[j] * cos(theta_1);

        float pt_pi0 = sqrt( (p0x+p1x)*(p0x+p1x) + (p0y+p1y)*(p0y+p1y));
        if (pt_pi0 < selePtPi0_)continue;
        float m_inv = sqrt ( (eClus[i] + eClus[j])*(eClus[i] + eClus[j]) - (p0x+p1x)*(p0x+p1x) - (p0y+p1y)*(p0y+p1y) - (p0z+p1z)*(p0z+p1z) );
        //std::cout << ">> m("<<i<<","<<j<<"):" << m_inv << std::endl;
        if ( (m_inv<seleMinvMaxPi0_) && (m_inv>seleMinvMinPi0_) ){

          //New Loop on cluster to measure isolation:
          vector<int> IsoClus;
          IsoClus.clear();
          float Iso = 0;
          TVector3 pi0vect = TVector3((p0x+p1x), (p0y+p1y), (p0z+p1z));
          for(Int_t k=0 ; k<nClus ; k++){
            if(k==i || k==j)continue;
            TVector3 Clusvect = TVector3(eClus[k] * sin(2. * atan(exp(-etaClus[k]))) * cos(phiClus[k]), eClus[k] * sin(2. * atan(exp(-etaClus[k]))) * sin(phiClus[k]) , eClus[k] * cos(2. * atan(exp(-etaClus[k]))));
            float dretaclpi0 = fabs(etaClus[k] - pi0vect.Eta());
            float drclpi0 = Clusvect.DeltaR(pi0vect);

            if((drclpi0<selePi0BeltDR_) && (dretaclpi0<selePi0BeltDeta_) ){

              Iso = Iso + etClus[k];
              IsoClus.push_back(k);
            }
          }
          if(Iso/pt_pi0<selePi0Iso_){

            //vPho_mass3b3_.push_back(m_inv);
            //hMinvPi0EB_->Fill(m_inv);
            //hPt1Pi0EB_->Fill(etClus[i]);
            //hPt2Pi0EB_->Fill(etClus[j]);
            //hPtPi0EB_->Fill(pt_pi0);
            //hIsoPi0EB_->Fill(Iso/pt_pi0);

            npi0_s++;
            //std::cout << ">> m, iso:" << m_inv << std::endl;
          } // iso cut
          if(npi0_s == MAXPI0S) {
            //hNgg_->Fill(npi0_s);
            return;
          }

          // Fill vars
          //std::cout << ">> m("<<i<<","<<j<<"):" << m_inv << std::endl;
          //hMinvPi0EB_noIso->Fill(m_inv);
          //vPho_mass3b3_.push_back( m_inv );

          // Store mass if pi0 cand matched to reco photon
          // Require one of the clusters in pi0 cand to be dR matched
          for ( unsigned int iP = 0; iP < vRegressPhoIdxs_.size(); iP++ ) {
            // Get dR(photon, clus_*)
            dR_i = reco::deltaR( etaClus[i], phiClus[i], vPho_eta_[iP],vPho_phi_[iP] );
            dR_j = reco::deltaR( etaClus[j], phiClus[j], vPho_eta_[iP],vPho_phi_[iP] );
            // At least one of the clusters must be dR matched to photon
            if ( dR_i > 2.*0.04 && dR_j > 2.*0.04 ) continue;
            // Prefer most energetic cand => dont overwrite previous cands
            if ( vPho_mass3b3_[iP] != default_mass ) continue;
            vPho_mass3b3_[iP] = m_inv;
            vPho_dR3b3_[iP] = dR;
            hPho_mass3b3->Fill(m_inv);
            // Store ieta,iphi of 3x3 cluster seeds
            if ( dR_i < dR_j ) {
              // core seed
              ieta_ = max_hit[i].ieta() > 0 ? max_hit[i].ieta()-1 : max_hit[i].ieta(); // [-85,...,-1,1,...,85]
              ieta_ += EBDetId::MAX_IETA; // [0,...,169]
              iphi_ = max_hit[i].iphi()-1; // [0,...,359]
              vPho_ieta03b3_[iP] = ieta_;
              vPho_iphi03b3_[iP] = iphi_;
              // satellite seed
              ieta_ = max_hit[j].ieta() > 0 ? max_hit[j].ieta()-1 : max_hit[j].ieta(); // [-85,...,-1,1,...,85]
              ieta_ += EBDetId::MAX_IETA; // [0,...,169]
              iphi_ = max_hit[j].iphi()-1; // [0,...,359]
              vPho_ieta13b3_[iP] = ieta_;
              vPho_iphi13b3_[iP] = iphi_;
            } else {
              // satellite seed
              ieta_ = max_hit[i].ieta() > 0 ? max_hit[i].ieta()-1 : max_hit[i].ieta(); // [-85,...,-1,1,...,85]
              ieta_ += EBDetId::MAX_IETA; // [0,...,169]
              iphi_ = max_hit[i].iphi()-1; // [0,...,359]
              vPho_ieta13b3_[iP] = ieta_;
              vPho_iphi13b3_[iP] = iphi_;
              // core seed
              ieta_ = max_hit[j].ieta() > 0 ? max_hit[j].ieta()-1 : max_hit[j].ieta(); // [-85,...,-1,1,...,85]
              ieta_ += EBDetId::MAX_IETA; // [0,...,169]
              iphi_ = max_hit[j].iphi()-1; // [0,...,359]
              vPho_ieta03b3_[iP] = ieta_;
              vPho_iphi03b3_[iP] = iphi_;
            }
            //std::cout << "   .. m:"<< m_inv << std::endl;
            nPi0Cands++;
            nJ++;
            break; // dont match to more than one photon, take most energetic photon
          } // LOOP: vRegressPhoIdxs_

          vClusPairIdxs.push_back( i );
          vClusPairIdxs.push_back( j );

        } // IF: m_inv window
      } // IF: pt && s4s9 cut

    } // LOOP: jth cluster
  } // LOOP: ith cluster

  //if ( nPi0Cands == 2 ) hNregress_->Fill(1.);
  //else return;
  //hmA0vmA1->Fill( vPho_mass3b3_[0], vPho_mass3b3_[1] );
  /*
     hmRecovmGen->Fill( vMpi0Gen_[0], vPho_mass3b3_[0] );
     hmRecovmGen->Fill( vMpi0Gen_[1], vPho_mass3b3_[1] );
     */
  //std::cout << ">> N pi0 cands:"<< nPi0Cands << std::endl;
  //std::cout << "   .. m[0]:" << vPho_mass3b3_[0] << " m[1]:" << vPho_mass3b3_[1]  << std::endl;
  //for ( unsigned int iP = 0; iP < vPho_mass3b3.size(); iP++ ) {
  //  vPho_mass3b3_.push_back(vPho_mass3b3[iP]);
  //}
  //hNgg_->Fill(npi0_s);
  //std::cout << ">> nPi0, iso+roi:" << npi0_s << std::endl;
  /*
     if ( npi0_s == 2 ) {
     hMinvPi0EB_1->Fill(vM0[0]);
     hMinvPi0EB_1->Fill(vM0[1]);
     }
     */
  //nEvtsReg++;

} // fill3b3
