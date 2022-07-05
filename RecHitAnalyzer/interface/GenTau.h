#ifndef HepMCCandidate_GenTau_h
#define HepMCCandidate_GenTau_h
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


namespace reco {

  class GenTau : public GenParticle {
  public:
    using GenParticle::GenParticle;

    // set p4 for charged component of tau
    void set_charge_p4(LorentzVector s) { charge_p4_ = s; }
    // set p4 for neutral component of tau (summing together pi0's)
    void set_neutral_p4(LorentzVector s) { neutral_p4_ = s; }
    // set p4 for the leading pi0
    void set_lead_pi0_p4(LorentzVector s) { lead_pi0_p4_ = s; }
    // set p4 for the neutrino
    void set_nu_p4(LorentzVector s) { nu_p4_ = s; }
    // set tau decay mode
    void set_decay_mode(int s) { dm_ = s; }

    // get p4 for charged component of tau
    LorentzVector charge_p4() { return charge_p4_; }
    // get p4 for neutral component of tau (from summing together pi0's)
    LorentzVector neutral_p4() { return neutral_p4_; }
    // get p4 for the leading pi0
    LorentzVector lead_pi0_p4() { return lead_pi0_p4_; }
    // get p4 for the neutrino
    LorentzVector nu_p4() { return nu_p4_; }
    // get p4 for the visible component of the tau
    LorentzVector vis_p4() { return this->p4()-nu_p4_; }
    // get tau decay mode
    int decay_mode() { return dm_; }


  private:

    LorentzVector charge_p4_ = LorentzVector(0,0,0,0);
    LorentzVector neutral_p4_ = LorentzVector(0,0,0,0);
    LorentzVector lead_pi0_p4_= LorentzVector(0,0,0,0);
    LorentzVector nu_p4_= LorentzVector(0,0,0,0);
    int dm_=-1;

  };

}  // namespace reco

#endif
