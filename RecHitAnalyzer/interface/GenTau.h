#ifndef HepMCCandidate_GenTau_h
#define HepMCCandidate_GenTau_h
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


namespace reco {

  class GenTau : public GenParticle {
  public:
    using GenParticle::GenParticle;
    //GenTau() {}
    //GenTau(Charge q, const LorentzVector &p4, const Point &vtx, int pdgId, int status, bool integerCharge);
    //~GenTau() override;

    //void set_charge_p4(LorentzVector s) { charge_p4_ = s; }

  private:

    LorentzVector charge_p4_;
    LorentzVector neutral_p4_;
    LorentzVector lead_pi0_p4_;
    LorentzVector nu_p4_;

  };

}  // namespace reco

#endif
