#include <string>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

class MtdEleIsoHarvester : public DQMEDHarvester {
public:
  explicit MtdEleIsoHarvester(const edm::ParameterSet& iConfig);
  ~MtdEleIsoHarvester() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

protected:
  void dqmEndJob(DQMStore::IBooker&, DQMStore::IGetter&) override;

private:
  void computeEfficiency1D(MonitorElement* num, MonitorElement* den, MonitorElement* result);

  const std::string folder_;

  // --- Histograms
  MonitorElement* mePtEffnoMTD_Sig_EB_;
  MonitorElement* meEtaEffnoMTD_Sig_EB_;
  MonitorElement* mePhiEffnoMTD_Sig_EB_;
  
  MonitorElement* mePtEffMTD_1_Sig_EB_;
  MonitorElement* meEtaEffMTD_1_Sig_EB_;
  MonitorElement* mePhiEffMTD_1_Sig_EB_;

  MonitorElement* mePtEffMTD_2_Sig_EB_;
  MonitorElement* meEtaEffMTD_2_Sig_EB_;
  MonitorElement* mePhiEffMTD_2_Sig_EB_;

  MonitorElement* mePtEffMTD_3_Sig_EB_;
  MonitorElement* meEtaEffMTD_3_Sig_EB_;
  MonitorElement* mePhiEffMTD_3_Sig_EB_;

  MonitorElement* mePtEffMTD_4_Sig_EB_;
  MonitorElement* meEtaEffMTD_4_Sig_EB_;
  MonitorElement* mePhiEffMTD_4_Sig_EB_;

  MonitorElement* mePtEffMTD_5_Sig_EB_;
  MonitorElement* meEtaEffMTD_5_Sig_EB_;
  MonitorElement* mePhiEffMTD_5_Sig_EB_;

  MonitorElement* mePtEffMTD_6_Sig_EB_;
  MonitorElement* meEtaEffMTD_6_Sig_EB_;
  MonitorElement* mePhiEffMTD_6_Sig_EB_;

  MonitorElement* mePtEffMTD_7_Sig_EB_;
  MonitorElement* meEtaEffMTD_7_Sig_EB_;
  MonitorElement* mePhiEffMTD_7_Sig_EB_;


  MonitorElement* mePtEffnoMTD_Sig_EE_;
  MonitorElement* meEtaEffnoMTD_Sig_EE_;
  MonitorElement* mePhiEffnoMTD_Sig_EE_;
  
  MonitorElement* mePtEffMTD_1_Sig_EE_;
  MonitorElement* meEtaEffMTD_1_Sig_EE_;
  MonitorElement* mePhiEffMTD_1_Sig_EE_;

  MonitorElement* mePtEffMTD_2_Sig_EE_;
  MonitorElement* meEtaEffMTD_2_Sig_EE_;
  MonitorElement* mePhiEffMTD_2_Sig_EE_;

  MonitorElement* mePtEffMTD_3_Sig_EE_;
  MonitorElement* meEtaEffMTD_3_Sig_EE_;
  MonitorElement* mePhiEffMTD_3_Sig_EE_;

  MonitorElement* mePtEffMTD_4_Sig_EE_;
  MonitorElement* meEtaEffMTD_4_Sig_EE_;
  MonitorElement* mePhiEffMTD_4_Sig_EE_;

  MonitorElement* mePtEffMTD_5_Sig_EE_;
  MonitorElement* meEtaEffMTD_5_Sig_EE_;
  MonitorElement* mePhiEffMTD_5_Sig_EE_;

  MonitorElement* mePtEffMTD_6_Sig_EE_;
  MonitorElement* meEtaEffMTD_6_Sig_EE_;
  MonitorElement* mePhiEffMTD_6_Sig_EE_;

  MonitorElement* mePtEffMTD_7_Sig_EE_;
  MonitorElement* meEtaEffMTD_7_Sig_EE_;
  MonitorElement* mePhiEffMTD_7_Sig_EE_;




  MonitorElement* mePtEffnoMTD_Bkg_EB_;
  MonitorElement* meEtaEffnoMTD_Bkg_EB_;
  MonitorElement* mePhiEffnoMTD_Bkg_EB_;
  
  MonitorElement* mePtEffMTD_1_Bkg_EB_;
  MonitorElement* meEtaEffMTD_1_Bkg_EB_;
  MonitorElement* mePhiEffMTD_1_Bkg_EB_;

  MonitorElement* mePtEffMTD_2_Bkg_EB_;
  MonitorElement* meEtaEffMTD_2_Bkg_EB_;
  MonitorElement* mePhiEffMTD_2_Bkg_EB_;

  MonitorElement* mePtEffMTD_3_Bkg_EB_;
  MonitorElement* meEtaEffMTD_3_Bkg_EB_;
  MonitorElement* mePhiEffMTD_3_Bkg_EB_;

  MonitorElement* mePtEffMTD_4_Bkg_EB_;
  MonitorElement* meEtaEffMTD_4_Bkg_EB_;
  MonitorElement* mePhiEffMTD_4_Bkg_EB_;

  MonitorElement* mePtEffMTD_5_Bkg_EB_;
  MonitorElement* meEtaEffMTD_5_Bkg_EB_;
  MonitorElement* mePhiEffMTD_5_Bkg_EB_;

  MonitorElement* mePtEffMTD_6_Bkg_EB_;
  MonitorElement* meEtaEffMTD_6_Bkg_EB_;
  MonitorElement* mePhiEffMTD_6_Bkg_EB_;

  MonitorElement* mePtEffMTD_7_Bkg_EB_;
  MonitorElement* meEtaEffMTD_7_Bkg_EB_;
  MonitorElement* mePhiEffMTD_7_Bkg_EB_;


  MonitorElement* mePtEffnoMTD_Bkg_EE_;
  MonitorElement* meEtaEffnoMTD_Bkg_EE_;
  MonitorElement* mePhiEffnoMTD_Bkg_EE_;
  
  MonitorElement* mePtEffMTD_1_Bkg_EE_;
  MonitorElement* meEtaEffMTD_1_Bkg_EE_;
  MonitorElement* mePhiEffMTD_1_Bkg_EE_;

  MonitorElement* mePtEffMTD_2_Bkg_EE_;
  MonitorElement* meEtaEffMTD_2_Bkg_EE_;
  MonitorElement* mePhiEffMTD_2_Bkg_EE_;

  MonitorElement* mePtEffMTD_3_Bkg_EE_;
  MonitorElement* meEtaEffMTD_3_Bkg_EE_;
  MonitorElement* mePhiEffMTD_3_Bkg_EE_;

  MonitorElement* mePtEffMTD_4_Bkg_EE_;
  MonitorElement* meEtaEffMTD_4_Bkg_EE_;
  MonitorElement* mePhiEffMTD_4_Bkg_EE_;

  MonitorElement* mePtEffMTD_5_Bkg_EE_;
  MonitorElement* meEtaEffMTD_5_Bkg_EE_;
  MonitorElement* mePhiEffMTD_5_Bkg_EE_;

  MonitorElement* mePtEffMTD_6_Bkg_EE_;
  MonitorElement* meEtaEffMTD_6_Bkg_EE_;
  MonitorElement* mePhiEffMTD_6_Bkg_EE_;

  MonitorElement* mePtEffMTD_7_Bkg_EE_;
  MonitorElement* meEtaEffMTD_7_Bkg_EE_;
  MonitorElement* mePhiEffMTD_7_Bkg_EE_;


};

// ------------ constructor and destructor --------------
MtdEleIsoHarvester::MtdEleIsoHarvester(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")) {}

MtdEleIsoHarvester::~MtdEleIsoHarvester() {}

// auxiliary method to compute efficiency from the ratio of two 1D MonitorElement
void MtdEleIsoHarvester::computeEfficiency1D(MonitorElement* num, MonitorElement* den, MonitorElement* result) {
  for (int ibin = 1; ibin <= den->getNbinsX(); ibin++) {
    double eff = num->getBinContent(ibin) / den->getBinContent(ibin);
    double bin_err = sqrt((num->getBinContent(ibin) * (den->getBinContent(ibin) - num->getBinContent(ibin))) /
                          pow(den->getBinContent(ibin), 3));
    if (den->getBinContent(ibin) == 0) {
      eff = 0;
      bin_err = 0;
    }
    result->setBinContent(ibin, eff);
    result->setBinError(ibin, bin_err);
  }
}

// ------------ endjob tasks ----------------------------
void MtdEleIsoHarvester::dqmEndJob(DQMStore::IBooker& ibook, DQMStore::IGetter& igetter) {
  // --- Get the monitoring histograms

  // For promt (signal)
  MonitorElement* meEle_pt_tot_Sig_EB_ = igetter.get(folder_ + "Ele_pT_tot_Sig_EB");
  MonitorElement* meEle_pt_noMTD_Sig_EB_ = igetter.get(folder_ + "Ele_pT_noMTD_Sig_EB");

  MonitorElement* meEle_eta_tot_Sig_EB_ = igetter.get(folder_ + "Ele_eta_tot_Sig_EB");
  MonitorElement* meEle_eta_noMTD_Sig_EB_ = igetter.get(folder_ + "Ele_eta_noMTD_Sig_EB");

  MonitorElement* meEle_phi_tot_Sig_EB_ = igetter.get(folder_ + "Ele_phi_tot_Sig_EB");
  MonitorElement* meEle_phi_noMTD_Sig_EB_ = igetter.get(folder_ + "Ele_phi_noMTD_Sig_EB");

  MonitorElement* meEle_pt_MTD_1_Sig_EB_ = igetter.get(folder_ + "Ele_pT_MTD_1_Sig_EB");
  MonitorElement* meEle_eta_MTD_1_Sig_EB_ = igetter.get(folder_ + "Ele_eta_MTD_1_Sig_EB");
  MonitorElement* meEle_phi_MTD_1_Sig_EB_ = igetter.get(folder_ + "Ele_phi_MTD_1_Sig_EB");

  MonitorElement* meEle_pt_MTD_2_Sig_EB_ = igetter.get(folder_ + "Ele_pT_MTD_2_Sig_EB");
  MonitorElement* meEle_eta_MTD_2_Sig_EB_ = igetter.get(folder_ + "Ele_eta_MTD_2_Sig_EB");
  MonitorElement* meEle_phi_MTD_2_Sig_EB_ = igetter.get(folder_ + "Ele_phi_MTD_2_Sig_EB");

  MonitorElement* meEle_pt_MTD_3_Sig_EB_ = igetter.get(folder_ + "Ele_pT_MTD_3_Sig_EB");
  MonitorElement* meEle_eta_MTD_3_Sig_EB_ = igetter.get(folder_ + "Ele_eta_MTD_3_Sig_EB");
  MonitorElement* meEle_phi_MTD_3_Sig_EB_ = igetter.get(folder_ + "Ele_phi_MTD_3_Sig_EB");

  MonitorElement* meEle_pt_MTD_4_Sig_EB_ = igetter.get(folder_ + "Ele_pT_MTD_4_Sig_EB");
  MonitorElement* meEle_eta_MTD_4_Sig_EB_ = igetter.get(folder_ + "Ele_eta_MTD_4_Sig_EB");
  MonitorElement* meEle_phi_MTD_4_Sig_EB_ = igetter.get(folder_ + "Ele_phi_MTD_4_Sig_EB");

  MonitorElement* meEle_pt_MTD_5_Sig_EB_ = igetter.get(folder_ + "Ele_pT_MTD_5_Sig_EB");
  MonitorElement* meEle_eta_MTD_5_Sig_EB_ = igetter.get(folder_ + "Ele_eta_MTD_5_Sig_EB");
  MonitorElement* meEle_phi_MTD_5_Sig_EB_ = igetter.get(folder_ + "Ele_phi_MTD_5_Sig_EB");

  MonitorElement* meEle_pt_MTD_6_Sig_EB_ = igetter.get(folder_ + "Ele_pT_MTD_6_Sig_EB");
  MonitorElement* meEle_eta_MTD_6_Sig_EB_ = igetter.get(folder_ + "Ele_eta_MTD_6_Sig_EB");
  MonitorElement* meEle_phi_MTD_6_Sig_EB_ = igetter.get(folder_ + "Ele_phi_MTD_6_Sig_EB");

  MonitorElement* meEle_pt_MTD_7_Sig_EB_ = igetter.get(folder_ + "Ele_pT_MTD_7_Sig_EB");
  MonitorElement* meEle_eta_MTD_7_Sig_EB_ = igetter.get(folder_ + "Ele_eta_MTD_7_Sig_EB");
  MonitorElement* meEle_phi_MTD_7_Sig_EB_ = igetter.get(folder_ + "Ele_phi_MTD_7_Sig_EB");


  MonitorElement* meEle_pt_tot_Sig_EE_ = igetter.get(folder_ + "Ele_pT_tot_Sig_EE");
  MonitorElement* meEle_pt_noMTD_Sig_EE_ = igetter.get(folder_ + "Ele_pT_noMTD_Sig_EE");

  MonitorElement* meEle_eta_tot_Sig_EE_ = igetter.get(folder_ + "Ele_eta_tot_Sig_EE");
  MonitorElement* meEle_eta_noMTD_Sig_EE_ = igetter.get(folder_ + "Ele_eta_noMTD_Sig_EE");

  MonitorElement* meEle_phi_tot_Sig_EE_ = igetter.get(folder_ + "Ele_phi_tot_Sig_EE");
  MonitorElement* meEle_phi_noMTD_Sig_EE_ = igetter.get(folder_ + "Ele_phi_noMTD_Sig_EE");

  MonitorElement* meEle_pt_MTD_1_Sig_EE_ = igetter.get(folder_ + "Ele_pT_MTD_1_Sig_EE");
  MonitorElement* meEle_eta_MTD_1_Sig_EE_ = igetter.get(folder_ + "Ele_eta_MTD_1_Sig_EE");
  MonitorElement* meEle_phi_MTD_1_Sig_EE_ = igetter.get(folder_ + "Ele_phi_MTD_1_Sig_EE");

  MonitorElement* meEle_pt_MTD_2_Sig_EE_ = igetter.get(folder_ + "Ele_pT_MTD_2_Sig_EE");
  MonitorElement* meEle_eta_MTD_2_Sig_EE_ = igetter.get(folder_ + "Ele_eta_MTD_2_Sig_EE");
  MonitorElement* meEle_phi_MTD_2_Sig_EE_ = igetter.get(folder_ + "Ele_phi_MTD_2_Sig_EE");

  MonitorElement* meEle_pt_MTD_3_Sig_EE_ = igetter.get(folder_ + "Ele_pT_MTD_3_Sig_EE");
  MonitorElement* meEle_eta_MTD_3_Sig_EE_ = igetter.get(folder_ + "Ele_eta_MTD_3_Sig_EE");
  MonitorElement* meEle_phi_MTD_3_Sig_EE_ = igetter.get(folder_ + "Ele_phi_MTD_3_Sig_EE");

  MonitorElement* meEle_pt_MTD_4_Sig_EE_ = igetter.get(folder_ + "Ele_pT_MTD_4_Sig_EE");
  MonitorElement* meEle_eta_MTD_4_Sig_EE_ = igetter.get(folder_ + "Ele_eta_MTD_4_Sig_EE");
  MonitorElement* meEle_phi_MTD_4_Sig_EE_ = igetter.get(folder_ + "Ele_phi_MTD_4_Sig_EE");

  MonitorElement* meEle_pt_MTD_5_Sig_EE_ = igetter.get(folder_ + "Ele_pT_MTD_5_Sig_EE");
  MonitorElement* meEle_eta_MTD_5_Sig_EE_ = igetter.get(folder_ + "Ele_eta_MTD_5_Sig_EE");
  MonitorElement* meEle_phi_MTD_5_Sig_EE_ = igetter.get(folder_ + "Ele_phi_MTD_5_Sig_EE");

  MonitorElement* meEle_pt_MTD_6_Sig_EE_ = igetter.get(folder_ + "Ele_pT_MTD_6_Sig_EE");
  MonitorElement* meEle_eta_MTD_6_Sig_EE_ = igetter.get(folder_ + "Ele_eta_MTD_6_Sig_EE");
  MonitorElement* meEle_phi_MTD_6_Sig_EE_ = igetter.get(folder_ + "Ele_phi_MTD_6_Sig_EE");

  MonitorElement* meEle_pt_MTD_7_Sig_EE_ = igetter.get(folder_ + "Ele_pT_MTD_7_Sig_EE");
  MonitorElement* meEle_eta_MTD_7_Sig_EE_ = igetter.get(folder_ + "Ele_eta_MTD_7_Sig_EE");
  MonitorElement* meEle_phi_MTD_7_Sig_EE_ = igetter.get(folder_ + "Ele_phi_MTD_7_Sig_EE");





  // for non-promt (background)
  MonitorElement* meEle_pt_tot_Bkg_EB_ = igetter.get(folder_ + "Ele_pT_tot_Bkg_EB");
  MonitorElement* meEle_pt_noMTD_Bkg_EB_ = igetter.get(folder_ + "Ele_pT_noMTD_Bkg_EB");

  MonitorElement* meEle_eta_tot_Bkg_EB_ = igetter.get(folder_ + "Ele_eta_tot_Bkg_EB");
  MonitorElement* meEle_eta_noMTD_Bkg_EB_ = igetter.get(folder_ + "Ele_eta_noMTD_Bkg_EB");

  MonitorElement* meEle_phi_tot_Bkg_EB_ = igetter.get(folder_ + "Ele_phi_tot_Bkg_EB");
  MonitorElement* meEle_phi_noMTD_Bkg_EB_ = igetter.get(folder_ + "Ele_phi_noMTD_Bkg_EB");

  MonitorElement* meEle_pt_MTD_1_Bkg_EB_ = igetter.get(folder_ + "Ele_pT_MTD_1_Bkg_EB");
  MonitorElement* meEle_eta_MTD_1_Bkg_EB_ = igetter.get(folder_ + "Ele_eta_MTD_1_Bkg_EB");
  MonitorElement* meEle_phi_MTD_1_Bkg_EB_ = igetter.get(folder_ + "Ele_phi_MTD_1_Bkg_EB");

  MonitorElement* meEle_pt_MTD_2_Bkg_EB_ = igetter.get(folder_ + "Ele_pT_MTD_2_Bkg_EB");
  MonitorElement* meEle_eta_MTD_2_Bkg_EB_ = igetter.get(folder_ + "Ele_eta_MTD_2_Bkg_EB");
  MonitorElement* meEle_phi_MTD_2_Bkg_EB_ = igetter.get(folder_ + "Ele_phi_MTD_2_Bkg_EB");

  MonitorElement* meEle_pt_MTD_3_Bkg_EB_ = igetter.get(folder_ + "Ele_pT_MTD_3_Bkg_EB");
  MonitorElement* meEle_eta_MTD_3_Bkg_EB_ = igetter.get(folder_ + "Ele_eta_MTD_3_Bkg_EB");
  MonitorElement* meEle_phi_MTD_3_Bkg_EB_ = igetter.get(folder_ + "Ele_phi_MTD_3_Bkg_EB");

  MonitorElement* meEle_pt_MTD_4_Bkg_EB_ = igetter.get(folder_ + "Ele_pT_MTD_4_Bkg_EB");
  MonitorElement* meEle_eta_MTD_4_Bkg_EB_ = igetter.get(folder_ + "Ele_eta_MTD_4_Bkg_EB");
  MonitorElement* meEle_phi_MTD_4_Bkg_EB_ = igetter.get(folder_ + "Ele_phi_MTD_4_Bkg_EB");

  MonitorElement* meEle_pt_MTD_5_Bkg_EB_ = igetter.get(folder_ + "Ele_pT_MTD_5_Bkg_EB");
  MonitorElement* meEle_eta_MTD_5_Bkg_EB_ = igetter.get(folder_ + "Ele_eta_MTD_5_Bkg_EB");
  MonitorElement* meEle_phi_MTD_5_Bkg_EB_ = igetter.get(folder_ + "Ele_phi_MTD_5_Bkg_EB");

  MonitorElement* meEle_pt_MTD_6_Bkg_EB_ = igetter.get(folder_ + "Ele_pT_MTD_6_Bkg_EB");
  MonitorElement* meEle_eta_MTD_6_Bkg_EB_ = igetter.get(folder_ + "Ele_eta_MTD_6_Bkg_EB");
  MonitorElement* meEle_phi_MTD_6_Bkg_EB_ = igetter.get(folder_ + "Ele_phi_MTD_6_Bkg_EB");

  MonitorElement* meEle_pt_MTD_7_Bkg_EB_ = igetter.get(folder_ + "Ele_pT_MTD_7_Bkg_EB");
  MonitorElement* meEle_eta_MTD_7_Bkg_EB_ = igetter.get(folder_ + "Ele_eta_MTD_7_Bkg_EB");
  MonitorElement* meEle_phi_MTD_7_Bkg_EB_ = igetter.get(folder_ + "Ele_phi_MTD_7_Bkg_EB");


  MonitorElement* meEle_pt_tot_Bkg_EE_ = igetter.get(folder_ + "Ele_pT_tot_Bkg_EE");
  MonitorElement* meEle_pt_noMTD_Bkg_EE_ = igetter.get(folder_ + "Ele_pT_noMTD_Bkg_EE");

  MonitorElement* meEle_eta_tot_Bkg_EE_ = igetter.get(folder_ + "Ele_eta_tot_Bkg_EE");
  MonitorElement* meEle_eta_noMTD_Bkg_EE_ = igetter.get(folder_ + "Ele_eta_noMTD_Bkg_EE");

  MonitorElement* meEle_phi_tot_Bkg_EE_ = igetter.get(folder_ + "Ele_phi_tot_Bkg_EE");
  MonitorElement* meEle_phi_noMTD_Bkg_EE_ = igetter.get(folder_ + "Ele_phi_noMTD_Bkg_EE");

  MonitorElement* meEle_pt_MTD_1_Bkg_EE_ = igetter.get(folder_ + "Ele_pT_MTD_1_Bkg_EE");
  MonitorElement* meEle_eta_MTD_1_Bkg_EE_ = igetter.get(folder_ + "Ele_eta_MTD_1_Bkg_EE");
  MonitorElement* meEle_phi_MTD_1_Bkg_EE_ = igetter.get(folder_ + "Ele_phi_MTD_1_Bkg_EE");

  MonitorElement* meEle_pt_MTD_2_Bkg_EE_ = igetter.get(folder_ + "Ele_pT_MTD_2_Bkg_EE");
  MonitorElement* meEle_eta_MTD_2_Bkg_EE_ = igetter.get(folder_ + "Ele_eta_MTD_2_Bkg_EE");
  MonitorElement* meEle_phi_MTD_2_Bkg_EE_ = igetter.get(folder_ + "Ele_phi_MTD_2_Bkg_EE");

  MonitorElement* meEle_pt_MTD_3_Bkg_EE_ = igetter.get(folder_ + "Ele_pT_MTD_3_Bkg_EE");
  MonitorElement* meEle_eta_MTD_3_Bkg_EE_ = igetter.get(folder_ + "Ele_eta_MTD_3_Bkg_EE");
  MonitorElement* meEle_phi_MTD_3_Bkg_EE_ = igetter.get(folder_ + "Ele_phi_MTD_3_Bkg_EE");

  MonitorElement* meEle_pt_MTD_4_Bkg_EE_ = igetter.get(folder_ + "Ele_pT_MTD_4_Bkg_EE");
  MonitorElement* meEle_eta_MTD_4_Bkg_EE_ = igetter.get(folder_ + "Ele_eta_MTD_4_Bkg_EE");
  MonitorElement* meEle_phi_MTD_4_Bkg_EE_ = igetter.get(folder_ + "Ele_phi_MTD_4_Bkg_EE");

  MonitorElement* meEle_pt_MTD_5_Bkg_EE_ = igetter.get(folder_ + "Ele_pT_MTD_5_Bkg_EE");
  MonitorElement* meEle_eta_MTD_5_Bkg_EE_ = igetter.get(folder_ + "Ele_eta_MTD_5_Bkg_EE");
  MonitorElement* meEle_phi_MTD_5_Bkg_EE_ = igetter.get(folder_ + "Ele_phi_MTD_5_Bkg_EE");

  MonitorElement* meEle_pt_MTD_6_Bkg_EE_ = igetter.get(folder_ + "Ele_pT_MTD_6_Bkg_EE");
  MonitorElement* meEle_eta_MTD_6_Bkg_EE_ = igetter.get(folder_ + "Ele_eta_MTD_6_Bkg_EE");
  MonitorElement* meEle_phi_MTD_6_Bkg_EE_ = igetter.get(folder_ + "Ele_phi_MTD_6_Bkg_EE");

  MonitorElement* meEle_pt_MTD_7_Bkg_EE_ = igetter.get(folder_ + "Ele_pT_MTD_7_Bkg_EE");
  MonitorElement* meEle_eta_MTD_7_Bkg_EE_ = igetter.get(folder_ + "Ele_eta_MTD_7_Bkg_EE");
  MonitorElement* meEle_phi_MTD_7_Bkg_EE_ = igetter.get(folder_ + "Ele_phi_MTD_7_Bkg_EE");



  if (!meEle_pt_tot_Sig_EB_ || !meEle_pt_MTD_1_Sig_EB_ || !meEle_pt_MTD_2_Sig_EB_ ||
      !meEle_pt_MTD_3_Sig_EB_ || !meEle_pt_MTD_4_Sig_EB_ || !meEle_pt_MTD_5_Sig_EB_ || !meEle_pt_MTD_6_Sig_EB_ || !meEle_pt_MTD_7_Sig_EB_ ||
      !meEle_pt_noMTD_Sig_EB_ || !meEle_eta_tot_Sig_EB_ || !meEle_eta_MTD_1_Sig_EB_ || !meEle_eta_MTD_2_Sig_EB_ || !meEle_eta_MTD_3_Sig_EB_ ||
      !meEle_eta_MTD_4_Sig_EB_ || !meEle_eta_MTD_5_Sig_EB_ || !meEle_eta_MTD_6_Sig_EB_ || !meEle_eta_MTD_7_Sig_EB_ || !meEle_eta_noMTD_Sig_EB_ || 
      !meEle_phi_tot_Sig_EB_ || !meEle_phi_MTD_1_Sig_EB_ || !meEle_phi_MTD_2_Sig_EB_ || !meEle_phi_MTD_3_Sig_EB_ || !meEle_phi_MTD_4_Sig_EB_ ||
      !meEle_phi_MTD_5_Sig_EB_ || !meEle_phi_MTD_6_Sig_EB_ || !meEle_phi_MTD_7_Sig_EB_ || !meEle_phi_noMTD_Sig_EB_ || 
      !meEle_pt_tot_Sig_EE_ || !meEle_pt_MTD_1_Sig_EE_ || !meEle_pt_MTD_2_Sig_EE_ ||
      !meEle_pt_MTD_3_Sig_EE_ || !meEle_pt_MTD_4_Sig_EE_ || !meEle_pt_MTD_5_Sig_EE_ || !meEle_pt_MTD_6_Sig_EE_ || !meEle_pt_MTD_7_Sig_EE_ ||
      !meEle_pt_noMTD_Sig_EE_ || !meEle_eta_tot_Sig_EE_ || !meEle_eta_MTD_1_Sig_EE_ || !meEle_eta_MTD_2_Sig_EE_ || !meEle_eta_MTD_3_Sig_EE_ ||
      !meEle_eta_MTD_4_Sig_EE_ || !meEle_eta_MTD_5_Sig_EE_ || !meEle_eta_MTD_6_Sig_EE_ || !meEle_eta_MTD_7_Sig_EE_ || !meEle_eta_noMTD_Sig_EE_ || 
      !meEle_phi_tot_Sig_EE_ || !meEle_phi_MTD_1_Sig_EE_ || !meEle_phi_MTD_2_Sig_EE_ || !meEle_phi_MTD_3_Sig_EE_ || !meEle_phi_MTD_4_Sig_EE_ ||
      !meEle_phi_MTD_5_Sig_EE_ || !meEle_phi_MTD_6_Sig_EE_ || !meEle_phi_MTD_7_Sig_EE_ || !meEle_phi_noMTD_Sig_EE_ ||
      !meEle_pt_tot_Bkg_EB_ || !meEle_pt_MTD_1_Bkg_EB_ || !meEle_pt_MTD_2_Bkg_EB_ ||
      !meEle_pt_MTD_3_Bkg_EB_ || !meEle_pt_MTD_4_Bkg_EB_ || !meEle_pt_MTD_5_Bkg_EB_ || !meEle_pt_MTD_6_Bkg_EB_ || !meEle_pt_MTD_7_Bkg_EB_ ||
      !meEle_pt_noMTD_Bkg_EB_ || !meEle_eta_tot_Bkg_EB_ || !meEle_eta_MTD_1_Bkg_EB_ || !meEle_eta_MTD_2_Bkg_EB_ || !meEle_eta_MTD_3_Bkg_EB_ ||
      !meEle_eta_MTD_4_Bkg_EB_ || !meEle_eta_MTD_5_Bkg_EB_ || !meEle_eta_MTD_6_Bkg_EB_ || !meEle_eta_MTD_7_Bkg_EB_ || !meEle_eta_noMTD_Bkg_EB_ || 
      !meEle_phi_tot_Bkg_EB_ || !meEle_phi_MTD_1_Bkg_EB_ || !meEle_phi_MTD_2_Bkg_EB_ || !meEle_phi_MTD_3_Bkg_EB_ || !meEle_phi_MTD_4_Bkg_EB_ ||
      !meEle_phi_MTD_5_Bkg_EB_ || !meEle_phi_MTD_6_Bkg_EB_ || !meEle_phi_MTD_7_Bkg_EB_ || !meEle_phi_noMTD_Bkg_EB_ || 
      !meEle_pt_tot_Bkg_EE_ || !meEle_pt_MTD_1_Bkg_EE_ || !meEle_pt_MTD_2_Bkg_EE_ ||
      !meEle_pt_MTD_3_Bkg_EE_ || !meEle_pt_MTD_4_Bkg_EE_ || !meEle_pt_MTD_5_Bkg_EE_ || !meEle_pt_MTD_6_Bkg_EE_ || !meEle_pt_MTD_7_Bkg_EE_ ||
      !meEle_pt_noMTD_Bkg_EE_ || !meEle_eta_tot_Bkg_EE_ || !meEle_eta_MTD_1_Bkg_EE_ || !meEle_eta_MTD_2_Bkg_EE_ || !meEle_eta_MTD_3_Bkg_EE_ ||
      !meEle_eta_MTD_4_Bkg_EE_ || !meEle_eta_MTD_5_Bkg_EE_ || !meEle_eta_MTD_6_Bkg_EE_ || !meEle_eta_MTD_7_Bkg_EE_ || !meEle_eta_noMTD_Bkg_EE_ || 
      !meEle_phi_tot_Bkg_EE_ || !meEle_phi_MTD_1_Bkg_EE_ || !meEle_phi_MTD_2_Bkg_EE_ || !meEle_phi_MTD_3_Bkg_EE_ || !meEle_phi_MTD_4_Bkg_EE_ ||
      !meEle_phi_MTD_5_Bkg_EE_ || !meEle_phi_MTD_6_Bkg_EE_ || !meEle_phi_MTD_7_Bkg_EE_ || !meEle_phi_noMTD_Bkg_EE_) {
    edm::LogError("MtdEleIsoHarvester") << "Monitoring histograms not found!" << std::endl;
    return;
  }

  // --- Book  histograms
  ibook.cd(folder_);
  // ele iso addition starts; MTD vs noMTD case
  /////////////////////////////////////////////////////// For promt (signal)

  mePtEffMTD_1_Sig_EB_ = ibook.book1D("pTeffMTD_1_Sig_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EB_->getNbinsX(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_1_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_1_Sig_EB_, meEle_pt_tot_Sig_EB_, mePtEffMTD_1_Sig_EB_);

  mePtEffMTD_2_Sig_EB_ = ibook.book1D("pTeffMTD_2_Sig_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EB_->getNbinsX(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_2_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_2_Sig_EB_, meEle_pt_tot_Sig_EB_, mePtEffMTD_2_Sig_EB_);

  mePtEffMTD_3_Sig_EB_ = ibook.book1D("pTeffMTD_3_Sig_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EB_->getNbinsX(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_3_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_3_Sig_EB_, meEle_pt_tot_Sig_EB_, mePtEffMTD_3_Sig_EB_);

  mePtEffMTD_4_Sig_EB_ = ibook.book1D("pTeffMTD_4_Sig_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EB_->getNbinsX(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_4_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_4_Sig_EB_, meEle_pt_tot_Sig_EB_, mePtEffMTD_4_Sig_EB_);

  mePtEffMTD_5_Sig_EB_ = ibook.book1D("pTeffMTD_5_Sig_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EB_->getNbinsX(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_5_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_5_Sig_EB_, meEle_pt_tot_Sig_EB_, mePtEffMTD_5_Sig_EB_);

  mePtEffMTD_6_Sig_EB_ = ibook.book1D("pTeffMTD_6_Sig_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EB_->getNbinsX(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_6_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_6_Sig_EB_, meEle_pt_tot_Sig_EB_, mePtEffMTD_6_Sig_EB_);

  mePtEffMTD_7_Sig_EB_ = ibook.book1D("pTeffMTD_7_Sig_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EB_->getNbinsX(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_7_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_7_Sig_EB_, meEle_pt_tot_Sig_EB_, mePtEffMTD_7_Sig_EB_);






  meEtaEffMTD_1_Sig_EB_ = ibook.book1D("EtaEffMTD_1_Sig_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EB_->getNbinsX(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_1_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_1_Sig_EB_, meEle_eta_tot_Sig_EB_, meEtaEffMTD_1_Sig_EB_);

  meEtaEffMTD_2_Sig_EB_ = ibook.book1D("EtaEffMTD_2_Sig_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EB_->getNbinsX(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_2_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_2_Sig_EB_, meEle_eta_tot_Sig_EB_, meEtaEffMTD_2_Sig_EB_);

  meEtaEffMTD_3_Sig_EB_ = ibook.book1D("EtaEffMTD_3_Sig_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EB_->getNbinsX(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_3_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_3_Sig_EB_, meEle_eta_tot_Sig_EB_, meEtaEffMTD_3_Sig_EB_);

  meEtaEffMTD_4_Sig_EB_ = ibook.book1D("EtaEffMTD_4_Sig_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EB_->getNbinsX(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_4_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_4_Sig_EB_, meEle_eta_tot_Sig_EB_, meEtaEffMTD_4_Sig_EB_);

  meEtaEffMTD_5_Sig_EB_ = ibook.book1D("EtaEffMTD_5_Sig_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EB_->getNbinsX(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_5_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_5_Sig_EB_, meEle_eta_tot_Sig_EB_, meEtaEffMTD_5_Sig_EB_);

  meEtaEffMTD_6_Sig_EB_ = ibook.book1D("EtaEffMTD_6_Sig_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EB_->getNbinsX(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_6_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_6_Sig_EB_, meEle_eta_tot_Sig_EB_, meEtaEffMTD_6_Sig_EB_);

  meEtaEffMTD_7_Sig_EB_ = ibook.book1D("EtaEffMTD_7_Sig_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EB_->getNbinsX(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_7_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_7_Sig_EB_, meEle_eta_tot_Sig_EB_, meEtaEffMTD_7_Sig_EB_);





  mePhiEffMTD_1_Sig_EB_ = ibook.book1D("PhiEffMTD_1_Sig_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EB_->getNbinsX(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_1_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_1_Sig_EB_, meEle_phi_tot_Sig_EB_, mePhiEffMTD_1_Sig_EB_);

  mePhiEffMTD_2_Sig_EB_ = ibook.book1D("PhiEffMTD_2_Sig_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EB_->getNbinsX(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_2_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_2_Sig_EB_, meEle_phi_tot_Sig_EB_, mePhiEffMTD_2_Sig_EB_);

  mePhiEffMTD_3_Sig_EB_ = ibook.book1D("PhiEffMTD_3_Sig_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EB_->getNbinsX(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_3_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_3_Sig_EB_, meEle_phi_tot_Sig_EB_, mePhiEffMTD_3_Sig_EB_);

  mePhiEffMTD_4_Sig_EB_ = ibook.book1D("PhiEffMTD_4_Sig_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EB_->getNbinsX(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_4_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_4_Sig_EB_, meEle_phi_tot_Sig_EB_, mePhiEffMTD_4_Sig_EB_);

  mePhiEffMTD_5_Sig_EB_ = ibook.book1D("PhiEffMTD_5_Sig_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EB_->getNbinsX(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_5_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_5_Sig_EB_, meEle_phi_tot_Sig_EB_, mePhiEffMTD_5_Sig_EB_);

  mePhiEffMTD_6_Sig_EB_ = ibook.book1D("PhiEffMTD_6_Sig_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EB_->getNbinsX(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_6_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_6_Sig_EB_, meEle_phi_tot_Sig_EB_, mePhiEffMTD_6_Sig_EB_);

  mePhiEffMTD_7_Sig_EB_ = ibook.book1D("PhiEffMTD_7_Sig_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EB_->getNbinsX(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_7_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_7_Sig_EB_, meEle_phi_tot_Sig_EB_, mePhiEffMTD_7_Sig_EB_);





  mePtEffnoMTD_Sig_EB_ = ibook.book1D("pTeffnoMTD_Sig_EB",
                              " noMTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EB_->getNbinsX(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffnoMTD_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_noMTD_Sig_EB_, meEle_pt_tot_Sig_EB_, mePtEffnoMTD_Sig_EB_);

  meEtaEffnoMTD_Sig_EB_ = ibook.book1D("EtaEffnoMTD_Sig_EB",
                              " noMTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EB_->getNbinsX(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffnoMTD_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_noMTD_Sig_EB_, meEle_eta_tot_Sig_EB_, meEtaEffnoMTD_Sig_EB_);

  mePhiEffnoMTD_Sig_EB_ = ibook.book1D("PhiEffnoMTD_Sig_EB",
                              " noMTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EB_->getNbinsX(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffnoMTD_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_noMTD_Sig_EB_, meEle_phi_tot_Sig_EB_, mePhiEffnoMTD_Sig_EB_);


  // Ele iso addition ends
  // For endcap now

  mePtEffMTD_1_Sig_EE_ = ibook.book1D("pTeffMTD_1_Sig_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EE_->getNbinsX(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_1_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_1_Sig_EE_, meEle_pt_tot_Sig_EE_, mePtEffMTD_1_Sig_EE_);

  mePtEffMTD_2_Sig_EE_ = ibook.book1D("pTeffMTD_2_Sig_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EE_->getNbinsX(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_2_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_2_Sig_EE_, meEle_pt_tot_Sig_EE_, mePtEffMTD_2_Sig_EE_);

  mePtEffMTD_3_Sig_EE_ = ibook.book1D("pTeffMTD_3_Sig_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EE_->getNbinsX(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_3_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_3_Sig_EE_, meEle_pt_tot_Sig_EE_, mePtEffMTD_3_Sig_EE_);

  mePtEffMTD_4_Sig_EE_ = ibook.book1D("pTeffMTD_4_Sig_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EE_->getNbinsX(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_4_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_4_Sig_EE_, meEle_pt_tot_Sig_EE_, mePtEffMTD_4_Sig_EE_);

  mePtEffMTD_5_Sig_EE_ = ibook.book1D("pTeffMTD_5_Sig_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EE_->getNbinsX(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_5_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_5_Sig_EE_, meEle_pt_tot_Sig_EE_, mePtEffMTD_5_Sig_EE_);

  mePtEffMTD_6_Sig_EE_ = ibook.book1D("pTeffMTD_6_Sig_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EE_->getNbinsX(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_6_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_6_Sig_EE_, meEle_pt_tot_Sig_EE_, mePtEffMTD_6_Sig_EE_);

  mePtEffMTD_7_Sig_EE_ = ibook.book1D("pTeffMTD_7_Sig_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EE_->getNbinsX(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_7_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_7_Sig_EE_, meEle_pt_tot_Sig_EE_, mePtEffMTD_7_Sig_EE_);






  meEtaEffMTD_1_Sig_EE_ = ibook.book1D("EtaEffMTD_1_Sig_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EE_->getNbinsX(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_1_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_1_Sig_EE_, meEle_eta_tot_Sig_EE_, meEtaEffMTD_1_Sig_EE_);

  meEtaEffMTD_2_Sig_EE_ = ibook.book1D("EtaEffMTD_2_Sig_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EE_->getNbinsX(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_2_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_2_Sig_EE_, meEle_eta_tot_Sig_EE_, meEtaEffMTD_2_Sig_EE_);

  meEtaEffMTD_3_Sig_EE_ = ibook.book1D("EtaEffMTD_3_Sig_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EE_->getNbinsX(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_3_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_3_Sig_EE_, meEle_eta_tot_Sig_EE_, meEtaEffMTD_3_Sig_EE_);

  meEtaEffMTD_4_Sig_EE_ = ibook.book1D("EtaEffMTD_4_Sig_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EE_->getNbinsX(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_4_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_4_Sig_EE_, meEle_eta_tot_Sig_EE_, meEtaEffMTD_4_Sig_EE_);

  meEtaEffMTD_5_Sig_EE_ = ibook.book1D("EtaEffMTD_5_Sig_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EE_->getNbinsX(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_5_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_5_Sig_EE_, meEle_eta_tot_Sig_EE_, meEtaEffMTD_5_Sig_EE_);

  meEtaEffMTD_6_Sig_EE_ = ibook.book1D("EtaEffMTD_6_Sig_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EE_->getNbinsX(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_6_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_6_Sig_EE_, meEle_eta_tot_Sig_EE_, meEtaEffMTD_6_Sig_EE_);

  meEtaEffMTD_7_Sig_EE_ = ibook.book1D("EtaEffMTD_7_Sig_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EE_->getNbinsX(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_7_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_7_Sig_EE_, meEle_eta_tot_Sig_EE_, meEtaEffMTD_7_Sig_EE_);





  mePhiEffMTD_1_Sig_EE_ = ibook.book1D("PhiEffMTD_1_Sig_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EE_->getNbinsX(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_1_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_1_Sig_EE_, meEle_phi_tot_Sig_EE_, mePhiEffMTD_1_Sig_EE_);

  mePhiEffMTD_2_Sig_EE_ = ibook.book1D("PhiEffMTD_2_Sig_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EE_->getNbinsX(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_2_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_2_Sig_EE_, meEle_phi_tot_Sig_EE_, mePhiEffMTD_2_Sig_EE_);

  mePhiEffMTD_3_Sig_EE_ = ibook.book1D("PhiEffMTD_3_Sig_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EE_->getNbinsX(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_3_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_3_Sig_EE_, meEle_phi_tot_Sig_EE_, mePhiEffMTD_3_Sig_EE_);

  mePhiEffMTD_4_Sig_EE_ = ibook.book1D("PhiEffMTD_4_Sig_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EE_->getNbinsX(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_4_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_4_Sig_EE_, meEle_phi_tot_Sig_EE_, mePhiEffMTD_4_Sig_EE_);

  mePhiEffMTD_5_Sig_EE_ = ibook.book1D("PhiEffMTD_5_Sig_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EE_->getNbinsX(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_5_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_5_Sig_EE_, meEle_phi_tot_Sig_EE_, mePhiEffMTD_5_Sig_EE_);

  mePhiEffMTD_6_Sig_EE_ = ibook.book1D("PhiEffMTD_6_Sig_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EE_->getNbinsX(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_6_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_6_Sig_EE_, meEle_phi_tot_Sig_EE_, mePhiEffMTD_6_Sig_EE_);

  mePhiEffMTD_7_Sig_EE_ = ibook.book1D("PhiEffMTD_7_Sig_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EE_->getNbinsX(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_7_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_7_Sig_EE_, meEle_phi_tot_Sig_EE_, mePhiEffMTD_7_Sig_EE_);





  mePtEffnoMTD_Sig_EE_ = ibook.book1D("pTeffnoMTD_Sig_EE",
                              " noMTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Sig_EE_->getNbinsX(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffnoMTD_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_noMTD_Sig_EE_, meEle_pt_tot_Sig_EE_, mePtEffnoMTD_Sig_EE_);

  meEtaEffnoMTD_Sig_EE_ = ibook.book1D("EtaEffnoMTD_Sig_EE",
                              " noMTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Sig_EE_->getNbinsX(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffnoMTD_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_noMTD_Sig_EE_, meEle_eta_tot_Sig_EE_, meEtaEffnoMTD_Sig_EE_);

  mePhiEffnoMTD_Sig_EE_ = ibook.book1D("PhiEffnoMTD_Sig_EE",
                              " noMTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Sig_EE_->getNbinsX(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffnoMTD_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_noMTD_Sig_EE_, meEle_phi_tot_Sig_EE_, mePhiEffnoMTD_Sig_EE_);


  /////////////////////////////////////////////////////// For non-promt (background)

  mePtEffMTD_1_Bkg_EB_ = ibook.book1D("pTeffMTD_1_Bkg_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EB_->getNbinsX(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_1_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_1_Bkg_EB_, meEle_pt_tot_Bkg_EB_, mePtEffMTD_1_Bkg_EB_);

  mePtEffMTD_2_Bkg_EB_ = ibook.book1D("pTeffMTD_2_Bkg_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EB_->getNbinsX(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_2_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_2_Bkg_EB_, meEle_pt_tot_Bkg_EB_, mePtEffMTD_2_Bkg_EB_);

  mePtEffMTD_3_Bkg_EB_ = ibook.book1D("pTeffMTD_3_Bkg_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EB_->getNbinsX(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_3_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_3_Bkg_EB_, meEle_pt_tot_Bkg_EB_, mePtEffMTD_3_Bkg_EB_);

  mePtEffMTD_4_Bkg_EB_ = ibook.book1D("pTeffMTD_4_Bkg_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EB_->getNbinsX(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_4_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_4_Bkg_EB_, meEle_pt_tot_Bkg_EB_, mePtEffMTD_4_Bkg_EB_);

  mePtEffMTD_5_Bkg_EB_ = ibook.book1D("pTeffMTD_5_Bkg_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EB_->getNbinsX(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_5_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_5_Bkg_EB_, meEle_pt_tot_Bkg_EB_, mePtEffMTD_5_Bkg_EB_);

  mePtEffMTD_6_Bkg_EB_ = ibook.book1D("pTeffMTD_6_Bkg_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EB_->getNbinsX(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_6_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_6_Bkg_EB_, meEle_pt_tot_Bkg_EB_, mePtEffMTD_6_Bkg_EB_);

  mePtEffMTD_7_Bkg_EB_ = ibook.book1D("pTeffMTD_7_Bkg_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EB_->getNbinsX(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_7_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_7_Bkg_EB_, meEle_pt_tot_Bkg_EB_, mePtEffMTD_7_Bkg_EB_);






  meEtaEffMTD_1_Bkg_EB_ = ibook.book1D("EtaEffMTD_1_Bkg_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EB_->getNbinsX(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_1_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_1_Bkg_EB_, meEle_eta_tot_Bkg_EB_, meEtaEffMTD_1_Bkg_EB_);

  meEtaEffMTD_2_Bkg_EB_ = ibook.book1D("EtaEffMTD_2_Bkg_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EB_->getNbinsX(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_2_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_2_Bkg_EB_, meEle_eta_tot_Bkg_EB_, meEtaEffMTD_2_Bkg_EB_);

  meEtaEffMTD_3_Bkg_EB_ = ibook.book1D("EtaEffMTD_3_Bkg_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EB_->getNbinsX(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_3_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_3_Bkg_EB_, meEle_eta_tot_Bkg_EB_, meEtaEffMTD_3_Bkg_EB_);

  meEtaEffMTD_4_Bkg_EB_ = ibook.book1D("EtaEffMTD_4_Bkg_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EB_->getNbinsX(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_4_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_4_Bkg_EB_, meEle_eta_tot_Bkg_EB_, meEtaEffMTD_4_Bkg_EB_);

  meEtaEffMTD_5_Bkg_EB_ = ibook.book1D("EtaEffMTD_5_Bkg_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EB_->getNbinsX(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_5_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_5_Bkg_EB_, meEle_eta_tot_Bkg_EB_, meEtaEffMTD_5_Bkg_EB_);

  meEtaEffMTD_6_Bkg_EB_ = ibook.book1D("EtaEffMTD_6_Bkg_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EB_->getNbinsX(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_6_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_6_Bkg_EB_, meEle_eta_tot_Bkg_EB_, meEtaEffMTD_6_Bkg_EB_);

  meEtaEffMTD_7_Bkg_EB_ = ibook.book1D("EtaEffMTD_7_Bkg_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EB_->getNbinsX(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_7_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_7_Bkg_EB_, meEle_eta_tot_Bkg_EB_, meEtaEffMTD_7_Bkg_EB_);





  mePhiEffMTD_1_Bkg_EB_ = ibook.book1D("PhiEffMTD_1_Bkg_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EB_->getNbinsX(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_1_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_1_Bkg_EB_, meEle_phi_tot_Bkg_EB_, mePhiEffMTD_1_Bkg_EB_);

  mePhiEffMTD_2_Bkg_EB_ = ibook.book1D("PhiEffMTD_2_Bkg_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EB_->getNbinsX(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_2_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_2_Bkg_EB_, meEle_phi_tot_Bkg_EB_, mePhiEffMTD_2_Bkg_EB_);

  mePhiEffMTD_3_Bkg_EB_ = ibook.book1D("PhiEffMTD_3_Bkg_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EB_->getNbinsX(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_3_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_3_Bkg_EB_, meEle_phi_tot_Bkg_EB_, mePhiEffMTD_3_Bkg_EB_);

  mePhiEffMTD_4_Bkg_EB_ = ibook.book1D("PhiEffMTD_4_Bkg_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EB_->getNbinsX(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_4_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_4_Bkg_EB_, meEle_phi_tot_Bkg_EB_, mePhiEffMTD_4_Bkg_EB_);

  mePhiEffMTD_5_Bkg_EB_ = ibook.book1D("PhiEffMTD_5_Bkg_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EB_->getNbinsX(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_5_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_5_Bkg_EB_, meEle_phi_tot_Bkg_EB_, mePhiEffMTD_5_Bkg_EB_);

  mePhiEffMTD_6_Bkg_EB_ = ibook.book1D("PhiEffMTD_6_Bkg_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EB_->getNbinsX(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_6_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_6_Bkg_EB_, meEle_phi_tot_Bkg_EB_, mePhiEffMTD_6_Bkg_EB_);

  mePhiEffMTD_7_Bkg_EB_ = ibook.book1D("PhiEffMTD_7_Bkg_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EB_->getNbinsX(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_7_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_7_Bkg_EB_, meEle_phi_tot_Bkg_EB_, mePhiEffMTD_7_Bkg_EB_);





  mePtEffnoMTD_Bkg_EB_ = ibook.book1D("pTeffnoMTD_Bkg_EB",
                              " noMTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EB_->getNbinsX(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffnoMTD_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_noMTD_Bkg_EB_, meEle_pt_tot_Bkg_EB_, mePtEffnoMTD_Bkg_EB_);

  meEtaEffnoMTD_Bkg_EB_ = ibook.book1D("EtaEffnoMTD_Bkg_EB",
                              " noMTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EB_->getNbinsX(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffnoMTD_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_noMTD_Bkg_EB_, meEle_eta_tot_Bkg_EB_, meEtaEffnoMTD_Bkg_EB_);

  mePhiEffnoMTD_Bkg_EB_ = ibook.book1D("PhiEffnoMTD_Bkg_EB",
                              " noMTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EB_->getNbinsX(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffnoMTD_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_noMTD_Bkg_EB_, meEle_phi_tot_Bkg_EB_, mePhiEffnoMTD_Bkg_EB_);


  // Ele iso addition ends
  // For endcap now

  mePtEffMTD_1_Bkg_EE_ = ibook.book1D("pTeffMTD_1_Bkg_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EE_->getNbinsX(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_1_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_1_Bkg_EE_, meEle_pt_tot_Bkg_EE_, mePtEffMTD_1_Bkg_EE_);

  mePtEffMTD_2_Bkg_EE_ = ibook.book1D("pTeffMTD_2_Bkg_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EE_->getNbinsX(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_2_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_2_Bkg_EE_, meEle_pt_tot_Bkg_EE_, mePtEffMTD_2_Bkg_EE_);

  mePtEffMTD_3_Bkg_EE_ = ibook.book1D("pTeffMTD_3_Bkg_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EE_->getNbinsX(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_3_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_3_Bkg_EE_, meEle_pt_tot_Bkg_EE_, mePtEffMTD_3_Bkg_EE_);

  mePtEffMTD_4_Bkg_EE_ = ibook.book1D("pTeffMTD_4_Bkg_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EE_->getNbinsX(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_4_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_4_Bkg_EE_, meEle_pt_tot_Bkg_EE_, mePtEffMTD_4_Bkg_EE_);

  mePtEffMTD_5_Bkg_EE_ = ibook.book1D("pTeffMTD_5_Bkg_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EE_->getNbinsX(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_5_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_5_Bkg_EE_, meEle_pt_tot_Bkg_EE_, mePtEffMTD_5_Bkg_EE_);

  mePtEffMTD_6_Bkg_EE_ = ibook.book1D("pTeffMTD_6_Bkg_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EE_->getNbinsX(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_6_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_6_Bkg_EE_, meEle_pt_tot_Bkg_EE_, mePtEffMTD_6_Bkg_EE_);

  mePtEffMTD_7_Bkg_EE_ = ibook.book1D("pTeffMTD_7_Bkg_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EE_->getNbinsX(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_7_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_7_Bkg_EE_, meEle_pt_tot_Bkg_EE_, mePtEffMTD_7_Bkg_EE_);






  meEtaEffMTD_1_Bkg_EE_ = ibook.book1D("EtaEffMTD_1_Bkg_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EE_->getNbinsX(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_1_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_1_Bkg_EE_, meEle_eta_tot_Bkg_EE_, meEtaEffMTD_1_Bkg_EE_);

  meEtaEffMTD_2_Bkg_EE_ = ibook.book1D("EtaEffMTD_2_Bkg_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EE_->getNbinsX(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_2_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_2_Bkg_EE_, meEle_eta_tot_Bkg_EE_, meEtaEffMTD_2_Bkg_EE_);

  meEtaEffMTD_3_Bkg_EE_ = ibook.book1D("EtaEffMTD_3_Bkg_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EE_->getNbinsX(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_3_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_3_Bkg_EE_, meEle_eta_tot_Bkg_EE_, meEtaEffMTD_3_Bkg_EE_);

  meEtaEffMTD_4_Bkg_EE_ = ibook.book1D("EtaEffMTD_4_Bkg_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EE_->getNbinsX(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_4_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_4_Bkg_EE_, meEle_eta_tot_Bkg_EE_, meEtaEffMTD_4_Bkg_EE_);

  meEtaEffMTD_5_Bkg_EE_ = ibook.book1D("EtaEffMTD_5_Bkg_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EE_->getNbinsX(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_5_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_5_Bkg_EE_, meEle_eta_tot_Bkg_EE_, meEtaEffMTD_5_Bkg_EE_);

  meEtaEffMTD_6_Bkg_EE_ = ibook.book1D("EtaEffMTD_6_Bkg_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EE_->getNbinsX(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_6_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_6_Bkg_EE_, meEle_eta_tot_Bkg_EE_, meEtaEffMTD_6_Bkg_EE_);

  meEtaEffMTD_7_Bkg_EE_ = ibook.book1D("EtaEffMTD_7_Bkg_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EE_->getNbinsX(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_7_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_7_Bkg_EE_, meEle_eta_tot_Bkg_EE_, meEtaEffMTD_7_Bkg_EE_);





  mePhiEffMTD_1_Bkg_EE_ = ibook.book1D("PhiEffMTD_1_Bkg_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EE_->getNbinsX(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_1_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_1_Bkg_EE_, meEle_phi_tot_Bkg_EE_, mePhiEffMTD_1_Bkg_EE_);

  mePhiEffMTD_2_Bkg_EE_ = ibook.book1D("PhiEffMTD_2_Bkg_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EE_->getNbinsX(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_2_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_2_Bkg_EE_, meEle_phi_tot_Bkg_EE_, mePhiEffMTD_2_Bkg_EE_);

  mePhiEffMTD_3_Bkg_EE_ = ibook.book1D("PhiEffMTD_3_Bkg_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EE_->getNbinsX(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_3_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_3_Bkg_EE_, meEle_phi_tot_Bkg_EE_, mePhiEffMTD_3_Bkg_EE_);

  mePhiEffMTD_4_Bkg_EE_ = ibook.book1D("PhiEffMTD_4_Bkg_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EE_->getNbinsX(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_4_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_4_Bkg_EE_, meEle_phi_tot_Bkg_EE_, mePhiEffMTD_4_Bkg_EE_);

  mePhiEffMTD_5_Bkg_EE_ = ibook.book1D("PhiEffMTD_5_Bkg_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EE_->getNbinsX(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_5_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_5_Bkg_EE_, meEle_phi_tot_Bkg_EE_, mePhiEffMTD_5_Bkg_EE_);

  mePhiEffMTD_6_Bkg_EE_ = ibook.book1D("PhiEffMTD_6_Bkg_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EE_->getNbinsX(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_6_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_6_Bkg_EE_, meEle_phi_tot_Bkg_EE_, mePhiEffMTD_6_Bkg_EE_);

  mePhiEffMTD_7_Bkg_EE_ = ibook.book1D("PhiEffMTD_7_Bkg_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EE_->getNbinsX(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_7_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_7_Bkg_EE_, meEle_phi_tot_Bkg_EE_, mePhiEffMTD_7_Bkg_EE_);





  mePtEffnoMTD_Bkg_EE_ = ibook.book1D("pTeffnoMTD_Bkg_EE",
                              " noMTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_Bkg_EE_->getNbinsX(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffnoMTD_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_noMTD_Bkg_EE_, meEle_pt_tot_Bkg_EE_, mePtEffnoMTD_Bkg_EE_);

  meEtaEffnoMTD_Bkg_EE_ = ibook.book1D("EtaEffnoMTD_Bkg_EE",
                              " noMTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_Bkg_EE_->getNbinsX(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffnoMTD_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_noMTD_Bkg_EE_, meEle_eta_tot_Bkg_EE_, meEtaEffnoMTD_Bkg_EE_);

  mePhiEffnoMTD_Bkg_EE_ = ibook.book1D("PhiEffnoMTD_Bkg_EE",
                              " noMTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_Bkg_EE_->getNbinsX(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffnoMTD_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_noMTD_Bkg_EE_, meEle_phi_tot_Bkg_EE_, mePhiEffnoMTD_Bkg_EE_);

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ----------
void MtdEleIsoHarvester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/ElectronIso/");

  descriptions.add("MtdEleIsoPostProcessor", desc);
}

DEFINE_FWK_MODULE(MtdEleIsoHarvester);


