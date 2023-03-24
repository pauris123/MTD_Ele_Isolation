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
  MonitorElement* mePtEffnoMTD_EB_;
  MonitorElement* meEtaEffnoMTD_EB_;
  MonitorElement* mePhiEffnoMTD_EB_;
  
  MonitorElement* mePtEffMTD_1_EB_;
  MonitorElement* meEtaEffMTD_1_EB_;
  MonitorElement* mePhiEffMTD_1_EB_;

  MonitorElement* mePtEffMTD_2_EB_;
  MonitorElement* meEtaEffMTD_2_EB_;
  MonitorElement* mePhiEffMTD_2_EB_;

  MonitorElement* mePtEffMTD_3_EB_;
  MonitorElement* meEtaEffMTD_3_EB_;
  MonitorElement* mePhiEffMTD_3_EB_;

  MonitorElement* mePtEffMTD_4_EB_;
  MonitorElement* meEtaEffMTD_4_EB_;
  MonitorElement* mePhiEffMTD_4_EB_;

  MonitorElement* mePtEffMTD_5_EB_;
  MonitorElement* meEtaEffMTD_5_EB_;
  MonitorElement* mePhiEffMTD_5_EB_;

  MonitorElement* mePtEffMTD_6_EB_;
  MonitorElement* meEtaEffMTD_6_EB_;
  MonitorElement* mePhiEffMTD_6_EB_;

  MonitorElement* mePtEffMTD_7_EB_;
  MonitorElement* meEtaEffMTD_7_EB_;
  MonitorElement* mePhiEffMTD_7_EB_;


  MonitorElement* mePtEffnoMTD_EE_;
  MonitorElement* meEtaEffnoMTD_EE_;
  MonitorElement* mePhiEffnoMTD_EE_;
  
  MonitorElement* mePtEffMTD_1_EE_;
  MonitorElement* meEtaEffMTD_1_EE_;
  MonitorElement* mePhiEffMTD_1_EE_;

  MonitorElement* mePtEffMTD_2_EE_;
  MonitorElement* meEtaEffMTD_2_EE_;
  MonitorElement* mePhiEffMTD_2_EE_;

  MonitorElement* mePtEffMTD_3_EE_;
  MonitorElement* meEtaEffMTD_3_EE_;
  MonitorElement* mePhiEffMTD_3_EE_;

  MonitorElement* mePtEffMTD_4_EE_;
  MonitorElement* meEtaEffMTD_4_EE_;
  MonitorElement* mePhiEffMTD_4_EE_;

  MonitorElement* mePtEffMTD_5_EE_;
  MonitorElement* meEtaEffMTD_5_EE_;
  MonitorElement* mePhiEffMTD_5_EE_;

  MonitorElement* mePtEffMTD_6_EE_;
  MonitorElement* meEtaEffMTD_6_EE_;
  MonitorElement* mePhiEffMTD_6_EE_;

  MonitorElement* mePtEffMTD_7_EE_;
  MonitorElement* meEtaEffMTD_7_EE_;
  MonitorElement* mePhiEffMTD_7_EE_;


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

  MonitorElement* meEle_pt_tot_EB_ = igetter.get(folder_ + "Ele_pT_tot_EB");
  MonitorElement* meEle_pt_noMTD_EB_ = igetter.get(folder_ + "Ele_pT_noMTD_EB");

  MonitorElement* meEle_eta_tot_EB_ = igetter.get(folder_ + "Ele_eta_tot_EB");
  MonitorElement* meEle_eta_noMTD_EB_ = igetter.get(folder_ + "Ele_eta_noMTD_EB");

  MonitorElement* meEle_phi_tot_EB_ = igetter.get(folder_ + "Ele_phi_tot_EB");
  MonitorElement* meEle_phi_noMTD_EB_ = igetter.get(folder_ + "Ele_phi_noMTD_EB");

  MonitorElement* meEle_pt_MTD_1_EB_ = igetter.get(folder_ + "Ele_pT_MTD_1_EB");
  MonitorElement* meEle_eta_MTD_1_EB_ = igetter.get(folder_ + "Ele_eta_MTD_1_EB");
  MonitorElement* meEle_phi_MTD_1_EB_ = igetter.get(folder_ + "Ele_phi_MTD_1_EB");

  MonitorElement* meEle_pt_MTD_2_EB_ = igetter.get(folder_ + "Ele_pT_MTD_2_EB");
  MonitorElement* meEle_eta_MTD_2_EB_ = igetter.get(folder_ + "Ele_eta_MTD_2_EB");
  MonitorElement* meEle_phi_MTD_2_EB_ = igetter.get(folder_ + "Ele_phi_MTD_2_EB");

  MonitorElement* meEle_pt_MTD_3_EB_ = igetter.get(folder_ + "Ele_pT_MTD_3_EB");
  MonitorElement* meEle_eta_MTD_3_EB_ = igetter.get(folder_ + "Ele_eta_MTD_3_EB");
  MonitorElement* meEle_phi_MTD_3_EB_ = igetter.get(folder_ + "Ele_phi_MTD_3_EB");

  MonitorElement* meEle_pt_MTD_4_EB_ = igetter.get(folder_ + "Ele_pT_MTD_4_EB");
  MonitorElement* meEle_eta_MTD_4_EB_ = igetter.get(folder_ + "Ele_eta_MTD_4_EB");
  MonitorElement* meEle_phi_MTD_4_EB_ = igetter.get(folder_ + "Ele_phi_MTD_4_EB");

  MonitorElement* meEle_pt_MTD_5_EB_ = igetter.get(folder_ + "Ele_pT_MTD_5_EB");
  MonitorElement* meEle_eta_MTD_5_EB_ = igetter.get(folder_ + "Ele_eta_MTD_5_EB");
  MonitorElement* meEle_phi_MTD_5_EB_ = igetter.get(folder_ + "Ele_phi_MTD_5_EB");

  MonitorElement* meEle_pt_MTD_6_EB_ = igetter.get(folder_ + "Ele_pT_MTD_6_EB");
  MonitorElement* meEle_eta_MTD_6_EB_ = igetter.get(folder_ + "Ele_eta_MTD_6_EB");
  MonitorElement* meEle_phi_MTD_6_EB_ = igetter.get(folder_ + "Ele_phi_MTD_6_EB");

  MonitorElement* meEle_pt_MTD_7_EB_ = igetter.get(folder_ + "Ele_pT_MTD_7_EB");
  MonitorElement* meEle_eta_MTD_7_EB_ = igetter.get(folder_ + "Ele_eta_MTD_7_EB");
  MonitorElement* meEle_phi_MTD_7_EB_ = igetter.get(folder_ + "Ele_phi_MTD_7_EB");


  MonitorElement* meEle_pt_tot_EE_ = igetter.get(folder_ + "Ele_pT_tot_EE");
  MonitorElement* meEle_pt_noMTD_EE_ = igetter.get(folder_ + "Ele_pT_noMTD_EE");

  MonitorElement* meEle_eta_tot_EE_ = igetter.get(folder_ + "Ele_eta_tot_EE");
  MonitorElement* meEle_eta_noMTD_EE_ = igetter.get(folder_ + "Ele_eta_noMTD_EE");

  MonitorElement* meEle_phi_tot_EE_ = igetter.get(folder_ + "Ele_phi_tot_EE");
  MonitorElement* meEle_phi_noMTD_EE_ = igetter.get(folder_ + "Ele_phi_noMTD_EE");

  MonitorElement* meEle_pt_MTD_1_EE_ = igetter.get(folder_ + "Ele_pT_MTD_1_EE");
  MonitorElement* meEle_eta_MTD_1_EE_ = igetter.get(folder_ + "Ele_eta_MTD_1_EE");
  MonitorElement* meEle_phi_MTD_1_EE_ = igetter.get(folder_ + "Ele_phi_MTD_1_EE");

  MonitorElement* meEle_pt_MTD_2_EE_ = igetter.get(folder_ + "Ele_pT_MTD_2_EE");
  MonitorElement* meEle_eta_MTD_2_EE_ = igetter.get(folder_ + "Ele_eta_MTD_2_EE");
  MonitorElement* meEle_phi_MTD_2_EE_ = igetter.get(folder_ + "Ele_phi_MTD_2_EE");

  MonitorElement* meEle_pt_MTD_3_EE_ = igetter.get(folder_ + "Ele_pT_MTD_3_EE");
  MonitorElement* meEle_eta_MTD_3_EE_ = igetter.get(folder_ + "Ele_eta_MTD_3_EE");
  MonitorElement* meEle_phi_MTD_3_EE_ = igetter.get(folder_ + "Ele_phi_MTD_3_EE");

  MonitorElement* meEle_pt_MTD_4_EE_ = igetter.get(folder_ + "Ele_pT_MTD_4_EE");
  MonitorElement* meEle_eta_MTD_4_EE_ = igetter.get(folder_ + "Ele_eta_MTD_4_EE");
  MonitorElement* meEle_phi_MTD_4_EE_ = igetter.get(folder_ + "Ele_phi_MTD_4_EE");

  MonitorElement* meEle_pt_MTD_5_EE_ = igetter.get(folder_ + "Ele_pT_MTD_5_EE");
  MonitorElement* meEle_eta_MTD_5_EE_ = igetter.get(folder_ + "Ele_eta_MTD_5_EE");
  MonitorElement* meEle_phi_MTD_5_EE_ = igetter.get(folder_ + "Ele_phi_MTD_5_EE");

  MonitorElement* meEle_pt_MTD_6_EE_ = igetter.get(folder_ + "Ele_pT_MTD_6_EE");
  MonitorElement* meEle_eta_MTD_6_EE_ = igetter.get(folder_ + "Ele_eta_MTD_6_EE");
  MonitorElement* meEle_phi_MTD_6_EE_ = igetter.get(folder_ + "Ele_phi_MTD_6_EE");

  MonitorElement* meEle_pt_MTD_7_EE_ = igetter.get(folder_ + "Ele_pT_MTD_7_EE");
  MonitorElement* meEle_eta_MTD_7_EE_ = igetter.get(folder_ + "Ele_eta_MTD_7_EE");
  MonitorElement* meEle_phi_MTD_7_EE_ = igetter.get(folder_ + "Ele_phi_MTD_7_EE");


  if (!meEle_pt_tot_EB_ || !meEle_pt_MTD_1_EB_ || !meEle_pt_MTD_2_EB_ ||
      !meEle_pt_MTD_3_EB_ || !meEle_pt_MTD_4_EB_ || !meEle_pt_MTD_5_EB_ || !meEle_pt_MTD_6_EB_ || !meEle_pt_MTD_7_EB_ ||
      !meEle_pt_noMTD_EB_ || !meEle_eta_tot_EB_ || !meEle_eta_MTD_1_EB_ || !meEle_eta_MTD_2_EB_ || !meEle_eta_MTD_3_EB_ ||
      !meEle_eta_MTD_4_EB_ || !meEle_eta_MTD_5_EB_ || !meEle_eta_MTD_6_EB_ || !meEle_eta_MTD_7_EB_ || !meEle_eta_noMTD_EB_ || 
      !meEle_phi_tot_EB_ || !meEle_phi_MTD_1_EB_ || !meEle_phi_MTD_2_EB_ || !meEle_phi_MTD_3_EB_ || !meEle_phi_MTD_4_EB_ ||
      !meEle_phi_MTD_5_EB_ || !meEle_phi_MTD_6_EB_ || !meEle_phi_MTD_7_EB_ || !meEle_phi_noMTD_EB_ || 
      !meEle_pt_tot_EE_ || !meEle_pt_MTD_1_EE_ || !meEle_pt_MTD_2_EE_ ||
      !meEle_pt_MTD_3_EE_ || !meEle_pt_MTD_4_EE_ || !meEle_pt_MTD_5_EE_ || !meEle_pt_MTD_6_EE_ || !meEle_pt_MTD_7_EE_ ||
      !meEle_pt_noMTD_EE_ || !meEle_eta_tot_EE_ || !meEle_eta_MTD_1_EE_ || !meEle_eta_MTD_2_EE_ || !meEle_eta_MTD_3_EE_ ||
      !meEle_eta_MTD_4_EE_ || !meEle_eta_MTD_5_EE_ || !meEle_eta_MTD_6_EE_ || !meEle_eta_MTD_7_EE_ || !meEle_eta_noMTD_EE_ || 
      !meEle_phi_tot_EE_ || !meEle_phi_MTD_1_EE_ || !meEle_phi_MTD_2_EE_ || !meEle_phi_MTD_3_EE_ || !meEle_phi_MTD_4_EE_ ||
      !meEle_phi_MTD_5_EE_ || !meEle_phi_MTD_6_EE_ || !meEle_phi_MTD_7_EE_ || !meEle_phi_noMTD_EE_) {
    edm::LogError("MtdEleIsoHarvester") << "Monitoring histograms not found!" << std::endl;
    return;
  }

  // --- Book  histograms
  ibook.cd(folder_);
  // ele iso addition starts; MTD vs noMTD case

  mePtEffMTD_1_EB_ = ibook.book1D("pTeffMTD_1_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EB_->getNbinsX(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_1_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_1_EB_, meEle_pt_tot_EB_, mePtEffMTD_1_EB_);

  mePtEffMTD_2_EB_ = ibook.book1D("pTeffMTD_2_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EB_->getNbinsX(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_2_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_2_EB_, meEle_pt_tot_EB_, mePtEffMTD_2_EB_);

  mePtEffMTD_3_EB_ = ibook.book1D("pTeffMTD_3_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EB_->getNbinsX(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_3_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_3_EB_, meEle_pt_tot_EB_, mePtEffMTD_3_EB_);

  mePtEffMTD_4_EB_ = ibook.book1D("pTeffMTD_4_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EB_->getNbinsX(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_4_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_4_EB_, meEle_pt_tot_EB_, mePtEffMTD_4_EB_);

  mePtEffMTD_5_EB_ = ibook.book1D("pTeffMTD_5_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EB_->getNbinsX(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_5_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_5_EB_, meEle_pt_tot_EB_, mePtEffMTD_5_EB_);

  mePtEffMTD_6_EB_ = ibook.book1D("pTeffMTD_6_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EB_->getNbinsX(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_6_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_6_EB_, meEle_pt_tot_EB_, mePtEffMTD_6_EB_);

  mePtEffMTD_7_EB_ = ibook.book1D("pTeffMTD_7_EB",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EB_->getNbinsX(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_7_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_7_EB_, meEle_pt_tot_EB_, mePtEffMTD_7_EB_);






  meEtaEffMTD_1_EB_ = ibook.book1D("EtaEffMTD_1_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EB_->getNbinsX(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_1_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_1_EB_, meEle_eta_tot_EB_, meEtaEffMTD_1_EB_);

  meEtaEffMTD_2_EB_ = ibook.book1D("EtaEffMTD_2_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EB_->getNbinsX(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_2_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_2_EB_, meEle_eta_tot_EB_, meEtaEffMTD_2_EB_);

  meEtaEffMTD_3_EB_ = ibook.book1D("EtaEffMTD_3_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EB_->getNbinsX(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_3_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_3_EB_, meEle_eta_tot_EB_, meEtaEffMTD_3_EB_);

  meEtaEffMTD_4_EB_ = ibook.book1D("EtaEffMTD_4_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EB_->getNbinsX(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_4_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_4_EB_, meEle_eta_tot_EB_, meEtaEffMTD_4_EB_);

  meEtaEffMTD_5_EB_ = ibook.book1D("EtaEffMTD_5_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EB_->getNbinsX(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_5_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_5_EB_, meEle_eta_tot_EB_, meEtaEffMTD_5_EB_);

  meEtaEffMTD_6_EB_ = ibook.book1D("EtaEffMTD_6_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EB_->getNbinsX(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_6_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_6_EB_, meEle_eta_tot_EB_, meEtaEffMTD_6_EB_);

  meEtaEffMTD_7_EB_ = ibook.book1D("EtaEffMTD_7_EB",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EB_->getNbinsX(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_7_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_7_EB_, meEle_eta_tot_EB_, meEtaEffMTD_7_EB_);





  mePhiEffMTD_1_EB_ = ibook.book1D("PhiEffMTD_1_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EB_->getNbinsX(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_1_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_1_EB_, meEle_phi_tot_EB_, mePhiEffMTD_1_EB_);

  mePhiEffMTD_2_EB_ = ibook.book1D("PhiEffMTD_2_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EB_->getNbinsX(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_2_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_2_EB_, meEle_phi_tot_EB_, mePhiEffMTD_2_EB_);

  mePhiEffMTD_3_EB_ = ibook.book1D("PhiEffMTD_3_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EB_->getNbinsX(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_3_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_3_EB_, meEle_phi_tot_EB_, mePhiEffMTD_3_EB_);

  mePhiEffMTD_4_EB_ = ibook.book1D("PhiEffMTD_4_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EB_->getNbinsX(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_4_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_4_EB_, meEle_phi_tot_EB_, mePhiEffMTD_4_EB_);

  mePhiEffMTD_5_EB_ = ibook.book1D("PhiEffMTD_5_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EB_->getNbinsX(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_5_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_5_EB_, meEle_phi_tot_EB_, mePhiEffMTD_5_EB_);

  mePhiEffMTD_6_EB_ = ibook.book1D("PhiEffMTD_6_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EB_->getNbinsX(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_6_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_6_EB_, meEle_phi_tot_EB_, mePhiEffMTD_6_EB_);

  mePhiEffMTD_7_EB_ = ibook.book1D("PhiEffMTD_7_EB",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EB_->getNbinsX(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_7_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_7_EB_, meEle_phi_tot_EB_, mePhiEffMTD_7_EB_);





  mePtEffnoMTD_EB_ = ibook.book1D("pTeffnoMTD_EB",
                              " noMTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EB_->getNbinsX(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffnoMTD_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_noMTD_EB_, meEle_pt_tot_EB_, mePtEffnoMTD_EB_);

  meEtaEffnoMTD_EB_ = ibook.book1D("EtaEffnoMTD_EB",
                              " noMTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EB_->getNbinsX(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffnoMTD_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_noMTD_EB_, meEle_eta_tot_EB_, meEtaEffnoMTD_EB_);

  mePhiEffnoMTD_EB_ = ibook.book1D("PhiEffnoMTD_EB",
                              " noMTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EB_->getNbinsX(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffnoMTD_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_noMTD_EB_, meEle_phi_tot_EB_, mePhiEffnoMTD_EB_);


  // Ele iso addition ends
  // For endcap now

  mePtEffMTD_1_EE_ = ibook.book1D("pTeffMTD_1_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EE_->getNbinsX(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_1_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_1_EE_, meEle_pt_tot_EE_, mePtEffMTD_1_EE_);

  mePtEffMTD_2_EE_ = ibook.book1D("pTeffMTD_2_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EE_->getNbinsX(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_2_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_2_EE_, meEle_pt_tot_EE_, mePtEffMTD_2_EE_);

  mePtEffMTD_3_EE_ = ibook.book1D("pTeffMTD_3_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EE_->getNbinsX(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_3_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_3_EE_, meEle_pt_tot_EE_, mePtEffMTD_3_EE_);

  mePtEffMTD_4_EE_ = ibook.book1D("pTeffMTD_4_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EE_->getNbinsX(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_4_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_4_EE_, meEle_pt_tot_EE_, mePtEffMTD_4_EE_);

  mePtEffMTD_5_EE_ = ibook.book1D("pTeffMTD_5_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EE_->getNbinsX(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_5_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_5_EE_, meEle_pt_tot_EE_, mePtEffMTD_5_EE_);

  mePtEffMTD_6_EE_ = ibook.book1D("pTeffMTD_6_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EE_->getNbinsX(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_6_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_6_EE_, meEle_pt_tot_EE_, mePtEffMTD_6_EE_);

  mePtEffMTD_7_EE_ = ibook.book1D("pTeffMTD_7_EE",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EE_->getNbinsX(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_7_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_7_EE_, meEle_pt_tot_EE_, mePtEffMTD_7_EE_);






  meEtaEffMTD_1_EE_ = ibook.book1D("EtaEffMTD_1_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EE_->getNbinsX(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_1_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_1_EE_, meEle_eta_tot_EE_, meEtaEffMTD_1_EE_);

  meEtaEffMTD_2_EE_ = ibook.book1D("EtaEffMTD_2_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EE_->getNbinsX(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_2_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_2_EE_, meEle_eta_tot_EE_, meEtaEffMTD_2_EE_);

  meEtaEffMTD_3_EE_ = ibook.book1D("EtaEffMTD_3_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EE_->getNbinsX(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_3_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_3_EE_, meEle_eta_tot_EE_, meEtaEffMTD_3_EE_);

  meEtaEffMTD_4_EE_ = ibook.book1D("EtaEffMTD_4_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EE_->getNbinsX(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_4_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_4_EE_, meEle_eta_tot_EE_, meEtaEffMTD_4_EE_);

  meEtaEffMTD_5_EE_ = ibook.book1D("EtaEffMTD_5_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EE_->getNbinsX(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_5_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_5_EE_, meEle_eta_tot_EE_, meEtaEffMTD_5_EE_);

  meEtaEffMTD_6_EE_ = ibook.book1D("EtaEffMTD_6_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EE_->getNbinsX(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_6_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_6_EE_, meEle_eta_tot_EE_, meEtaEffMTD_6_EE_);

  meEtaEffMTD_7_EE_ = ibook.book1D("EtaEffMTD_7_EE",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EE_->getNbinsX(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_7_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_7_EE_, meEle_eta_tot_EE_, meEtaEffMTD_7_EE_);





  mePhiEffMTD_1_EE_ = ibook.book1D("PhiEffMTD_1_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EE_->getNbinsX(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_1_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_1_EE_, meEle_phi_tot_EE_, mePhiEffMTD_1_EE_);

  mePhiEffMTD_2_EE_ = ibook.book1D("PhiEffMTD_2_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EE_->getNbinsX(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_2_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_2_EE_, meEle_phi_tot_EE_, mePhiEffMTD_2_EE_);

  mePhiEffMTD_3_EE_ = ibook.book1D("PhiEffMTD_3_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EE_->getNbinsX(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_3_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_3_EE_, meEle_phi_tot_EE_, mePhiEffMTD_3_EE_);

  mePhiEffMTD_4_EE_ = ibook.book1D("PhiEffMTD_4_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EE_->getNbinsX(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_4_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_4_EE_, meEle_phi_tot_EE_, mePhiEffMTD_4_EE_);

  mePhiEffMTD_5_EE_ = ibook.book1D("PhiEffMTD_5_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EE_->getNbinsX(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_5_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_5_EE_, meEle_phi_tot_EE_, mePhiEffMTD_5_EE_);

  mePhiEffMTD_6_EE_ = ibook.book1D("PhiEffMTD_6_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EE_->getNbinsX(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_6_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_6_EE_, meEle_phi_tot_EE_, mePhiEffMTD_6_EE_);

  mePhiEffMTD_7_EE_ = ibook.book1D("PhiEffMTD_7_EE",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EE_->getNbinsX(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_7_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_7_EE_, meEle_phi_tot_EE_, mePhiEffMTD_7_EE_);





  mePtEffnoMTD_EE_ = ibook.book1D("pTeffnoMTD_EE",
                              " noMTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot_EE_->getNbinsX(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffnoMTD_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_noMTD_EE_, meEle_pt_tot_EE_, mePtEffnoMTD_EE_);

  meEtaEffnoMTD_EE_ = ibook.book1D("EtaEffnoMTD_EE",
                              " noMTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot_EE_->getNbinsX(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffnoMTD_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_noMTD_EE_, meEle_eta_tot_EE_, meEtaEffnoMTD_EE_);

  mePhiEffnoMTD_EE_ = ibook.book1D("PhiEffnoMTD_EE",
                              " noMTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot_EE_->getNbinsX(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffnoMTD_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_noMTD_EE_, meEle_phi_tot_EE_, mePhiEffnoMTD_EE_);

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ----------
void MtdEleIsoHarvester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/Tracks/");

  descriptions.add("MtdEleIsoPostProcessor", desc);
}

DEFINE_FWK_MODULE(MtdEleIsoHarvester);


