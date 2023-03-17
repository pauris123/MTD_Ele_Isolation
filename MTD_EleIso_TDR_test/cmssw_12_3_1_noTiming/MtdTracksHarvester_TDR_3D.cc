#include <string>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

class MtdTracksHarvester : public DQMEDHarvester {
public:
  explicit MtdTracksHarvester(const edm::ParameterSet& iConfig);
  ~MtdTracksHarvester() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

protected:
  void dqmEndJob(DQMStore::IBooker&, DQMStore::IGetter&) override;

private:
  void computeEfficiency1D(MonitorElement* num, MonitorElement* den, MonitorElement* result);

  const std::string folder_;

  // --- Histograms
  MonitorElement* meBtlEtaEff_;
  MonitorElement* meBtlPhiEff_;
  MonitorElement* meBtlPtEff_;
  MonitorElement* meEtlEtaEff_[2];
  MonitorElement* meEtlPhiEff_[2];
  MonitorElement* meEtlPtEff_[2];
  
  MonitorElement* meMVAPtSelEff_;
  MonitorElement* meMVAEtaSelEff_;
  MonitorElement* meMVAPtMatchEff_;
  MonitorElement* meMVAEtaMatchEff_;

  MonitorElement* mePtEffnoMTD_;
  MonitorElement* meEtaEffnoMTD_;
  MonitorElement* mePhiEffnoMTD_;
  
  MonitorElement* mePtEffMTD_1_;
  MonitorElement* meEtaEffMTD_1_;
  MonitorElement* mePhiEffMTD_1_;

  MonitorElement* mePtEffMTD_2_;
  MonitorElement* meEtaEffMTD_2_;
  MonitorElement* mePhiEffMTD_2_;

  MonitorElement* mePtEffMTD_3_;
  MonitorElement* meEtaEffMTD_3_;
  MonitorElement* mePhiEffMTD_3_;

  MonitorElement* mePtEffMTD_4_;
  MonitorElement* meEtaEffMTD_4_;
  MonitorElement* mePhiEffMTD_4_;

  MonitorElement* mePtEffMTD_5_;
  MonitorElement* meEtaEffMTD_5_;
  MonitorElement* mePhiEffMTD_5_;

  MonitorElement* mePtEffMTD_6_;
  MonitorElement* meEtaEffMTD_6_;
  MonitorElement* mePhiEffMTD_6_;

  MonitorElement* mePtEffMTD_7_;
  MonitorElement* meEtaEffMTD_7_;
  MonitorElement* mePhiEffMTD_7_;

  
};

// ------------ constructor and destructor --------------
MtdTracksHarvester::MtdTracksHarvester(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")) {}

MtdTracksHarvester::~MtdTracksHarvester() {}

// auxiliary method to compute efficiency from the ratio of two 1D MonitorElement
void MtdTracksHarvester::computeEfficiency1D(MonitorElement* num, MonitorElement* den, MonitorElement* result) {
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
void MtdTracksHarvester::dqmEndJob(DQMStore::IBooker& ibook, DQMStore::IGetter& igetter) {
  // --- Get the monitoring histograms
  MonitorElement* meBTLTrackEffEtaTot = igetter.get(folder_ + "TrackBTLEffEtaTot");
  MonitorElement* meBTLTrackEffPhiTot = igetter.get(folder_ + "TrackBTLEffPhiTot");
  MonitorElement* meBTLTrackEffPtTot = igetter.get(folder_ + "TrackBTLEffPtTot");
  MonitorElement* meBTLTrackEffEtaMtd = igetter.get(folder_ + "TrackBTLEffEtaMtd");
  MonitorElement* meBTLTrackEffPhiMtd = igetter.get(folder_ + "TrackBTLEffPhiMtd");
  MonitorElement* meBTLTrackEffPtMtd = igetter.get(folder_ + "TrackBTLEffPtMtd");
  MonitorElement* meETLTrackEffEtaTotZneg = igetter.get(folder_ + "TrackETLEffEtaTotZneg");
  MonitorElement* meETLTrackEffPhiTotZneg = igetter.get(folder_ + "TrackETLEffPhiTotZneg");
  MonitorElement* meETLTrackEffPtTotZneg = igetter.get(folder_ + "TrackETLEffPtTotZneg");
  MonitorElement* meETLTrackEffEtaMtdZneg = igetter.get(folder_ + "TrackETLEffEtaMtdZneg");
  MonitorElement* meETLTrackEffPhiMtdZneg = igetter.get(folder_ + "TrackETLEffPhiMtdZneg");
  MonitorElement* meETLTrackEffPtMtdZneg = igetter.get(folder_ + "TrackETLEffPtMtdZneg");
  
  MonitorElement* meETLTrackEffEtaTotZpos = igetter.get(folder_ + "TrackETLEffEtaTotZpos");
  MonitorElement* meETLTrackEffPhiTotZpos = igetter.get(folder_ + "TrackETLEffPhiTotZpos");
  MonitorElement* meETLTrackEffPtTotZpos = igetter.get(folder_ + "TrackETLEffPtTotZpos");
  MonitorElement* meETLTrackEffEtaMtdZpos = igetter.get(folder_ + "TrackETLEffEtaMtdZpos");
  MonitorElement* meETLTrackEffPhiMtdZpos = igetter.get(folder_ + "TrackETLEffPhiMtdZpos");
  MonitorElement* meETLTrackEffPtMtdZpos = igetter.get(folder_ + "TrackETLEffPtMtdZpos");
  
  MonitorElement* meMVATrackEffPtTot = igetter.get(folder_ + "MVAEffPtTot");
  MonitorElement* meMVATrackMatchedEffPtTot = igetter.get(folder_ + "MVAMatchedEffPtTot");
  MonitorElement* meMVATrackMatchedEffPtMtd = igetter.get(folder_ + "MVAMatchedEffPtMtd");
  
  MonitorElement* meMVATrackEffEtaTot = igetter.get(folder_ + "MVAEffEtaTot");
  MonitorElement* meMVATrackMatchedEffEtaTot = igetter.get(folder_ + "MVAMatchedEffEtaTot");
  MonitorElement* meMVATrackMatchedEffEtaMtd = igetter.get(folder_ + "MVAMatchedEffEtaMtd");

  MonitorElement* meEle_pt_tot = igetter.get(folder_ + "Ele_pT_tot");
  MonitorElement* meEle_pt_noMTD = igetter.get(folder_ + "Ele_pT_noMTD");

  MonitorElement* meEle_eta_tot = igetter.get(folder_ + "Ele_eta_tot");
  MonitorElement* meEle_eta_noMTD = igetter.get(folder_ + "Ele_eta_noMTD");

  MonitorElement* meEle_phi_tot = igetter.get(folder_ + "Ele_phi_tot");
  MonitorElement* meEle_phi_noMTD = igetter.get(folder_ + "Ele_phi_noMTD");

  MonitorElement* meEle_pt_MTD_1_ = igetter.get(folder_ + "Ele_pT_MTD_1");
  MonitorElement* meEle_eta_MTD_1_ = igetter.get(folder_ + "Ele_eta_MTD_1");
  MonitorElement* meEle_phi_MTD_1_ = igetter.get(folder_ + "Ele_phi_MTD_1");

  MonitorElement* meEle_pt_MTD_2_ = igetter.get(folder_ + "Ele_pT_MTD_2");
  MonitorElement* meEle_eta_MTD_2_ = igetter.get(folder_ + "Ele_eta_MTD_2");
  MonitorElement* meEle_phi_MTD_2_ = igetter.get(folder_ + "Ele_phi_MTD_2");

  MonitorElement* meEle_pt_MTD_3_ = igetter.get(folder_ + "Ele_pT_MTD_3");
  MonitorElement* meEle_eta_MTD_3_ = igetter.get(folder_ + "Ele_eta_MTD_3");
  MonitorElement* meEle_phi_MTD_3_ = igetter.get(folder_ + "Ele_phi_MTD_3");

  MonitorElement* meEle_pt_MTD_4_ = igetter.get(folder_ + "Ele_pT_MTD_4");
  MonitorElement* meEle_eta_MTD_4_ = igetter.get(folder_ + "Ele_eta_MTD_4");
  MonitorElement* meEle_phi_MTD_4_ = igetter.get(folder_ + "Ele_phi_MTD_4");

  MonitorElement* meEle_pt_MTD_5_ = igetter.get(folder_ + "Ele_pT_MTD_5");
  MonitorElement* meEle_eta_MTD_5_ = igetter.get(folder_ + "Ele_eta_MTD_5");
  MonitorElement* meEle_phi_MTD_5_ = igetter.get(folder_ + "Ele_phi_MTD_5");

  MonitorElement* meEle_pt_MTD_6_ = igetter.get(folder_ + "Ele_pT_MTD_6");
  MonitorElement* meEle_eta_MTD_6_ = igetter.get(folder_ + "Ele_eta_MTD_6");
  MonitorElement* meEle_phi_MTD_6_ = igetter.get(folder_ + "Ele_phi_MTD_6");

  MonitorElement* meEle_pt_MTD_7_ = igetter.get(folder_ + "Ele_pT_MTD_7");
  MonitorElement* meEle_eta_MTD_7_ = igetter.get(folder_ + "Ele_eta_MTD_7");
  MonitorElement* meEle_phi_MTD_7_ = igetter.get(folder_ + "Ele_phi_MTD_7");


  if (!meBTLTrackEffEtaTot || !meBTLTrackEffPhiTot || !meBTLTrackEffPtTot || !meBTLTrackEffEtaMtd ||
      !meBTLTrackEffPhiMtd || !meBTLTrackEffPtMtd || !meETLTrackEffEtaTotZneg || !meETLTrackEffPhiTotZneg ||
      !meETLTrackEffPtTotZneg || !meETLTrackEffEtaMtdZneg || !meETLTrackEffPhiMtdZneg || !meETLTrackEffPtMtdZneg ||
      !meETLTrackEffEtaTotZpos || !meETLTrackEffPhiTotZpos || !meETLTrackEffPtTotZpos || !meETLTrackEffEtaMtdZpos || !meETLTrackEffPhiMtdZpos ||
      !meETLTrackEffPtMtdZpos ||  !meMVATrackEffPtTot || !meMVATrackMatchedEffPtTot || !meMVATrackMatchedEffPtMtd || !meMVATrackEffEtaTot ||
      !meMVATrackMatchedEffEtaTot || !meMVATrackMatchedEffEtaMtd ||!meEle_pt_tot || !meEle_pt_MTD_1_ || !meEle_pt_MTD_2_ ||
      !meEle_pt_MTD_3_ || !meEle_pt_MTD_4_ || !meEle_pt_MTD_5_ || !meEle_pt_MTD_6_ || !meEle_pt_MTD_7_ ||
      !meEle_pt_noMTD || !meEle_eta_tot || !meEle_eta_MTD_1_ || !meEle_eta_MTD_2_ || !meEle_eta_MTD_3_ ||
      !meEle_eta_MTD_4_ || !meEle_eta_MTD_5_ || !meEle_eta_MTD_6_ || !meEle_eta_MTD_7_ || !meEle_eta_noMTD || 
      !meEle_phi_tot || !meEle_phi_MTD_1_ || !meEle_phi_MTD_2_ || !meEle_phi_MTD_3_ || !meEle_phi_MTD_4_ ||
      !meEle_phi_MTD_5_ || !meEle_phi_MTD_6_ || !meEle_phi_MTD_7_ || !meEle_phi_noMTD) {
    edm::LogError("MtdTracksHarvester") << "Monitoring histograms not found!" << std::endl;
    return;
  }

  // --- Book  histograms
  ibook.cd(folder_);
  meBtlEtaEff_ = ibook.book1D("BtlEtaEff",
                              " Track Efficiency VS Eta;#eta;Efficiency",
                              meBTLTrackEffEtaTot->getNbinsX(),
                              meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                              meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBtlEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackEffEtaMtd, meBTLTrackEffEtaTot, meBtlEtaEff_);

  // ele iso addition starts; MTD vs noMTD case

  mePtEffMTD_1_ = ibook.book1D("pTeffMTD_1",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot->getNbinsX(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_1_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_1_, meEle_pt_tot, mePtEffMTD_1_);

  mePtEffMTD_2_ = ibook.book1D("pTeffMTD_2",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot->getNbinsX(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_2_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_2_, meEle_pt_tot, mePtEffMTD_2_);

  mePtEffMTD_3_ = ibook.book1D("pTeffMTD_3",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot->getNbinsX(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_3_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_3_, meEle_pt_tot, mePtEffMTD_3_);

  mePtEffMTD_4_ = ibook.book1D("pTeffMTD_4",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot->getNbinsX(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_4_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_4_, meEle_pt_tot, mePtEffMTD_4_);

  mePtEffMTD_5_ = ibook.book1D("pTeffMTD_5",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot->getNbinsX(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_5_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_5_, meEle_pt_tot, mePtEffMTD_5_);

  mePtEffMTD_6_ = ibook.book1D("pTeffMTD_6",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot->getNbinsX(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_6_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_6_, meEle_pt_tot, mePtEffMTD_6_);

  mePtEffMTD_7_ = ibook.book1D("pTeffMTD_7",
                              " MTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot->getNbinsX(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_7_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_MTD_7_, meEle_pt_tot, mePtEffMTD_7_);






  meEtaEffMTD_1_ = ibook.book1D("EtaEffMTD_1",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot->getNbinsX(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_1_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_1_, meEle_eta_tot, meEtaEffMTD_1_);

  meEtaEffMTD_2_ = ibook.book1D("EtaEffMTD_2",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot->getNbinsX(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_2_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_2_, meEle_eta_tot, meEtaEffMTD_2_);

  meEtaEffMTD_3_ = ibook.book1D("EtaEffMTD_3",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot->getNbinsX(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_3_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_3_, meEle_eta_tot, meEtaEffMTD_3_);

  meEtaEffMTD_4_ = ibook.book1D("EtaEffMTD_4",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot->getNbinsX(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_4_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_4_, meEle_eta_tot, meEtaEffMTD_4_);

  meEtaEffMTD_5_ = ibook.book1D("EtaEffMTD_5",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot->getNbinsX(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_5_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_5_, meEle_eta_tot, meEtaEffMTD_5_);

  meEtaEffMTD_6_ = ibook.book1D("EtaEffMTD_6",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot->getNbinsX(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_6_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_6_, meEle_eta_tot, meEtaEffMTD_6_);

  meEtaEffMTD_7_ = ibook.book1D("EtaEffMTD_7",
                              " MTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot->getNbinsX(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_7_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_MTD_7_, meEle_eta_tot, meEtaEffMTD_7_);





  mePhiEffMTD_1_ = ibook.book1D("PhiEffMTD_1",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot->getNbinsX(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_1_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_1_, meEle_phi_tot, mePhiEffMTD_1_);

  mePhiEffMTD_2_ = ibook.book1D("PhiEffMTD_2",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot->getNbinsX(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_2_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_2_, meEle_phi_tot, mePhiEffMTD_2_);

  mePhiEffMTD_3_ = ibook.book1D("PhiEffMTD_3",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot->getNbinsX(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_3_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_3_, meEle_phi_tot, mePhiEffMTD_3_);

  mePhiEffMTD_4_ = ibook.book1D("PhiEffMTD_4",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot->getNbinsX(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_4_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_4_, meEle_phi_tot, mePhiEffMTD_4_);

  mePhiEffMTD_5_ = ibook.book1D("PhiEffMTD_5",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot->getNbinsX(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_5_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_5_, meEle_phi_tot, mePhiEffMTD_5_);

  mePhiEffMTD_6_ = ibook.book1D("PhiEffMTD_6",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot->getNbinsX(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_6_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_6_, meEle_phi_tot, mePhiEffMTD_6_);

  mePhiEffMTD_7_ = ibook.book1D("PhiEffMTD_7",
                              " MTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot->getNbinsX(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_7_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_MTD_7_, meEle_phi_tot, mePhiEffMTD_7_);





  mePtEffnoMTD_ = ibook.book1D("pTeffnoMTD",
                              " noMTD isolation Efficiency VS pT;#pT;Efficiency",
                              meEle_pt_tot->getNbinsX(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_pt_tot->getTH1()->GetXaxis()->GetXmax());
  mePtEffnoMTD_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_pt_noMTD, meEle_pt_tot, mePtEffnoMTD_);

  meEtaEffnoMTD_ = ibook.book1D("EtaEffnoMTD",
                              " noMTD isolation Efficiency VS Eta;#eta;Efficiency",
                              meEle_eta_tot->getNbinsX(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_eta_tot->getTH1()->GetXaxis()->GetXmax());
  meEtaEffnoMTD_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_eta_noMTD, meEle_eta_tot, meEtaEffnoMTD_);

  mePhiEffnoMTD_ = ibook.book1D("PhiEffnoMTD",
                              " noMTD isolation Efficiency VS Phi;#phi;Efficiency",
                              meEle_phi_tot->getNbinsX(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmin(),
                              meEle_phi_tot->getTH1()->GetXaxis()->GetXmax());
  mePhiEffnoMTD_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meEle_phi_noMTD, meEle_phi_tot, mePhiEffnoMTD_);


  // Ele iso addition ends



  meBtlPhiEff_ = ibook.book1D("BtlPhiEff",
                              "Track Efficiency VS Phi;#phi [rad];Efficiency",
                              meBTLTrackEffPhiTot->getNbinsX(),
                              meBTLTrackEffPhiTot->getTH1()->GetXaxis()->GetXmin(),
                              meBTLTrackEffPhiTot->getTH1()->GetXaxis()->GetXmax());
  meBtlPhiEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackEffPhiMtd, meBTLTrackEffPhiTot, meBtlPhiEff_);

  meBtlPtEff_ = ibook.book1D("BtlPtEff",
                             "Track Efficiency VS Pt;Pt [GeV];Efficiency",
                             meBTLTrackEffPtTot->getNbinsX(),
                             meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                             meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBtlPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackEffPtMtd, meBTLTrackEffPtTot, meBtlPtEff_);

  meEtlEtaEff_[0] = ibook.book1D("EtlEtaEffZneg",
                                 " Track Efficiency VS Eta (-Z);#eta;Efficiency",
                                 meETLTrackEffEtaTotZneg->getNbinsX(),
                                 meETLTrackEffEtaTotZneg->getTH1()->GetXaxis()->GetXmin(),
                                 meETLTrackEffEtaTotZneg->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEff_[0]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEtaMtdZneg, meETLTrackEffEtaTotZneg, meEtlEtaEff_[0]);

  meEtlPhiEff_[0] = ibook.book1D("EtlPhiEffZneg",
                                 "Track Efficiency VS Phi (-Z);#phi [rad];Efficiency",
                                 meETLTrackEffPhiTotZneg->getNbinsX(),
                                 meETLTrackEffPhiTotZneg->getTH1()->GetXaxis()->GetXmin(),
                                 meETLTrackEffPhiTotZneg->getTH1()->GetXaxis()->GetXmax());
  meEtlPhiEff_[0]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPhiMtdZneg, meETLTrackEffPhiTotZneg, meEtlPhiEff_[0]);

  meEtlPtEff_[0] = ibook.book1D("EtlPtEffZneg",
                                "Track Efficiency VS Pt (-Z);Pt [GeV];Efficiency",
                                meETLTrackEffPtTotZneg->getNbinsX(),
                                meETLTrackEffPtTotZneg->getTH1()->GetXaxis()->GetXmin(),
                                meETLTrackEffPtTotZneg->getTH1()->GetXaxis()->GetXmax());
  meEtlPtEff_[0]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPtMtdZneg, meETLTrackEffPtTotZneg, meEtlPtEff_[0]);

  meEtlEtaEff_[1] = ibook.book1D("EtlEtaEffZpos",
                                 " Track Efficiency VS Eta (+Z);#eta;Efficiency",
                                 meETLTrackEffEtaTotZpos->getNbinsX(),
                                 meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                 meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEff_[1]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEtaMtdZpos, meETLTrackEffEtaTotZpos, meEtlEtaEff_[1]);

  meEtlPhiEff_[1] = ibook.book1D("EtlPhiEffZpos",
                                 "Track Efficiency VS Phi (+Z);#phi [rad];Efficiency",
                                 meETLTrackEffPhiTotZpos->getNbinsX(),
                                 meETLTrackEffPhiTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                 meETLTrackEffPhiTotZpos->getTH1()->GetXaxis()->GetXmax());
  meEtlPhiEff_[1]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPhiMtdZpos, meETLTrackEffPhiTotZpos, meEtlPhiEff_[1]);

  meEtlPtEff_[1] = ibook.book1D("EtlPtEffZpos",
                                "Track Efficiency VS Pt (+Z);Pt [GeV];Efficiency",
                                meETLTrackEffPtTotZpos->getNbinsX(),
                                meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmax());
  meEtlPtEff_[1]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPtMtdZpos, meETLTrackEffPtTotZpos, meEtlPtEff_[1]);

  

  meMVAPtSelEff_ = ibook.book1D("MVAPtSelEff",
                                "Track selected efficiency VS Pt;Pt [GeV];Efficiency",
                                meMVATrackEffPtTot->getNbinsX(),
                                meMVATrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                meMVATrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meMVAPtSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMVATrackMatchedEffPtTot, meMVATrackEffPtTot, meMVAPtSelEff_);

  meMVAEtaSelEff_ = ibook.book1D("MVAEtaSelEff",
                                 "Track selected efficiency VS Eta;Eta;Efficiency",
                                 meMVATrackEffEtaTot->getNbinsX(),
                                 meMVATrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                 meMVATrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meMVAEtaSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMVATrackMatchedEffEtaTot, meMVATrackEffEtaTot, meMVAEtaSelEff_);

  meMVAPtMatchEff_ = ibook.book1D("MVAPtMatchEff",
                                  "Track matched to GEN efficiency VS Pt;Pt [GeV];Efficiency",
                                  meMVATrackMatchedEffPtTot->getNbinsX(),
                                  meMVATrackMatchedEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                  meMVATrackMatchedEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meMVAPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMVATrackMatchedEffPtMtd, meMVATrackMatchedEffPtTot, meMVAPtMatchEff_);

  meMVAEtaMatchEff_ = ibook.book1D("MVAEtaMatchEff",
                                   "Track matched to GEN efficiency VS Eta;Eta;Efficiency",
                                   meMVATrackMatchedEffEtaTot->getNbinsX(),
                                   meMVATrackMatchedEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                   meMVATrackMatchedEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meMVAEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMVATrackMatchedEffEtaMtd, meMVATrackMatchedEffEtaTot, meMVAEtaMatchEff_);

  meBtlEtaEff_->getTH1()->SetMinimum(0.);
  meBtlPhiEff_->getTH1()->SetMinimum(0.);
  meBtlPtEff_->getTH1()->SetMinimum(0.);
  for (int i = 0; i < 2; i++) {
    meEtlEtaEff_[i]->getTH1()->SetMinimum(0.);
    meEtlPhiEff_[i]->getTH1()->SetMinimum(0.);
    meEtlPtEff_[i]->getTH1()->SetMinimum(0.);
  }
  meMVAPtSelEff_->getTH1()->SetMinimum(0.);
  meMVAEtaSelEff_->getTH1()->SetMinimum(0.);
  meMVAPtMatchEff_->getTH1()->SetMinimum(0.);
  meMVAEtaMatchEff_->getTH1()->SetMinimum(0.);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ----------
void MtdTracksHarvester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/Tracks/");

  descriptions.add("MtdTracksPostProcessor", desc);
}

DEFINE_FWK_MODULE(MtdTracksHarvester);
