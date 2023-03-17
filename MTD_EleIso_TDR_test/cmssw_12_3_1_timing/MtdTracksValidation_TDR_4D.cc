#include <string>
#include <tuple>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/Common/interface/ValidHandle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/GeantUnits.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"
#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "HepMC/GenRanges.h"
#include "CLHEP/Units/PhysicalConstants.h"


#include "DataFormats/EgammaCandidates/interface/GsfElectron.h" // Adding header files for electrons
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h" // Adding header files for electrons
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "MTD_ele_iso.h"
#include "MTD_ele_iso_TDR_3D.h"

class MtdTracksValidation : public DQMEDAnalyzer {
public:
  explicit MtdTracksValidation(const edm::ParameterSet&);
  ~MtdTracksValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  const bool mvaGenSel(const HepMC::GenParticle&, const float&);
  const bool mvaRecSel(const reco::TrackBase&, const reco::Vertex&, const double&, const double&);
  const bool mvaGenRecMatch(const HepMC::GenParticle&, const double&, const reco::TrackBase&);

  // ------------ member data ------------

  const std::string folder_;
  const float trackMinPt_;
  const float trackMinEta_;
  const float trackMaxEta_;

  static constexpr double etacutGEN_ = 4.;     // |eta| < 4;
  static constexpr double etacutREC_ = 3.;     // |eta| < 3;
  static constexpr double pTcut_ = 0.7;        // PT > 0.7 GeV
  static constexpr double deltaZcut_ = 0.1;    // dz separation 1 mm
  static constexpr double deltaPTcut_ = 0.05;  // dPT < 5%
  static constexpr double deltaDRcut_ = 0.03;  // DeltaR separation

  bool electron_iso_calc_;

  edm::EDGetTokenT<reco::TrackCollection> GenRecTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> RecTrackToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> RecVertexToken_;

  edm::EDGetTokenT<reco::GsfElectronCollection> GsfElectronToken_; // Adding token for electron collection 
  edm::EDGetTokenT<reco::GenParticleCollection> GenParticleToken_;


  edm::EDGetTokenT<edm::HepMCProduct> HepMCProductToken_;

  edm::EDGetTokenT<edm::ValueMap<int>> trackAssocToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> pathLengthToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmatmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> trackMVAQualToken_;

  edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdtopoToken_;
  edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particleTableToken_;

  MonitorElement* meBTLTrackRPTime_;
  MonitorElement* meBTLTrackEffEtaTot_;
  MonitorElement* meBTLTrackEffPhiTot_;
  MonitorElement* meBTLTrackEffPtTot_;
  MonitorElement* meBTLTrackEffEtaMtd_;
  MonitorElement* meBTLTrackEffPhiMtd_;
  MonitorElement* meBTLTrackEffPtMtd_;
  MonitorElement* meBTLTrackPtRes_;

  MonitorElement* meEleISO_Ntracks_; // Adding histograms for electrons (isolation stuff) no MTD case
  MonitorElement* meEleISO_chIso_;
  MonitorElement* meEleISO_rel_chIso_;

  MonitorElement* meEleISO_Ntracks_MTD_1_; // Adding histograms for electrons (isolation stuff) MTD case
  MonitorElement* meEleISO_chIso_MTD_1_;
  MonitorElement* meEleISO_rel_chIso_MTD_1_;

  MonitorElement* meEleISO_Ntracks_MTD_2_; 
  MonitorElement* meEleISO_chIso_MTD_2_;
  MonitorElement* meEleISO_rel_chIso_MTD_2_;

  MonitorElement* meEleISO_Ntracks_MTD_3_; 
  MonitorElement* meEleISO_chIso_MTD_3_;
  MonitorElement* meEleISO_rel_chIso_MTD_3_;

  MonitorElement* meEleISO_Ntracks_MTD_4_; 
  MonitorElement* meEleISO_chIso_MTD_4_;
  MonitorElement* meEleISO_rel_chIso_MTD_4_;

  MonitorElement* meEleISO_Ntracks_MTD_5_; 
  MonitorElement* meEleISO_chIso_MTD_5_;
  MonitorElement* meEleISO_rel_chIso_MTD_5_;

  MonitorElement* meEleISO_Ntracks_MTD_6_; 
  MonitorElement* meEleISO_chIso_MTD_6_;
  MonitorElement* meEleISO_rel_chIso_MTD_6_;

  MonitorElement* meEleISO_Ntracks_MTD_7_; 
  MonitorElement* meEleISO_chIso_MTD_7_;
  MonitorElement* meEleISO_rel_chIso_MTD_7_;

  MonitorElement* meEle_pt_tot_;
  MonitorElement* meEle_eta_tot_;
  MonitorElement* meEle_phi_tot_;
  MonitorElement* meEle_test_;
  MonitorElement* meEle_test_GenM_;
  MonitorElement* meEle_test_GenM_dR_;
  MonitorElement* meEle_test_ConeTracks_;
  MonitorElement* meEle_test_iVtx_;
  MonitorElement* meEle_test_nTracks_vtx_;
  MonitorElement* meEle_test_nEle_GenP_Z_perEvent_;
  MonitorElement* meEle_mother_test_;
  MonitorElement* meEle_Pvtx_match_;

  MonitorElement* meEle_Gen_Z_pt_;
  MonitorElement* meEle_Gen_Z_eta_;
  MonitorElement* meEle_Gen_Z_charge_;
  MonitorElement* meEle_Gen_Zmother_pt_;
  MonitorElement* meEle_Gen_Zmother_eta_;

  MonitorElement* meEle_Gen_Z_outsideEE_pt_;
  MonitorElement* meEle_Gen_Z_outsideEE_eta_;
  MonitorElement* meEle_Gen_Zmother_outsideEE_pt_;
  MonitorElement* meEle_Gen_Zmother_outsideEE_eta_;

  MonitorElement* meEle_track_pt_;
  MonitorElement* meEle_track_eta_;
  MonitorElement* meEle_track_phi_;
  MonitorElement* meEle_track_vz_;

  MonitorElement* meEle_pt_MTD_1_;
  MonitorElement* meEle_eta_MTD_1_;
  MonitorElement* meEle_phi_MTD_1_;

  MonitorElement* meEle_pt_MTD_2_;
  MonitorElement* meEle_eta_MTD_2_;
  MonitorElement* meEle_phi_MTD_2_;

  MonitorElement* meEle_pt_MTD_3_;
  MonitorElement* meEle_eta_MTD_3_;
  MonitorElement* meEle_phi_MTD_3_;

  MonitorElement* meEle_pt_MTD_4_;
  MonitorElement* meEle_eta_MTD_4_;
  MonitorElement* meEle_phi_MTD_4_;

  MonitorElement* meEle_pt_MTD_5_;
  MonitorElement* meEle_eta_MTD_5_;
  MonitorElement* meEle_phi_MTD_5_;

  MonitorElement* meEle_pt_MTD_6_;
  MonitorElement* meEle_eta_MTD_6_;
  MonitorElement* meEle_phi_MTD_6_;

  MonitorElement* meEle_pt_MTD_7_;
  MonitorElement* meEle_eta_MTD_7_;
  MonitorElement* meEle_phi_MTD_7_;

  MonitorElement* meEle_pt_noMTD_;
  MonitorElement* meEle_eta_noMTD_;
  MonitorElement* meEle_phi_noMTD_; // extra histograms end here





  MonitorElement* meETLTrackRPTime_;
  MonitorElement* meETLTrackEffEtaTot_[2];
  MonitorElement* meETLTrackEffPhiTot_[2];
  MonitorElement* meETLTrackEffPtTot_[2];
  MonitorElement* meETLTrackEffEtaMtd_[2];
  MonitorElement* meETLTrackEffPhiMtd_[2];
  MonitorElement* meETLTrackEffPtMtd_[2];
  MonitorElement* meETLTrackPtRes_;

  MonitorElement* meTracktmtd_;
  MonitorElement* meTrackt0Src_;
  MonitorElement* meTrackSigmat0Src_;
  MonitorElement* meTrackt0Pid_;
  MonitorElement* meTrackSigmat0Pid_;
  MonitorElement* meTrackt0SafePid_;
  MonitorElement* meTrackSigmat0SafePid_;
  MonitorElement* meTrackNumHits_;
  MonitorElement* meTrackMVAQual_;
  MonitorElement* meTrackPathLenghtvsEta_;

  MonitorElement* meMVATrackEffPtTot_;
  MonitorElement* meMVATrackMatchedEffPtTot_;
  MonitorElement* meMVATrackMatchedEffPtMtd_;
  MonitorElement* meMVATrackEffEtaTot_;
  MonitorElement* meMVATrackMatchedEffEtaTot_;
  MonitorElement* meMVATrackMatchedEffEtaMtd_;
  MonitorElement* meMVATrackResTot_;
  MonitorElement* meMVATrackPullTot_;
  MonitorElement* meMVATrackZposResTot_;
};

// ------------ constructor and destructor --------------
MtdTracksValidation::MtdTracksValidation(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      trackMinPt_(iConfig.getParameter<double>("trackMinimumPt")),
      trackMinEta_(iConfig.getParameter<double>("trackMinimumEta")),
      trackMaxEta_(iConfig.getParameter<double>("trackMaximumEta")),
      electron_iso_calc_(iConfig.getUntrackedParameter<bool>("optionalEleIso")) {
  GenRecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagG"));
  RecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagT"));
  RecVertexToken_ = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("inputTagV"));

  GsfElectronToken_ = consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("inputEle")); // Adding electron information to the EleToken
  GenParticleToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("inputGenP"));

  HepMCProductToken_ = consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("inputTagH"));
  trackAssocToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("trackAssocSrc"));
  pathLengthToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pathLengthSrc"));
  tmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tmtd"));
  SigmatmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatmtd"));
  t0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"));
  Sigmat0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"));
  t0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0PID"));
  Sigmat0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0PID"));
  t0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0SafePID"));
  Sigmat0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0SafePID"));
  trackMVAQualToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMVAQual"));
  mtdtopoToken_ = esConsumes<MTDTopology, MTDTopologyRcd>();
  particleTableToken_ = esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>();
}

MtdTracksValidation::~MtdTracksValidation() {}

// ------------ method called for each event  ------------
void MtdTracksValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace geant_units::operators;
  using namespace std;

  auto topologyHandle = iSetup.getTransientHandle(mtdtopoToken_);
  const MTDTopology* topology = topologyHandle.product();

  bool topo1Dis = false;
  bool topo2Dis = false;
  if (topology->getMTDTopologyMode() <= static_cast<int>(MTDTopologyMode::Mode::barphiflat)) {
    topo1Dis = true;
  }
  if (topology->getMTDTopologyMode() > static_cast<int>(MTDTopologyMode::Mode::barphiflat)) {
    topo2Dis = true;
  }

  auto GenRecTrackHandle = makeValid(iEvent.getHandle(GenRecTrackToken_));
  auto RecVertexHandle = makeValid(iEvent.getHandle(RecVertexToken_));

  const auto& tMtd = iEvent.get(tmtdToken_);
  const auto& SigmatMtd = iEvent.get(SigmatmtdToken_);
  const auto& t0Src = iEvent.get(t0SrcToken_);
  const auto& Sigmat0Src = iEvent.get(Sigmat0SrcToken_);
  const auto& t0Pid = iEvent.get(t0PidToken_);
  const auto& Sigmat0Pid = iEvent.get(Sigmat0PidToken_);
  const auto& t0Safe = iEvent.get(t0SafePidToken_);
  const auto& Sigmat0Safe = iEvent.get(Sigmat0SafePidToken_);
  const auto& mtdQualMVA = iEvent.get(trackMVAQualToken_);
  const auto& trackAssoc = iEvent.get(trackAssocToken_);
  const auto& pathLength = iEvent.get(pathLengthToken_);

 
 
  if(electron_iso_calc_){
    
    float min_dR_cut = 0.01;
    float max_dR_cut = 0.3;
    //float min_pt_cut = 0.7; // 0.4 for endcap
    float min_pt_cut = 0.4;
    //float max_dz_cut = 0.1; // will check 0.15/0.20/0.25/0.30 for tracks in iso cone // 0.1 for barrel, 0.2 for endcap
    float max_dz_cut = 0.4;
    float max_dz_vtx_cut = 0.5;
    float max_dxy_vtx_cut = 0.2; 
    std::vector<float> max_dt_vtx_cut{0.30,0.20,0.15,0.10,0.08,0.06,0.04}; // has to be 7 values here!!! Code created for 7 dt values!! If change iso_object if number of dt values change.
    std::vector<float> max_dt_track_cut{9999};
    float min_strip_cut = 0.01;
    float min_track_mtd_mva_cut = 0.5;

    int n_Gen_matched = 0;

    meEle_test_->Fill(1);

    

    //iEvent.getHandle(GenRecTrackToken_)
    //GenRecTrackHandle
    MTD_Ele_iso iso_object(max_dR_cut,
                          min_dR_cut,
                          min_pt_cut,
                          max_dz_cut,
                          max_dz_vtx_cut,
                          max_dxy_vtx_cut,
                          max_dt_vtx_cut,
                          max_dt_track_cut,
                          min_strip_cut,
                          min_track_mtd_mva_cut,
                          iEvent.getHandle(GenRecTrackToken_),
                          iEvent.get(t0PidToken_),
                          iEvent.get(Sigmat0PidToken_),
                          iEvent.get(trackMVAQualToken_),
                          //Vtx_chosen,
                          iEvent.getHandle(RecVertexToken_));
                          //iEvent.getHandle(GenParticleToken_));
    
    auto eleHandle = makeValid(iEvent.getHandle(GsfElectronToken_));
    reco::GsfElectronCollection eleColl = *(eleHandle.product());

    auto GenPartHandle = makeValid(iEvent.getHandle(GenParticleToken_));
    reco::GenParticleCollection GenPartColl = *(GenPartHandle.product());


    int n_Gen_Z_electrons = 0;

    for(const auto& genParticle_check : GenPartColl){
      
      if(genParticle_check.pdgId() == 11 || genParticle_check.pdgId() == -11){
        
        if(genParticle_check.mother()->pdgId() == 23){

          meEle_Gen_Z_pt_->Fill(genParticle_check.pt());
          meEle_Gen_Z_eta_->Fill(genParticle_check.eta());
          meEle_Gen_Z_charge_->Fill(genParticle_check.charge());
          meEle_Gen_Zmother_pt_->Fill(genParticle_check.mother()->pt());
          meEle_Gen_Zmother_eta_->Fill(genParticle_check.mother()->eta());
          
          ++n_Gen_Z_electrons;

          if(fabs(genParticle_check.eta()) > 2.4){
            
            meEle_Gen_Z_outsideEE_pt_->Fill(genParticle_check.pt());
            meEle_Gen_Z_outsideEE_eta_->Fill(genParticle_check.eta());
            meEle_Gen_Zmother_outsideEE_pt_->Fill(genParticle_check.mother()->pt());
            meEle_Gen_Zmother_outsideEE_eta_->Fill(genParticle_check.mother()->eta());

          }

        }
      }
    }

    meEle_test_nEle_GenP_Z_perEvent_->Fill(n_Gen_Z_electrons);
    
    
    for (const auto& ele : eleColl){

      meEle_test_->Fill(2);
      
      if( ele.pt()> 10 && fabs(ele.eta()) < 2.4 ){

        math::XYZVector EleMomentum = ele.momentum();

        bool eleMatch_found = false;

        meEle_test_->Fill(3);

        // GenP matching
        for(const auto& genParticle : GenPartColl){
          
          if(genParticle.pdgId() == 11 || genParticle.pdgId() == -11){

            //meEle_mother_test_->Fill(genParticle.mother()->pdgId());

            if(genParticle.mother()->pdgId() == 23){

              math::XYZVector GenPartMomentum = genParticle.momentum();
              double dr_match = reco::deltaR(GenPartMomentum, EleMomentum);
              
              if((genParticle.pdgId() == 11 && ele.charge() == -1) || (genParticle.pdgId() == -11 && ele.charge() == 1)){

                meEle_test_GenM_dR_->Fill(dr_match);
                                   
                if(dr_match < 0.10){ // from Egamma 0.05 value is taken, will try 0.10
                  eleMatch_found = true;
                  ++n_Gen_matched;
                  break;
                }
                    
              }else{
                continue;
              }

            }else{
              continue;
            } 

          }else{
            continue;
          }        
        }

        if(eleMatch_found){


          meEle_test_->Fill(8);
        
        
          std::tuple<int,float,float,int,float,float,int,float,float,int,float,float,int,float,float,int,float,float,int,float,float,int,float,float,int,int,int> iso_return;
          
          iso_return = iso_object.ele_iso(&(*ele.gsfTrack()));

          meEle_pt_tot_->Fill(ele.pt());
          meEle_eta_tot_->Fill(ele.eta());
          meEle_phi_tot_->Fill(ele.phi());


          //meEle_track_pt_->Fill(eleTrack_test->pt());
          //meEle_track_eta_->Fill(eleTrack_test->eta());
          //meEle_track_phi_->Fill(eleTrack_test->phi());
          //meEle_track_vz_->Fill(eleTrack_test->vz());

          meEleISO_Ntracks_->Fill(std::get<0>(iso_return)); // Filling hists for Ntraks and chIso sums for noMTD case //
          meEleISO_chIso_->Fill(std::get<1>(iso_return));
          meEleISO_rel_chIso_->Fill(std::get<2>(iso_return));

          meEleISO_Ntracks_MTD_1_->Fill(std::get<3>(iso_return)); // Filling hists for Ntraks and chIso sums for MTD case//
          meEleISO_chIso_MTD_1_->Fill(std::get<4>(iso_return));
          meEleISO_rel_chIso_MTD_1_->Fill(std::get<5>(iso_return)); // 

          meEleISO_Ntracks_MTD_2_->Fill(std::get<6>(iso_return)); 
          meEleISO_chIso_MTD_2_->Fill(std::get<7>(iso_return));
          meEleISO_rel_chIso_MTD_2_->Fill(std::get<8>(iso_return));

          meEleISO_Ntracks_MTD_3_->Fill(std::get<9>(iso_return)); 
          meEleISO_chIso_MTD_3_->Fill(std::get<10>(iso_return));
          meEleISO_rel_chIso_MTD_3_->Fill(std::get<11>(iso_return));

          meEleISO_Ntracks_MTD_4_->Fill(std::get<12>(iso_return)); 
          meEleISO_chIso_MTD_4_->Fill(std::get<13>(iso_return));
          meEleISO_rel_chIso_MTD_4_->Fill(std::get<14>(iso_return));

          meEleISO_Ntracks_MTD_5_->Fill(std::get<15>(iso_return)); 
          meEleISO_chIso_MTD_5_->Fill(std::get<16>(iso_return));
          meEleISO_rel_chIso_MTD_5_->Fill(std::get<17>(iso_return));

          meEleISO_Ntracks_MTD_6_->Fill(std::get<18>(iso_return)); 
          meEleISO_chIso_MTD_6_->Fill(std::get<19>(iso_return));
          meEleISO_rel_chIso_MTD_6_->Fill(std::get<20>(iso_return));

          meEleISO_Ntracks_MTD_7_->Fill(std::get<21>(iso_return)); 
          meEleISO_chIso_MTD_7_->Fill(std::get<22>(iso_return));
          meEleISO_rel_chIso_MTD_7_->Fill(std::get<23>(iso_return));

          
          meEle_test_ConeTracks_->Fill(std::get<24>(iso_return));
          meEle_test_iVtx_->Fill(std::get<25>(iso_return));
          meEle_test_nTracks_vtx_->Fill(std::get<26>(iso_return));


          if(std::get<2>(iso_return) < 0.08){
            meEle_pt_noMTD_->Fill(ele.pt());
            meEle_eta_noMTD_->Fill(ele.eta());
            meEle_phi_noMTD_->Fill(ele.phi());
          }
          
          if(std::get<5>(iso_return) < 0.08){
            meEle_pt_MTD_1_->Fill(ele.pt());
            meEle_eta_MTD_1_->Fill(ele.eta());
            meEle_phi_MTD_1_->Fill(ele.phi());
          }

          if(std::get<8>(iso_return) < 0.08){
            meEle_pt_MTD_2_->Fill(ele.pt());
            meEle_eta_MTD_2_->Fill(ele.eta());
            meEle_phi_MTD_2_->Fill(ele.phi());
          }

          if(std::get<11>(iso_return) < 0.08){
            meEle_pt_MTD_3_->Fill(ele.pt());
            meEle_eta_MTD_3_->Fill(ele.eta());
            meEle_phi_MTD_3_->Fill(ele.phi());
          }

          if(std::get<14>(iso_return) < 0.08){
            meEle_pt_MTD_4_->Fill(ele.pt());
            meEle_eta_MTD_4_->Fill(ele.eta());
            meEle_phi_MTD_4_->Fill(ele.phi());
          }

          if(std::get<17>(iso_return) < 0.08){
            meEle_pt_MTD_5_->Fill(ele.pt());
            meEle_eta_MTD_5_->Fill(ele.eta());
            meEle_phi_MTD_5_->Fill(ele.phi());
          }

          if(std::get<20>(iso_return) < 0.08){
            meEle_pt_MTD_6_->Fill(ele.pt());
            meEle_eta_MTD_6_->Fill(ele.eta());
            meEle_phi_MTD_6_->Fill(ele.phi());
          }

          if(std::get<23>(iso_return) < 0.08){
            meEle_pt_MTD_7_->Fill(ele.pt());
            meEle_eta_MTD_7_->Fill(ele.eta());
            meEle_phi_MTD_7_->Fill(ele.phi());
          }

        } // genP matching cut

      } // ele pt and eta cut

    } // electron collection inside single event end

    meEle_test_GenM_->Fill(n_Gen_matched);
 
  } // Bool iso statement end
 
 
 
  unsigned int index = 0;
  // --- Loop over all RECO tracks ---
  for (const auto& trackGen : *GenRecTrackHandle) {
    const reco::TrackRef trackref(iEvent.getHandle(GenRecTrackToken_), index);
    index++;

    if (trackAssoc[trackref] == -1) {
      LogInfo("mtdTracks") << "Extended track not associated";
      continue;
    }

    const reco::TrackRef mtdTrackref = reco::TrackRef(iEvent.getHandle(RecTrackToken_), trackAssoc[trackref]);
    const reco::Track track = *mtdTrackref;

    if (track.pt() < trackMinPt_)
      continue;

    meTracktmtd_->Fill(tMtd[trackref]);
    if (std::round(SigmatMtd[trackref] - Sigmat0Pid[trackref]) != 0) {
      LogWarning("mtdTracks") << "TimeError associated to refitted track is different from TimeError stored in tofPID "
                                 "sigmat0 ValueMap: this should not happen";
    }

    meTrackt0Src_->Fill(t0Src[trackref]);
    meTrackSigmat0Src_->Fill(Sigmat0Src[trackref]);

    meTrackt0Pid_->Fill(t0Pid[trackref]);
    meTrackSigmat0Pid_->Fill(Sigmat0Pid[trackref]);
    meTrackt0SafePid_->Fill(t0Safe[trackref]);
    meTrackSigmat0SafePid_->Fill(Sigmat0Safe[trackref]);
    meTrackMVAQual_->Fill(mtdQualMVA[trackref]);

    meTrackPathLenghtvsEta_->Fill(std::abs(track.eta()), pathLength[trackref]);

    if (std::abs(track.eta()) < trackMinEta_) {
      // --- all BTL tracks (with and without hit in MTD) ---
      meBTLTrackEffEtaTot_->Fill(track.eta());
      meBTLTrackEffPhiTot_->Fill(track.phi());
      meBTLTrackEffPtTot_->Fill(track.pt());

      bool MTDBtl = false;
      int numMTDBtlvalidhits = 0;
      for (const auto hit : track.recHits()) {
        if (hit->isValid() == false)
          continue;
        MTDDetId Hit = hit->geographicalId();
        if ((Hit.det() == 6) && (Hit.subdetId() == 1) && (Hit.mtdSubDetector() == 1)) {
          MTDBtl = true;
          numMTDBtlvalidhits++;
        }
      }
      meTrackNumHits_->Fill(numMTDBtlvalidhits);

      // --- keeping only tracks with last hit in MTD ---
      if (MTDBtl == true) {
        meBTLTrackEffEtaMtd_->Fill(track.eta());
        meBTLTrackEffPhiMtd_->Fill(track.phi());
        meBTLTrackEffPtMtd_->Fill(track.pt());
        meBTLTrackRPTime_->Fill(track.t0());
        meBTLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
      }
    }  //loop over (geometrical) BTL tracks

    else {
      // --- all ETL tracks (with and without hit in MTD) ---
      if ((track.eta() < -trackMinEta_) && (track.eta() > -trackMaxEta_)) {
        meETLTrackEffEtaTot_[0]->Fill(track.eta());
        meETLTrackEffPhiTot_[0]->Fill(track.phi());
        meETLTrackEffPtTot_[0]->Fill(track.pt());
      }

      if ((track.eta() > trackMinEta_) && (track.eta() < trackMaxEta_)) {
        meETLTrackEffEtaTot_[1]->Fill(track.eta());
        meETLTrackEffPhiTot_[1]->Fill(track.phi());
        meETLTrackEffPtTot_[1]->Fill(track.pt());
      }

      bool MTDEtlZnegD1 = false;
      bool MTDEtlZnegD2 = false;
      bool MTDEtlZposD1 = false;
      bool MTDEtlZposD2 = false;
      int numMTDEtlvalidhits = 0;
      for (const auto hit : track.recHits()) {
        if (hit->isValid() == false)
          continue;
        MTDDetId Hit = hit->geographicalId();
        if ((Hit.det() == 6) && (Hit.subdetId() == 1) && (Hit.mtdSubDetector() == 2)) {
          ETLDetId ETLHit = hit->geographicalId();

          if (topo2Dis) {
            if ((ETLHit.zside() == -1) && (ETLHit.nDisc() == 1)) {
              MTDEtlZnegD1 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == -1) && (ETLHit.nDisc() == 2)) {
              MTDEtlZnegD2 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == 1) && (ETLHit.nDisc() == 1)) {
              MTDEtlZposD1 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == 1) && (ETLHit.nDisc() == 2)) {
              MTDEtlZposD2 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              numMTDEtlvalidhits++;
            }
          }

          if (topo1Dis) {
            if (ETLHit.zside() == -1) {
              MTDEtlZnegD1 = true;
              meETLTrackRPTime_->Fill(track.t0());
              numMTDEtlvalidhits++;
            }
            if (ETLHit.zside() == 1) {
              MTDEtlZposD1 = true;
              meETLTrackRPTime_->Fill(track.t0());
              numMTDEtlvalidhits++;
            }
          }
        }
      }
      meTrackNumHits_->Fill(-numMTDEtlvalidhits);

      // --- keeping only tracks with last hit in MTD ---
      if ((track.eta() < -trackMinEta_) && (track.eta() > -trackMaxEta_)) {
        if ((MTDEtlZnegD1 == true) || (MTDEtlZnegD2 == true)) {
          meETLTrackEffEtaMtd_[0]->Fill(track.eta());
          meETLTrackEffPhiMtd_[0]->Fill(track.phi());
          meETLTrackEffPtMtd_[0]->Fill(track.pt());
        }
      }
      if ((track.eta() > trackMinEta_) && (track.eta() < trackMaxEta_)) {
        if ((MTDEtlZposD1 == true) || (MTDEtlZposD2 == true)) {
          meETLTrackEffEtaMtd_[1]->Fill(track.eta());
          meETLTrackEffPhiMtd_[1]->Fill(track.phi());
          meETLTrackEffPtMtd_[1]->Fill(track.pt());
        }
      }
    }
  }  //RECO tracks loop

  // reco-gen matching used for MVA quality flag
  const auto& primRecoVtx = *(RecVertexHandle.product()->begin());

  auto GenEventHandle = makeValid(iEvent.getHandle(HepMCProductToken_));
  const HepMC::GenEvent* mc = GenEventHandle->GetEvent();
  double zsim = convertMmToCm((*(mc->vertices_begin()))->position().z());
  double tsim = (*(mc->vertices_begin()))->position().t() * CLHEP::mm / CLHEP::c_light;

  auto pdt = iSetup.getHandle(particleTableToken_);
  const HepPDT::ParticleDataTable* pdTable = pdt.product();

  // select events with reco vertex close to true simulated primary vertex
  if (std::abs(primRecoVtx.z() - zsim) < deltaZcut_) {
    index = 0;
    for (const auto& trackGen : *GenRecTrackHandle) {
      const reco::TrackRef trackref(iEvent.getHandle(GenRecTrackToken_), index);
      index++;

      // select the reconstructed track

      if (trackAssoc[trackref] == -1) {
        continue;
      }

      if (mvaRecSel(trackGen, primRecoVtx, t0Safe[trackref], Sigmat0Safe[trackref])) {
        meMVATrackEffPtTot_->Fill(trackGen.pt());
        meMVATrackEffEtaTot_->Fill(std::abs(trackGen.eta()));

        double dZ = trackGen.vz() - zsim;
        double dT(-9999.);
        double pullT(-9999.);
        if (Sigmat0Safe[trackref] != -1.) {
          dT = t0Safe[trackref] - tsim;
          pullT = dT / Sigmat0Safe[trackref];
        }
        for (const auto& genP : mc->particle_range()) {
          // select status 1 genParticles and match them to the reconstructed track

          float charge = pdTable->particle(HepPDT::ParticleID(genP->pdg_id())) != nullptr
                             ? pdTable->particle(HepPDT::ParticleID(genP->pdg_id()))->charge()
                             : 0.f;
          if (mvaGenSel(*genP, charge)) {
            if (mvaGenRecMatch(*genP, zsim, trackGen)) {
              meMVATrackZposResTot_->Fill(dZ);
              meMVATrackMatchedEffPtTot_->Fill(trackGen.pt());
              meMVATrackMatchedEffEtaTot_->Fill(std::abs(trackGen.eta()));
              if (pullT > -9999.) {
                meMVATrackResTot_->Fill(dT);
                meMVATrackPullTot_->Fill(pullT);
                meMVATrackMatchedEffPtMtd_->Fill(trackGen.pt());
                meMVATrackMatchedEffEtaMtd_->Fill(std::abs(trackGen.eta()));
              }
              break;
            }
          }
        }
      }
    }
  }
}

// ------------ method for histogram booking ------------
void MtdTracksValidation::bookHistograms(DQMStore::IBooker& ibook, edm::Run const& run, edm::EventSetup const& iSetup) {
  ibook.setCurrentFolder(folder_);

  // histogram booking
  meBTLTrackRPTime_ = ibook.book1D("TrackBTLRPTime", "Track t0 with respect to R.P.;t0 [ns]", 100, -1, 3);
  meBTLTrackEffEtaTot_ = ibook.book1D("TrackBTLEffEtaTot", "Track efficiency vs eta (Tot);#eta_{RECO}", 100, -1.6, 1.6);
  meBTLTrackEffPhiTot_ =
      ibook.book1D("TrackBTLEffPhiTot", "Track efficiency vs phi (Tot);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meBTLTrackEffPtTot_ = ibook.book1D("TrackBTLEffPtTot", "Track efficiency vs pt (Tot);pt_{RECO} [GeV]", 50, 0, 10);
  meBTLTrackEffEtaMtd_ = ibook.book1D("TrackBTLEffEtaMtd", "Track efficiency vs eta (Mtd);#eta_{RECO}", 100, -1.6, 1.6);
  meBTLTrackEffPhiMtd_ =
      ibook.book1D("TrackBTLEffPhiMtd", "Track efficiency vs phi (Mtd);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meBTLTrackEffPtMtd_ = ibook.book1D("TrackBTLEffPtMtd", "Track efficiency vs pt (Mtd);pt_{RECO} [GeV]", 50, 0, 10);
  


  meEleISO_Ntracks_ = ibook.book1D("Ele_Iso_Ntracks", "Tracks in isolation cone around electron track after basic cuts", 20, 0, 20); // hists for electrons
  meEleISO_chIso_ = ibook.book1D("Ele_chIso_sum", "Track pT sum in isolation cone around electron track after basic cuts", 400, 0, 20);
  meEleISO_rel_chIso_ = ibook.book1D("Ele_rel_chIso_sum", "Track relative pT sum in isolation cone around electron track after basic cuts", 200, 0, 4);

  meEleISO_Ntracks_MTD_1_ = ibook.book1D("Ele_Iso_Ntracks_MTD_1", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_1_ = ibook.book1D("Ele_chIso_sum_MTD_1", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_1_ = ibook.book1D("Ele_rel_chIso_sum_MTD_1", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_2_ = ibook.book1D("Ele_Iso_Ntracks_MTD_2", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_2_ = ibook.book1D("Ele_chIso_sum_MTD_2", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_2_ = ibook.book1D("Ele_rel_chIso_sum_MTD_2", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_3_ = ibook.book1D("Ele_Iso_Ntracks_MTD_3", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_3_ = ibook.book1D("Ele_chIso_sum_MTD_3", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_3_ = ibook.book1D("Ele_rel_chIso_sum_MTD_3", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_4_ = ibook.book1D("Ele_Iso_Ntracks_MTD_4", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_4_ = ibook.book1D("Ele_chIso_sum_MTD_4", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_4_ = ibook.book1D("Ele_rel_chIso_sum_MTD_4", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_5_ = ibook.book1D("Ele_Iso_Ntracks_MTD_5", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_5_ = ibook.book1D("Ele_chIso_sum_MTD_5", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_5_ = ibook.book1D("Ele_rel_chIso_sum_MTD_5", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_6_ = ibook.book1D("Ele_Iso_Ntracks_MTD_6", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_6_ = ibook.book1D("Ele_chIso_sum_MTD_6", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_6_ = ibook.book1D("Ele_rel_chIso_sum_MTD_6", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_7_ = ibook.book1D("Ele_Iso_Ntracks_MTD_7", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_7_ = ibook.book1D("Ele_chIso_sum_MTD_7", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_7_ = ibook.book1D("Ele_rel_chIso_sum_MTD_7", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);


  meEle_pt_tot_ = ibook.book1D("Ele_pT_tot", "Electron pT tot", 30, 10, 100); // hists for ele isto stuff start
  meEle_pt_noMTD_ = ibook.book1D("Ele_pT_noMTD", "Electron pT noMTD", 30, 10, 100);

  meEle_eta_tot_ = ibook.book1D("Ele_eta_tot", "Electron eta tot", 128, -3.2, 3.2);
  meEle_eta_noMTD_ = ibook.book1D("Ele_eta_noMTD", "Electron eta noMTD", 128, -3.2, 3.2);

  meEle_phi_tot_ = ibook.book1D("Ele_phi_tot", "Electron phi tot", 128, -3.2, 3.2);
  meEle_phi_noMTD_ = ibook.book1D("Ele_phi_noMTD", "Electron phi noMTD", 128, -3.2, 3.2);

  
  meEle_pt_MTD_1_ = ibook.book1D("Ele_pT_MTD_1", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_1_ = ibook.book1D("Ele_eta_MTD_1", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_1_ = ibook.book1D("Ele_phi_MTD_1", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_2_ = ibook.book1D("Ele_pT_MTD_2", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_2_ = ibook.book1D("Ele_eta_MTD_2", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_2_ = ibook.book1D("Ele_phi_MTD_2", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_3_ = ibook.book1D("Ele_pT_MTD_3", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_3_ = ibook.book1D("Ele_eta_MTD_3", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_3_ = ibook.book1D("Ele_phi_MTD_3", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_4_ = ibook.book1D("Ele_pT_MTD_4", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_4_ = ibook.book1D("Ele_eta_MTD_4", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_4_ = ibook.book1D("Ele_phi_MTD_4", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_5_ = ibook.book1D("Ele_pT_MTD_5", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_5_ = ibook.book1D("Ele_eta_MTD_5", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_5_ = ibook.book1D("Ele_phi_MTD_5", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_6_ = ibook.book1D("Ele_pT_MTD_6", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_6_ = ibook.book1D("Ele_eta_MTD_6", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_6_ = ibook.book1D("Ele_phi_MTD_6", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_7_ = ibook.book1D("Ele_pT_MTD_7", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_7_ = ibook.book1D("Ele_eta_MTD_7", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_7_ = ibook.book1D("Ele_phi_MTD_7", "Electron phi MTD", 128, -3.2, 3.2);

  //meEle_track_pt_ = ibook.book1D("Ele_track_pT", "Electron track pT ", 100, 0, 100);
  //meEle_track_eta_ = ibook.book1D("Ele_track_eta", "Electron track eta", 100, -3.2, 3.2);
  //meEle_track_phi_ = ibook.book1D("Ele_track_phi", "Electron track phi", 100, -3.2, 3.2);
  //meEle_track_vz_ = ibook.book1D("Ele_track_vz", "Electron track vz", 100, -10, 10);

  meEle_test_ = ibook.book1D("Ele_test", "Test values", 10, 0, 10);
  meEle_test_GenM_ = ibook.book1D("Ele_GenMatching_test","Checking how many electrons per event are GenMatched",5,0,5);
  meEle_test_GenM_dR_ = ibook.book1D("Ele_test_GenM_dR_test","Calculated dR between RECO and GenP electron, after the charge is the same and mother is Z",160,0,3.2);
  meEle_test_ConeTracks_ = ibook.book1D("Ele_ConeTracks_test","The amount of tracks inside electron isolation cone without any dz/dxy/dt/pT cuts",100,0,100);
  meEle_test_iVtx_ = ibook.book1D("Ele_test_iVtx","Which vertex by count was chosen to be PV from full vtx collection",10,0,10);
  meEle_test_nTracks_vtx_ = ibook.book1D("Ele_test_nTracks_vtx","How many tracks are used to reconstruct this vertex",100,0,100);
  meEle_test_nEle_GenP_Z_perEvent_ = ibook.book1D("Ele_test_nEle_GenP_Z_perEvent","How man GenP Z electrons are in the event",5,0,5);
  meEle_mother_test_ = ibook.book1D("Ele_gen_mother_test", "pdgId values",200,0,200);

  meEle_Gen_Z_pt_ = ibook.book1D("Ele_Gen_Z_pt","GenP Z electron pT",100,0,100);
  meEle_Gen_Z_eta_ = ibook.book1D("Ele_Gen_Z_eta","GenP Z electron eta",280,-7.0,7.0);
  meEle_Gen_Z_charge_ = ibook.book1D("Ele_Gen_Z_charge","GenP Z electron charge",3,-1,2);
  meEle_Gen_Zmother_pt_ = ibook.book1D("Ele_Gen_Zmother_pt","GenP Z pT",100,0,100);
  meEle_Gen_Zmother_eta_ = ibook.book1D("Ele_Gen_Zmother_eta","GenP Z eta",280,-7.0,7);

  meEle_Gen_Z_outsideEE_pt_ = ibook.book1D("Ele_Gen_Z_outsideEE_pt","GenP Z electron pT outside endcap acceptance",100,0,100);
  meEle_Gen_Z_outsideEE_eta_ = ibook.book1D("Ele_Gen_Z__outsideEEeta","GenP Z electron eta outside endcap acceptance",128,-3.2,3.2);
  meEle_Gen_Zmother_outsideEE_pt_ = ibook.book1D("Ele_Gen_Zmother_outsideEE_pt","GenP Z pT, if electron was outside endcap acceptance",100,0,100);
  meEle_Gen_Zmother_outsideEE_eta_ = ibook.book1D("Ele_Gen_Zmother_outsideEE_eta","GenP Z eta if electron was outside endcap acceptance",32,-3.2,3.2);
  //meEle_track_eta_ = ibook.book1D("Ele_track_eta", "Electron track eta", 100, -3.2, 3.2);
  //meEle_Pvtx_match_ = ibook.book1D("Ele_Vtx_match","Electron track match to primary vertex",10,0,10); // hists for ele is end
  
   
  
  meBTLTrackPtRes_ =
      ibook.book1D("TrackBTLPtRes", "Track pT resolution  ;pT_{Gentrack}-pT_{MTDtrack}/pT_{Gentrack} ", 100, -0.1, 0.1);
  meETLTrackRPTime_ = ibook.book1D("TrackETLRPTime", "Track t0 with respect to R.P.;t0 [ns]", 100, -1, 3);
  meETLTrackEffEtaTot_[0] =
      ibook.book1D("TrackETLEffEtaTotZneg", "Track efficiency vs eta (Tot) (-Z);#eta_{RECO}", 100, -3.2, -1.4);
  meETLTrackEffEtaTot_[1] =
      ibook.book1D("TrackETLEffEtaTotZpos", "Track efficiency vs eta (Tot) (+Z);#eta_{RECO}", 100, 1.4, 3.2);
  meETLTrackEffPhiTot_[0] =
      ibook.book1D("TrackETLEffPhiTotZneg", "Track efficiency vs phi (Tot) (-Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPhiTot_[1] =
      ibook.book1D("TrackETLEffPhiTotZpos", "Track efficiency vs phi (Tot) (+Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPtTot_[0] =
      ibook.book1D("TrackETLEffPtTotZneg", "Track efficiency vs pt (Tot) (-Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffPtTot_[1] =
      ibook.book1D("TrackETLEffPtTotZpos", "Track efficiency vs pt (Tot) (+Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffEtaMtd_[0] =
      ibook.book1D("TrackETLEffEtaMtdZneg", "Track efficiency vs eta (Mtd) (-Z);#eta_{RECO}", 100, -3.2, -1.4);
  meETLTrackEffEtaMtd_[1] =
      ibook.book1D("TrackETLEffEtaMtdZpos", "Track efficiency vs eta (Mtd) (+Z);#eta_{RECO}", 100, 1.4, 3.2);
  meETLTrackEffPhiMtd_[0] =
      ibook.book1D("TrackETLEffPhiMtdZneg", "Track efficiency vs phi (Mtd) (-Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPhiMtd_[1] =
      ibook.book1D("TrackETLEffPhiMtdZpos", "Track efficiency vs phi (Mtd) (+Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPtMtd_[0] =
      ibook.book1D("TrackETLEffPtMtdZneg", "Track efficiency vs pt (Mtd) (-Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffPtMtd_[1] =
      ibook.book1D("TrackETLEffPtMtdZpos", "Track efficiency vs pt (Mtd) (+Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackPtRes_ =
      ibook.book1D("TrackETLPtRes", "Track pT resolution;pT_{Gentrack}-pT_{MTDtrack}/pT_{Gentrack} ", 100, -0.1, 0.1);

  meTracktmtd_ = ibook.book1D("Tracktmtd", "Track time from TrackExtenderWithMTD;tmtd [ns]", 150, 1, 16);
  meTrackt0Src_ = ibook.book1D("Trackt0Src", "Track time from TrackExtenderWithMTD;t0Src [ns]", 100, -1.5, 1.5);
  meTrackSigmat0Src_ =
      ibook.book1D("TrackSigmat0Src", "Time Error from TrackExtenderWithMTD; #sigma_{t0Src} [ns]", 100, 0, 0.1);

  meTrackt0Pid_ = ibook.book1D("Trackt0Pid", "Track t0 as stored in TofPid;t0 [ns]", 100, -1, 1);
  meTrackSigmat0Pid_ = ibook.book1D("TrackSigmat0Pid", "Sigmat0 as stored in TofPid; #sigma_{t0} [ns]", 100, 0, 0.1);
  meTrackt0SafePid_ = ibook.book1D("Trackt0SafePID", "Track t0 Safe as stored in TofPid;t0 [ns]", 100, -1, 1);
  meTrackSigmat0SafePid_ =
      ibook.book1D("TrackSigmat0SafePID", "Sigmat0 Safe as stored in TofPid; #sigma_{t0} [ns]", 100, 0, 0.1);
  meTrackNumHits_ = ibook.book1D("TrackNumHits", "Number of valid MTD hits per track ; Number of hits", 10, -5, 5);
  meTrackMVAQual_ = ibook.book1D("TrackMVAQual", "Track MVA Quality as stored in Value Map ; MVAQual", 100, 0, 1);
  meTrackPathLenghtvsEta_ = ibook.bookProfile(
      "TrackPathLenghtvsEta", "MTD Track pathlength vs MTD track Eta;|#eta|;Pathlength", 100, 0, 3.2, 100.0, 400.0, "S");
  meMVATrackEffPtTot_ = ibook.book1D("MVAEffPtTot", "Pt of tracks associated to LV; track pt [GeV] ", 110, 0., 11.);
  meMVATrackMatchedEffPtTot_ =
      ibook.book1D("MVAMatchedEffPtTot", "Pt of tracks associated to LV matched to GEN; track pt [GeV] ", 110, 0., 11.);
  meMVATrackMatchedEffPtMtd_ = ibook.book1D(
      "MVAMatchedEffPtMtd", "Pt of tracks associated to LV matched to GEN with time; track pt [GeV] ", 110, 0., 11.);
  meMVATrackEffEtaTot_ = ibook.book1D("MVAEffEtaTot", "Pt of tracks associated to LV; track eta ", 66, 0., 3.3);
  meMVATrackMatchedEffEtaTot_ =
      ibook.book1D("MVAMatchedEffEtaTot", "Pt of tracks associated to LV matched to GEN; track eta ", 66, 0., 3.3);
  meMVATrackMatchedEffEtaMtd_ = ibook.book1D(
      "MVAMatchedEffEtaMtd", "Pt of tracks associated to LV matched to GEN with time; track eta ", 66, 0., 3.3);
  meMVATrackResTot_ = ibook.book1D(
      "MVATrackRes", "t_{rec} - t_{sim} for LV associated tracks; t_{rec} - t_{sim} [ns] ", 120, -0.15, 0.15);
  meMVATrackPullTot_ =
      ibook.book1D("MVATrackPull", "Pull for associated tracks; (t_{rec}-t_{sim})/#sigma_{t}", 50, -5., 5.);
  meMVATrackZposResTot_ = ibook.book1D(
      "MVATrackZposResTot", "Z_{PCA} - Z_{sim} for associated tracks;Z_{PCA} - Z_{sim} [cm] ", 100, -0.1, 0.1);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void MtdTracksValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/Tracks");
  desc.add<edm::InputTag>("inputTagG", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("inputTagT", edm::InputTag("trackExtenderWithMTD"));
  desc.add<edm::InputTag>("inputTagV", edm::InputTag("offlinePrimaryVertices4D")); //  offlinePrimaryVertices4D
  desc.add<edm::InputTag>("inputTagH", edm::InputTag("generatorSmeared"));

  //desc.add<edm::InputTag>("inputEle", edm::InputTag("gedGsfElectrons")); // Adding the electron collection name ! (At least I think I'm doing that) // barrel
  //desc.add<edm::InputTag>("inputEle", edm::InputTag("ecalDrivenGsfElectrons")); // barrel + endcap, but without track seeded electrons
  desc.add<edm::InputTag>("inputEle", edm::InputTag("ecalDrivenGsfElectronsHGC")); // only endcap electrons
  desc.add<edm::InputTag>("inputGenP", edm::InputTag("genParticles"));

  desc.add<edm::InputTag>("tmtd", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmatmtd", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("t0Src", edm::InputTag("trackExtenderWithMTD:generalTrackt0"));
  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("trackExtenderWithMTD:generalTracksigmat0"));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"))
      ->setComment("Association between General and MTD Extended tracks");
  desc.add<edm::InputTag>("pathLengthSrc", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("t0SafePID", edm::InputTag("tofPID:t0safe"));
  desc.add<edm::InputTag>("sigmat0SafePID", edm::InputTag("tofPID:sigmat0safe"));
  desc.add<edm::InputTag>("sigmat0PID", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("t0PID", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("trackMVAQual", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<double>("trackMinimumPt", 1.0);  // [GeV]
  desc.add<double>("trackMinimumEta", 1.5);
  desc.add<double>("trackMaximumEta", 3.2);
  desc.addUntracked<bool>("optionalEleIso", true);
  
  descriptions.add("mtdTracksValid", desc);
   // Added option to calculate ele iso
}

const bool MtdTracksValidation::mvaGenSel(const HepMC::GenParticle& gp, const float& charge) {
  bool match = false;
  if (gp.status() != 1) {
    return match;
  }
  match = charge != 0.f && gp.momentum().perp() > pTcut_ && std::abs(gp.momentum().eta()) < etacutGEN_;
  return match;
}

const bool MtdTracksValidation::mvaRecSel(const reco::TrackBase& trk,
                                          const reco::Vertex& vtx,
                                          const double& t0,
                                          const double& st0) {
  bool match = false;
  match = trk.pt() > pTcut_ && std::abs(trk.eta()) < etacutREC_ && std::abs(trk.vz() - vtx.z()) <= deltaZcut_;
  if (st0 > 0.) {
    match = match && std::abs(t0 - vtx.t()) < 3. * st0;
  }
  return match;
}

const bool MtdTracksValidation::mvaGenRecMatch(const HepMC::GenParticle& genP,
                                               const double& zsim,
                                               const reco::TrackBase& trk) {
  bool match = false;
  double dR = reco::deltaR(genP.momentum(), trk.momentum());
  double genPT = genP.momentum().perp();
  match =
      std::abs(genPT - trk.pt()) < trk.pt() * deltaPTcut_ && dR < deltaDRcut_ && std::abs(trk.vz() - zsim) < deltaZcut_;
  return match;
}

DEFINE_FWK_MODULE(MtdTracksValidation); 


