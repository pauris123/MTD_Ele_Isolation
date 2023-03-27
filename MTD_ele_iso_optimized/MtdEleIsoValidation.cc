/**/
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
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"

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


class MtdEleIsoValidation : public DQMEDAnalyzer {
public:
  explicit MtdEleIsoValidation(const edm::ParameterSet&);
  ~MtdEleIsoValidation() override;

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

  bool electron_iso_calc_;
  bool track_match_PV_;
  bool dt_sig_vtx_;
  bool dt_sig_track_;
  bool vertex_3D_;
  bool vertex_4D_;

  edm::EDGetTokenT<reco::TrackCollection> GenRecTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> RecTrackToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> RecVertexToken_4D_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> RecVertexToken_3D_;

  edm::EDGetTokenT<reco::GsfElectronCollection> GsfElectronToken_EB_; // Adding token for electron collection 
  edm::EDGetTokenT<reco::GsfElectronCollection> GsfElectronToken_EE_;
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

  edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particleTableToken_;

  MonitorElement* meEleISO_Ntracks_EB_; // Adding histograms for barrel electrons (isolation stuff)
  MonitorElement* meEleISO_chIso_EB_;
  MonitorElement* meEleISO_rel_chIso_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_1_EB_; 
  MonitorElement* meEleISO_chIso_MTD_1_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_1_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_2_EB_; 
  MonitorElement* meEleISO_chIso_MTD_2_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_2_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_3_EB_; 
  MonitorElement* meEleISO_chIso_MTD_3_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_3_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_4_EB_; 
  MonitorElement* meEleISO_chIso_MTD_4_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_4_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_5_EB_; 
  MonitorElement* meEleISO_chIso_MTD_5_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_5_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_6_EB_; 
  MonitorElement* meEleISO_chIso_MTD_6_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_6_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_7_EB_; 
  MonitorElement* meEleISO_chIso_MTD_7_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_7_EB_;

  MonitorElement* meEle_pt_tot_EB_;
  MonitorElement* meEle_eta_tot_EB_;
  MonitorElement* meEle_phi_tot_EB_;

  MonitorElement* meEle_pt_MTD_1_EB_;
  MonitorElement* meEle_eta_MTD_1_EB_;
  MonitorElement* meEle_phi_MTD_1_EB_;

  MonitorElement* meEle_pt_MTD_2_EB_;
  MonitorElement* meEle_eta_MTD_2_EB_;
  MonitorElement* meEle_phi_MTD_2_EB_;

  MonitorElement* meEle_pt_MTD_3_EB_;
  MonitorElement* meEle_eta_MTD_3_EB_;
  MonitorElement* meEle_phi_MTD_3_EB_;

  MonitorElement* meEle_pt_MTD_4_EB_;
  MonitorElement* meEle_eta_MTD_4_EB_;
  MonitorElement* meEle_phi_MTD_4_EB_;

  MonitorElement* meEle_pt_MTD_5_EB_;
  MonitorElement* meEle_eta_MTD_5_EB_;
  MonitorElement* meEle_phi_MTD_5_EB_;

  MonitorElement* meEle_pt_MTD_6_EB_;
  MonitorElement* meEle_eta_MTD_6_EB_;
  MonitorElement* meEle_phi_MTD_6_EB_;

  MonitorElement* meEle_pt_MTD_7_EB_;
  MonitorElement* meEle_eta_MTD_7_EB_;
  MonitorElement* meEle_phi_MTD_7_EB_;

  MonitorElement* meEle_pt_noMTD_EB_;
  MonitorElement* meEle_eta_noMTD_EB_;
  MonitorElement* meEle_phi_noMTD_EB_;


  MonitorElement* meEleISO_Ntracks_EE_; // Adding histograms for endcap electrons (isolation stuff)
  MonitorElement* meEleISO_chIso_EE_;
  MonitorElement* meEleISO_rel_chIso_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_1_EE_; 
  MonitorElement* meEleISO_chIso_MTD_1_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_1_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_2_EE_; 
  MonitorElement* meEleISO_chIso_MTD_2_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_2_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_3_EE_; 
  MonitorElement* meEleISO_chIso_MTD_3_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_3_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_4_EE_; 
  MonitorElement* meEleISO_chIso_MTD_4_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_4_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_5_EE_; 
  MonitorElement* meEleISO_chIso_MTD_5_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_5_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_6_EE_; 
  MonitorElement* meEleISO_chIso_MTD_6_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_6_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_7_EE_; 
  MonitorElement* meEleISO_chIso_MTD_7_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_7_EE_;

  MonitorElement* meEle_pt_tot_EE_;
  MonitorElement* meEle_eta_tot_EE_;
  MonitorElement* meEle_phi_tot_EE_;

  MonitorElement* meEle_pt_MTD_1_EE_;
  MonitorElement* meEle_eta_MTD_1_EE_;
  MonitorElement* meEle_phi_MTD_1_EE_;

  MonitorElement* meEle_pt_MTD_2_EE_;
  MonitorElement* meEle_eta_MTD_2_EE_;
  MonitorElement* meEle_phi_MTD_2_EE_;

  MonitorElement* meEle_pt_MTD_3_EE_;
  MonitorElement* meEle_eta_MTD_3_EE_;
  MonitorElement* meEle_phi_MTD_3_EE_;

  MonitorElement* meEle_pt_MTD_4_EE_;
  MonitorElement* meEle_eta_MTD_4_EE_;
  MonitorElement* meEle_phi_MTD_4_EE_;

  MonitorElement* meEle_pt_MTD_5_EE_;
  MonitorElement* meEle_eta_MTD_5_EE_;
  MonitorElement* meEle_phi_MTD_5_EE_;

  MonitorElement* meEle_pt_MTD_6_EE_;
  MonitorElement* meEle_eta_MTD_6_EE_;
  MonitorElement* meEle_phi_MTD_6_EE_;

  MonitorElement* meEle_pt_MTD_7_EE_;
  MonitorElement* meEle_eta_MTD_7_EE_;
  MonitorElement* meEle_phi_MTD_7_EE_;

  MonitorElement* meEle_pt_noMTD_EE_;
  MonitorElement* meEle_eta_noMTD_EE_;
  MonitorElement* meEle_phi_noMTD_EE_;

};

// ------------ constructor and destructor --------------
MtdEleIsoValidation::MtdEleIsoValidation(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      trackMinPt_(iConfig.getParameter<double>("trackMinimumPt")),
      trackMinEta_(iConfig.getParameter<double>("trackMinimumEta")),
      trackMaxEta_(iConfig.getParameter<double>("trackMaximumEta")),
      electron_iso_calc_(iConfig.getUntrackedParameter<bool>("optionalEleIso")),
      track_match_PV_(iConfig.getUntrackedParameter<bool>("optionTrackMatchToPV")),
      dt_sig_vtx_(iConfig.getUntrackedParameter<bool>("option_dtToPV")),
      dt_sig_track_(iConfig.getUntrackedParameter<bool>("option_dtToTrack")),
      vertex_3D_(iConfig.getUntrackedParameter<bool>("option_vtx_3D")),
      vertex_4D_(iConfig.getUntrackedParameter<bool>("option_vtx_4D")) {
  GenRecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagG"));
  RecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagT"));
  RecVertexToken_4D_ = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("inputTag_4D")); // Vtx 4D collection
  RecVertexToken_3D_ = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("inputTag_3D")); // Vtx 3D collection

  GsfElectronToken_EB_ = consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("inputEle_EB")); // Barrel electron collection input/token
  GsfElectronToken_EE_ = consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("inputEle_EE")); // Endcap electron collection input/token
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
  particleTableToken_ = esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>();
}

MtdEleIsoValidation::~MtdEleIsoValidation() {}

// ------------ method called for each event  ------------
void MtdEleIsoValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace geant_units::operators;
  using namespace std;

  //auto GenRecTrackHandle = makeValid(iEvent.getHandle(GenRecTrackToken_));

  edm::Handle <reco::TrackCollection> GenRecTrackHandle = iEvent.getHandle(GenRecTrackToken_);
  
  auto RecVertexHandle_4D = makeValid(iEvent.getHandle(RecVertexToken_4D_));
  auto RecVertexHandle_3D = makeValid(iEvent.getHandle(RecVertexToken_3D_));

  //const auto& tMtd = iEvent.get(tmtdToken_);
  //const auto& SigmatMtd = iEvent.get(SigmatmtdToken_);
  //const auto& t0Src = iEvent.get(t0SrcToken_);
  //const auto& Sigmat0Src = iEvent.get(Sigmat0SrcToken_);
  const auto& t0Pid = iEvent.get(t0PidToken_);
  const auto& Sigmat0Pid = iEvent.get(Sigmat0PidToken_);
  //const auto& t0Safe = iEvent.get(t0SafePidToken_);
  //const auto& Sigmat0Safe = iEvent.get(Sigmat0SafePidToken_);
  const auto& mtdQualMVA = iEvent.get(trackMVAQualToken_);
  //const auto& trackAssoc = iEvent.get(trackAssocToken_);
  //const auto& pathLength = iEvent.get(pathLengthToken_);

 
 
  if(electron_iso_calc_){
    
    float min_dR_cut = 0.01;
    float max_dR_cut = 0.3;
    float min_pt_cut_EB = 0.7; // 0.4 for endcap
    float min_pt_cut_EE = 0.4;
    float max_dz_cut_EB = 0.1; // will check 0.15/0.20/0.25/0.30 for tracks in iso cone // 0.1 for barrel, 0.2 for endcap
    float max_dz_cut_EE = 0.2;
    float max_dz_vtx_cut = 0.5;
    float max_dxy_vtx_cut = 0.2; 
    std::vector<double> max_dt_vtx_cut{0.30,0.20,0.15,0.10,0.08,0.06,0.04}; // has to be 7 values here!!! Code created for 7 dt values!! If change iso_object if number of dt values change.
    std::vector<double> max_dt_track_cut{0.30,0.20,0.15,0.10,0.08,0.06,0.04};
    float min_strip_cut = 0.01;
    float min_track_mtd_mva_cut = 0.5;
    
    
    auto eleHandle_EB = makeValid(iEvent.getHandle(GsfElectronToken_EB_));
    reco::GsfElectronCollection eleColl_EB = *(eleHandle_EB.product());

    auto eleHandle_EE = makeValid(iEvent.getHandle(GsfElectronToken_EE_));
    reco::GsfElectronCollection eleColl_EE = *(eleHandle_EE.product());

    auto GenPartHandle = makeValid(iEvent.getHandle(GenParticleToken_));
    reco::GenParticleCollection GenPartColl = *(GenPartHandle.product());

    auto VertexHandle_3D = iEvent.getHandle(RecVertexToken_3D_);
    auto VertexHandle_4D = iEvent.getHandle(RecVertexToken_4D_);

    // Creating combined electron collection

    std::vector<reco::GsfElectron> localEleCollection;

    for(const auto& ele_EB : eleColl_EB){
        if(ele_EB.isEB()){
            localEleCollection.push_back(ele_EB);
        }
    }

    for(const auto& ele_EE : eleColl_EE){
        if(ele_EE.isEE()){
            localEleCollection.push_back(ele_EE);
        }
    }  
    
    // Selecting the PV from 3D and 4D vertex collections
    reco::Vertex Vtx_chosen_3D;
    reco::Vertex Vtx_chosen_4D;
    reco::Vertex Vtx_chosen;
    std::vector <reco::Vertex> vertices_3D = *VertexHandle_3D;
    std::vector <reco::Vertex> vertices_4D = *VertexHandle_4D; 
    
    int nGoodVertex_3D = 0;
    int nGoodVertex_4D = 0;
    
    for(int iVtx = 0; iVtx < (int) vertices_3D.size(); iVtx++){
        const reco::Vertex &vertex = vertices_3D.at(iVtx);
      
        bool isGoodVertex = (!vertex.isFake() && vertex.ndof() >= 4);
      
        nGoodVertex_3D += (int) isGoodVertex;
      
        if(nGoodVertex_3D)
        {
            Vtx_chosen_3D = vertex;
            break;
        }
    }


    for(int iVtx = 0; iVtx < (int) vertices_4D.size(); iVtx++){
        const reco::Vertex &vertex = vertices_4D.at(iVtx);
      
        bool isGoodVertex = (!vertex.isFake() && vertex.ndof() >= 4);
      
        nGoodVertex_4D += (int) isGoodVertex;
      
        if(nGoodVertex_4D){
            Vtx_chosen_4D = vertex;
            break;
        }
    }

    if(vertex_3D_){
        Vtx_chosen = Vtx_chosen_3D;
    }else if(vertex_4D_){
        Vtx_chosen = Vtx_chosen_4D;
    }else{
        Vtx_chosen = Vtx_chosen_4D;
    }
    // Vertex selecton ends


    //for (const auto& ele : eleColl_EB){
    for (const auto& ele : localEleCollection){
      
        bool eleMatch_found = false;

        if( ele.pt()> 10 && fabs(ele.eta()) < 2.4 ){

            float ele_track_source_dz = 0;
            float ele_track_source_dxy = 0;

            ele_track_source_dz = fabs(ele.gsfTrack()->dz(Vtx_chosen.position()));
            ele_track_source_dxy = fabs(ele.gsfTrack()->dxy(Vtx_chosen.position()));
        
            if( ele_track_source_dz < max_dz_vtx_cut && ele_track_source_dxy < max_dxy_vtx_cut){
            
                math::XYZVector EleMomentum = ele.momentum(); 
            
                for(const auto& genParticle : GenPartColl){
                
                    if((genParticle.pdgId() == 11 || genParticle.pdgId() == -11) && genParticle.mother()->pdgId() == 23){

                        math::XYZVector GenPartMomentum = genParticle.momentum();
                        double dr_match = reco::deltaR(GenPartMomentum, EleMomentum);
                        
                        if( ((genParticle.pdgId() == 11 && ele.charge() == -1) || (genParticle.pdgId() == -11 && ele.charge() == 1)) && dr_match < 0.10 ){
                                            
                            eleMatch_found = true;
                            break;
                                
                        }
                    }        
                }
            }
        } 

        if(eleMatch_found){

            // Electron timing information check (whether electron signal track can be found in general track collection, for 99% we can.)
            int ele_SigTrkIdx = -1; // used for IDing the matched electron track from general track collection //
            int eletrkIdx = -1;

            math::XYZVector EleSigTrackMomentumAtVtx = ele.gsfTrack()->momentum(); // not sure if correct, but works
            double EleSigTrackEtaAtVtx = ele.gsfTrack()->eta();

            for (const auto& trackGen : *GenRecTrackHandle) { // looping over tracks to find the match to electron track
          
                eletrkIdx++;
                double dr = reco::deltaR(trackGen.momentum(), EleSigTrackMomentumAtVtx); // should do a pT cut aswell

                if (dr < min_dR_cut){
                    ele_SigTrkIdx = eletrkIdx;
                }
            }

            reco::TrackRef ele_sigTrkRef; // matched electron track ref variable
            double ele_sigTrkTime = -1; // electron signal track MTD information
            double ele_sigTrkTimeErr = -1;
            double ele_sigTrkMtdMva = -1;


            if (ele_SigTrkIdx >= 0) { // if we found a track match, we add MTD timing information for this matched track and do the isolation check
        
                const reco::TrackRef ele_sigTrkRef(GenRecTrackHandle, ele_SigTrkIdx);

                bool Barrel_ele = ele.isEB(); // for track pT/dz cuts (Different for EB and EE in TDR)
                
                float min_pt_cut = Barrel_ele ? min_pt_cut_EB : min_pt_cut_EE;
                float max_dz_cut = Barrel_ele ? max_dz_cut_EB : max_dz_cut_EE;

                
      
                ele_sigTrkTime = t0Pid[ele_sigTrkRef];
                ele_sigTrkTimeErr = Sigmat0Pid[ele_sigTrkRef];
                ele_sigTrkMtdMva = mtdQualMVA[ele_sigTrkRef];
                ele_sigTrkTimeErr = (ele_sigTrkMtdMva > min_track_mtd_mva_cut)? ele_sigTrkTimeErr: -1; // track MTD MVA cut/check

                if(Barrel_ele){
                    meEle_pt_tot_EB_->Fill(ele.pt()); // All selected electron information for efficiency plots later
                    meEle_eta_tot_EB_->Fill(ele.eta());
                    meEle_phi_tot_EB_->Fill(ele.phi());
                }else{
                    meEle_pt_tot_EE_->Fill(ele.pt()); // All selected electron information for efficiency plots later
                    meEle_eta_tot_EE_->Fill(ele.eta());
                    meEle_phi_tot_EE_->Fill(ele.phi());
                }
                               

                int N_tracks_noMTD = 0; // values for no MTD case
                double pT_sum_noMTD = 0;
                double rel_pT_sum_noMTD = 0;                
                std::vector<int> N_tracks_MTD{0,0,0,0,0,0,0}; // values for MTD case - 7 timing cuts
                std::vector<double> pT_sum_MTD{0,0,0,0,0,0,0};
                std::vector<double> rel_pT_sum_MTD{0,0,0,0,0,0,0};
                
                int general_index = 0;
                for(const auto& trackGen : *GenRecTrackHandle){
                    
                    const reco::TrackRef trackref_general(GenRecTrackHandle, general_index);
                    general_index++;

                    
                    if (trackGen.pt() < min_pt_cut ){ // track pT cut
                        continue;
                    }

                    if (fabs(trackGen.vz() - ele.gsfTrack()->vz()) > max_dz_cut ){ // general track vs signal track dz cut
                        continue;
                    }

                    if(track_match_PV_){
                        if(Vtx_chosen.trackWeight(trackref_general) < 0.5){ // cut for general track matching to PV, used as an extra check for 4D vertex collection !!!!!!!!!!!!!!!!!
                            continue;
                        }
                    }
                    
                    double dr_check = reco::deltaR(trackGen.momentum(), EleSigTrackMomentumAtVtx);
                    double deta = fabs(trackGen.eta() - EleSigTrackEtaAtVtx);

                    if(dr_check > min_dR_cut && dr_check < max_dR_cut && deta > min_strip_cut){ // checking if the track is inside isolation cone
                
                        ++N_tracks_noMTD;
                        pT_sum_noMTD += trackGen.pt();
                    }else{
                        continue;
                    }

                    // checking the MTD timing cuts                   
                     
                    double TrkMTDTime = t0Pid[trackref_general]; // MTD timing info for particular track fron general tracks
                    double TrkMTDTimeErr = Sigmat0Pid[trackref_general];
                    double TrkMTDMva = mtdQualMVA[trackref_general];
                    TrkMTDTimeErr = (TrkMTDMva > min_track_mtd_mva_cut)? TrkMTDTimeErr: -1; // track MTD MVA cut/check

                    if(dt_sig_track_){

                        double dt_sigTrk = 0; // dt regular track vs signal track (absolute value) 
                        
                        if(TrkMTDTimeErr > 0 && ele_sigTrkTimeErr > 0){

                            dt_sigTrk = fabs(TrkMTDTime - ele_sigTrkTime);
                            //dt_sigTrk_signif = dt_sigTrk/std::sqrt(TrkMTDTimeErr*TrkMTDTimeErr + ele_sigTrkTimeErr*ele_sigTrkTimeErr); 

                            for(long unsigned int i = 0; i < N_tracks_MTD.size() ; i++){
                                if(dt_sigTrk < max_dt_track_cut[i]){
                                  N_tracks_MTD[i] = N_tracks_MTD[i] + 1;
                                  pT_sum_MTD[i] = pT_sum_MTD[i] + trackGen.pt();
                                }
                            }
                        }else{ // if there is no error for MTD information, we count the MTD isolation case saem as noMTD
                            for(long unsigned int i = 0; i < N_tracks_MTD.size() ; i++){
                                N_tracks_MTD[i] = N_tracks_noMTD;
                                pT_sum_MTD[i] = pT_sum_noMTD;  
                            }
                        }
                    }

                    if(dt_sig_vtx_ && vertex_4D_){

                        double dt_vtx = 0;  // dt regular track vs vtx

                        if(TrkMTDTimeErr > 0 && Vtx_chosen.tError() > 0){

                            dt_vtx = TrkMTDTime - Vtx_chosen.t();
                            //dt_vtx_signif = dt_vtx/std::sqrt(TrkMTDTimeErr*TrkMTDTimeErr + Vtx_chosen.tError()*Vtx_chosen.tError());

                            for(long unsigned int i = 0; i < N_tracks_MTD.size() ; i++){
                                if(dt_vtx < max_dt_vtx_cut[i]){
                                  N_tracks_MTD[i] = N_tracks_MTD[i] + 1;
                                  pT_sum_MTD[i] = pT_sum_MTD[i] + trackGen.pt();
                                }
                            }
                        }else{
                            for(long unsigned int i = 0; i < N_tracks_MTD.size() ; i++){
                                N_tracks_MTD[i] = N_tracks_noMTD;
                                pT_sum_MTD[i] = pT_sum_noMTD;  
                            }
                        }
                    }
                }

                rel_pT_sum_noMTD = pT_sum_noMTD/ele.gsfTrack()->pt(); // rel_ch_iso calculation
                for(long unsigned int i = 0; i < N_tracks_MTD.size() ; i++){
                    rel_pT_sum_MTD[i] = pT_sum_MTD[i]/ele.gsfTrack()->pt(); 
                }

                std::vector<MonitorElement*> Ntracks_EB_list = {meEleISO_Ntracks_MTD_1_EB_,meEleISO_Ntracks_MTD_2_EB_,meEleISO_Ntracks_MTD_3_EB_,meEleISO_Ntracks_MTD_4_EB_,meEleISO_Ntracks_MTD_5_EB_,meEleISO_Ntracks_MTD_6_EB_,meEleISO_Ntracks_MTD_7_EB_};
                std::vector<MonitorElement*> ch_iso_EB_list = {meEleISO_chIso_MTD_1_EB_,meEleISO_chIso_MTD_2_EB_,meEleISO_chIso_MTD_3_EB_,meEleISO_chIso_MTD_4_EB_,meEleISO_chIso_MTD_5_EB_,meEleISO_chIso_MTD_6_EB_,meEleISO_chIso_MTD_7_EB_};
                std::vector<MonitorElement*> rel_ch_iso_EB_list = {meEleISO_rel_chIso_MTD_1_EB_,meEleISO_rel_chIso_MTD_2_EB_,meEleISO_rel_chIso_MTD_3_EB_,meEleISO_rel_chIso_MTD_4_EB_,meEleISO_rel_chIso_MTD_5_EB_,meEleISO_rel_chIso_MTD_6_EB_,meEleISO_rel_chIso_MTD_7_EB_};

                std::vector<MonitorElement*> Ntracks_EE_list = {meEleISO_Ntracks_MTD_1_EE_,meEleISO_Ntracks_MTD_2_EE_,meEleISO_Ntracks_MTD_3_EE_,meEleISO_Ntracks_MTD_4_EE_,meEleISO_Ntracks_MTD_5_EE_,meEleISO_Ntracks_MTD_6_EE_,meEleISO_Ntracks_MTD_7_EE_};
                std::vector<MonitorElement*> ch_iso_EE_list = {meEleISO_chIso_MTD_1_EE_,meEleISO_chIso_MTD_2_EE_,meEleISO_chIso_MTD_3_EE_,meEleISO_chIso_MTD_4_EE_,meEleISO_chIso_MTD_5_EE_,meEleISO_chIso_MTD_6_EE_,meEleISO_chIso_MTD_7_EE_};
                std::vector<MonitorElement*> rel_ch_iso_EE_list = {meEleISO_rel_chIso_MTD_1_EE_,meEleISO_rel_chIso_MTD_2_EE_,meEleISO_rel_chIso_MTD_3_EE_,meEleISO_rel_chIso_MTD_4_EE_,meEleISO_rel_chIso_MTD_5_EE_,meEleISO_rel_chIso_MTD_6_EE_,meEleISO_rel_chIso_MTD_7_EE_};

                std::vector<MonitorElement*> Ele_pT_MTD_EB_list = {meEle_pt_MTD_1_EB_,meEle_pt_MTD_2_EB_,meEle_pt_MTD_3_EB_,meEle_pt_MTD_4_EB_,meEle_pt_MTD_5_EB_,meEle_pt_MTD_6_EB_,meEle_pt_MTD_7_EB_};
                std::vector<MonitorElement*> Ele_eta_MTD_EB_list = {meEle_eta_MTD_1_EB_,meEle_eta_MTD_2_EB_,meEle_eta_MTD_3_EB_,meEle_eta_MTD_4_EB_,meEle_eta_MTD_5_EB_,meEle_eta_MTD_6_EB_,meEle_eta_MTD_7_EB_};
                std::vector<MonitorElement*> Ele_phi_MTD_EB_list = {meEle_phi_MTD_1_EB_,meEle_phi_MTD_2_EB_,meEle_phi_MTD_3_EB_,meEle_phi_MTD_4_EB_,meEle_phi_MTD_5_EB_,meEle_phi_MTD_6_EB_,meEle_phi_MTD_7_EB_};

                std::vector<MonitorElement*> Ele_pT_MTD_EE_list = {meEle_pt_MTD_1_EE_,meEle_pt_MTD_2_EE_,meEle_pt_MTD_3_EE_,meEle_pt_MTD_4_EE_,meEle_pt_MTD_5_EE_,meEle_pt_MTD_6_EE_,meEle_pt_MTD_7_EE_};
                std::vector<MonitorElement*> Ele_eta_MTD_EE_list = {meEle_eta_MTD_1_EE_,meEle_eta_MTD_2_EE_,meEle_eta_MTD_3_EE_,meEle_eta_MTD_4_EE_,meEle_eta_MTD_5_EE_,meEle_eta_MTD_6_EE_,meEle_eta_MTD_7_EE_};
                std::vector<MonitorElement*> Ele_phi_MTD_EE_list = {meEle_phi_MTD_1_EE_,meEle_phi_MTD_2_EE_,meEle_phi_MTD_3_EE_,meEle_phi_MTD_4_EE_,meEle_phi_MTD_5_EE_,meEle_phi_MTD_6_EE_,meEle_phi_MTD_7_EE_};

                if(Barrel_ele){

                    meEleISO_Ntracks_EB_->Fill(N_tracks_noMTD); // Filling hists for Ntraks and chIso sums for noMTD case //
                    meEleISO_chIso_EB_->Fill(pT_sum_noMTD);
                    meEleISO_rel_chIso_EB_->Fill(rel_pT_sum_noMTD);

                    for(long unsigned int j = 0 ; j < Ntracks_EB_list.size(); j++){
                        Ntracks_EB_list[j]->Fill(N_tracks_MTD[j]);
                        ch_iso_EB_list[j]->Fill(pT_sum_MTD[j]);
                        rel_ch_iso_EB_list[j]->Fill(rel_pT_sum_MTD[j]);
                    }

                    if(rel_pT_sum_noMTD < 0.08 ){
                        meEle_pt_noMTD_EB_->Fill(ele.pt());
                        meEle_eta_noMTD_EB_->Fill(ele.eta());
                        meEle_phi_noMTD_EB_->Fill(ele.phi());
                    }

                    for(long unsigned int k = 0 ; k < Ntracks_EB_list.size(); k++){
                        if(rel_pT_sum_MTD[k] < 0.08){
                            Ele_pT_MTD_EB_list[k]->Fill(ele.pt());
                            Ele_eta_MTD_EB_list[k]->Fill(ele.eta());
                            Ele_phi_MTD_EB_list[k]->Fill(ele.phi());
                        }   
                    }

                }else{ // for endcap
                    meEleISO_Ntracks_EE_->Fill(N_tracks_noMTD); // Filling hists for Ntraks and chIso sums for noMTD case //
                    meEleISO_chIso_EE_->Fill(pT_sum_noMTD);
                    meEleISO_rel_chIso_EE_->Fill(rel_pT_sum_noMTD);

                    for(long unsigned int j = 0 ; j < Ntracks_EE_list.size(); j++){
                        Ntracks_EE_list[j]->Fill(N_tracks_MTD[j]);
                        ch_iso_EE_list[j]->Fill(pT_sum_MTD[j]);
                        rel_ch_iso_EE_list[j]->Fill(rel_pT_sum_MTD[j]);
                    }

                    if(rel_pT_sum_noMTD < 0.08 ){
                        meEle_pt_noMTD_EE_->Fill(ele.pt());
                        meEle_eta_noMTD_EE_->Fill(ele.eta());
                        meEle_phi_noMTD_EE_->Fill(ele.phi());
                    }

                    for(long unsigned int k = 0 ; k < Ntracks_EE_list.size(); k++){
                        if(rel_pT_sum_MTD[k] < 0.08){
                            Ele_pT_MTD_EE_list[k]->Fill(ele.pt());
                            Ele_eta_MTD_EE_list[k]->Fill(ele.eta());
                            Ele_phi_MTD_EE_list[k]->Fill(ele.phi());
                        }   
                    }
                } 
            }
        } // genP matching cut
    } // electron collection inside single event
  } // Bool iso statement end

}

// ------------ method for histogram booking ------------
void MtdEleIsoValidation::bookHistograms(DQMStore::IBooker& ibook, edm::Run const& run, edm::EventSetup const& iSetup) {
  ibook.setCurrentFolder(folder_);

  // histogram booking

  meEleISO_Ntracks_EB_ = ibook.book1D("Ele_Iso_Ntracks_EB", "Tracks in isolation cone around electron track after basic cuts", 20, 0, 20); // hists for electrons
  meEleISO_chIso_EB_ = ibook.book1D("Ele_chIso_sum_EB", "Track pT sum in isolation cone around electron track after basic cuts", 400, 0, 20);
  meEleISO_rel_chIso_EB_ = ibook.book1D("Ele_rel_chIso_sum_EB", "Track relative pT sum in isolation cone around electron track after basic cuts", 200, 0, 4);

  meEleISO_Ntracks_MTD_1_EB_ = ibook.book1D("Ele_Iso_Ntracks_MTD_1_EB", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_1_EB_ = ibook.book1D("Ele_chIso_sum_MTD_1_EB", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_1_EB_ = ibook.book1D("Ele_rel_chIso_sum_MTD_1_EB", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_2_EB_ = ibook.book1D("Ele_Iso_Ntracks_MTD_2_EB", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_2_EB_ = ibook.book1D("Ele_chIso_sum_MTD_2_EB", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_2_EB_ = ibook.book1D("Ele_rel_chIso_sum_MTD_2_EB", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_3_EB_ = ibook.book1D("Ele_Iso_Ntracks_MTD_3_EB", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_3_EB_ = ibook.book1D("Ele_chIso_sum_MTD_3_EB", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_3_EB_ = ibook.book1D("Ele_rel_chIso_sum_MTD_3_EB", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_4_EB_ = ibook.book1D("Ele_Iso_Ntracks_MTD_4_EB", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_4_EB_ = ibook.book1D("Ele_chIso_sum_MTD_4_EB", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_4_EB_ = ibook.book1D("Ele_rel_chIso_sum_MTD_4_EB", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_5_EB_ = ibook.book1D("Ele_Iso_Ntracks_MTD_5_EB", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_5_EB_ = ibook.book1D("Ele_chIso_sum_MTD_5_EB", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_5_EB_ = ibook.book1D("Ele_rel_chIso_sum_MTD_5_EB", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_6_EB_ = ibook.book1D("Ele_Iso_Ntracks_MTD_6_EB", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_6_EB_ = ibook.book1D("Ele_chIso_sum_MTD_6_EB", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_6_EB_ = ibook.book1D("Ele_rel_chIso_sum_MTD_6_EB", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_7_EB_ = ibook.book1D("Ele_Iso_Ntracks_MTD_7_EB", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_7_EB_ = ibook.book1D("Ele_chIso_sum_MTD_7_EB", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_7_EB_ = ibook.book1D("Ele_rel_chIso_sum_MTD_7_EB", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);


  meEle_pt_tot_EB_ = ibook.book1D("Ele_pT_tot_EB", "Electron pT tot", 30, 10, 100); // hists for ele isto stuff start
  meEle_pt_noMTD_EB_ = ibook.book1D("Ele_pT_noMTD_EB", "Electron pT noMTD", 30, 10, 100);

  meEle_eta_tot_EB_ = ibook.book1D("Ele_eta_tot_EB", "Electron eta tot", 128, -3.2, 3.2);
  meEle_eta_noMTD_EB_ = ibook.book1D("Ele_eta_noMTD_EB", "Electron eta noMTD", 128, -3.2, 3.2);

  meEle_phi_tot_EB_ = ibook.book1D("Ele_phi_tot_EB", "Electron phi tot", 128, -3.2, 3.2);
  meEle_phi_noMTD_EB_ = ibook.book1D("Ele_phi_noMTD_EB", "Electron phi noMTD", 128, -3.2, 3.2);

  
  meEle_pt_MTD_1_EB_ = ibook.book1D("Ele_pT_MTD_1_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_1_EB_ = ibook.book1D("Ele_eta_MTD_1_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_1_EB_ = ibook.book1D("Ele_phi_MTD_1_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_2_EB_ = ibook.book1D("Ele_pT_MTD_2_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_2_EB_ = ibook.book1D("Ele_eta_MTD_2_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_2_EB_ = ibook.book1D("Ele_phi_MTD_2_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_3_EB_ = ibook.book1D("Ele_pT_MTD_3_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_3_EB_ = ibook.book1D("Ele_eta_MTD_3_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_3_EB_ = ibook.book1D("Ele_phi_MTD_3_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_4_EB_ = ibook.book1D("Ele_pT_MTD_4_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_4_EB_ = ibook.book1D("Ele_eta_MTD_4_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_4_EB_ = ibook.book1D("Ele_phi_MTD_4_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_5_EB_ = ibook.book1D("Ele_pT_MTD_5_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_5_EB_ = ibook.book1D("Ele_eta_MTD_5_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_5_EB_ = ibook.book1D("Ele_phi_MTD_5_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_6_EB_ = ibook.book1D("Ele_pT_MTD_6_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_6_EB_ = ibook.book1D("Ele_eta_MTD_6_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_6_EB_ = ibook.book1D("Ele_phi_MTD_6_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_7_EB_ = ibook.book1D("Ele_pT_MTD_7_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_7_EB_ = ibook.book1D("Ele_eta_MTD_7_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_7_EB_ = ibook.book1D("Ele_phi_MTD_7_EB", "Electron phi MTD", 128, -3.2, 3.2);



  meEleISO_Ntracks_EE_ = ibook.book1D("Ele_Iso_Ntracks_EE", "Tracks in isolation cone around electron track after basic cuts", 20, 0, 20); // hists for electrons
  meEleISO_chIso_EE_ = ibook.book1D("Ele_chIso_sum_EE", "Track pT sum in isolation cone around electron track after basic cuts", 400, 0, 20);
  meEleISO_rel_chIso_EE_ = ibook.book1D("Ele_rel_chIso_sum_EE", "Track relative pT sum in isolation cone around electron track after basic cuts", 200, 0, 4);

  meEleISO_Ntracks_MTD_1_EE_ = ibook.book1D("Ele_Iso_Ntracks_MTD_1_EE", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_1_EE_ = ibook.book1D("Ele_chIso_sum_MTD_1_EE", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_1_EE_ = ibook.book1D("Ele_rel_chIso_sum_MTD_1_EE", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_2_EE_ = ibook.book1D("Ele_Iso_Ntracks_MTD_2_EE", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_2_EE_ = ibook.book1D("Ele_chIso_sum_MTD_2_EE", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_2_EE_ = ibook.book1D("Ele_rel_chIso_sum_MTD_2_EE", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_3_EE_ = ibook.book1D("Ele_Iso_Ntracks_MTD_3_EE", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_3_EE_ = ibook.book1D("Ele_chIso_sum_MTD_3_EE", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_3_EE_ = ibook.book1D("Ele_rel_chIso_sum_MTD_3_EE", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_4_EE_ = ibook.book1D("Ele_Iso_Ntracks_MTD_4_EE", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_4_EE_ = ibook.book1D("Ele_chIso_sum_MTD_4_EE", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_4_EE_ = ibook.book1D("Ele_rel_chIso_sum_MTD_4_EE", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_5_EE_ = ibook.book1D("Ele_Iso_Ntracks_MTD_5_EE", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_5_EE_ = ibook.book1D("Ele_chIso_sum_MTD_5_EE", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_5_EE_ = ibook.book1D("Ele_rel_chIso_sum_MTD_5_EE", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_6_EE_ = ibook.book1D("Ele_Iso_Ntracks_MTD_6_EE", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_6_EE_ = ibook.book1D("Ele_chIso_sum_MTD_6_EE", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_6_EE_ = ibook.book1D("Ele_rel_chIso_sum_MTD_6_EE", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);

  meEleISO_Ntracks_MTD_7_EE_ = ibook.book1D("Ele_Iso_Ntracks_MTD_7_EE", "Tracks in isolation cone around electron track after basic cuts with MTD", 20, 0, 20); // hists for electrons
  meEleISO_chIso_MTD_7_EE_ = ibook.book1D("Ele_chIso_sum_MTD_7_EE", "Track pT sum in isolation cone around electron track after basic cuts with MTD", 400, 0, 20);
  meEleISO_rel_chIso_MTD_7_EE_ = ibook.book1D("Ele_rel_chIso_sum_MTD_7_EE", "Track relative pT sum in isolation cone around electron track after basic cuts with MTD", 200, 0, 4);


  meEle_pt_tot_EE_ = ibook.book1D("Ele_pT_tot_EE", "Electron pT tot", 30, 10, 100); // hists for ele isto stuff start
  meEle_pt_noMTD_EE_ = ibook.book1D("Ele_pT_noMTD_EE", "Electron pT noMTD", 30, 10, 100);

  meEle_eta_tot_EE_ = ibook.book1D("Ele_eta_tot_EE", "Electron eta tot", 128, -3.2, 3.2);
  meEle_eta_noMTD_EE_ = ibook.book1D("Ele_eta_noMTD_EE", "Electron eta noMTD", 128, -3.2, 3.2);

  meEle_phi_tot_EE_ = ibook.book1D("Ele_phi_tot_EE", "Electron phi tot", 128, -3.2, 3.2);
  meEle_phi_noMTD_EE_ = ibook.book1D("Ele_phi_noMTD_EE", "Electron phi noMTD", 128, -3.2, 3.2);

  
  meEle_pt_MTD_1_EE_ = ibook.book1D("Ele_pT_MTD_1_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_1_EE_ = ibook.book1D("Ele_eta_MTD_1_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_1_EE_ = ibook.book1D("Ele_phi_MTD_1_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_2_EE_ = ibook.book1D("Ele_pT_MTD_2_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_2_EE_ = ibook.book1D("Ele_eta_MTD_2_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_2_EE_ = ibook.book1D("Ele_phi_MTD_2_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_3_EE_ = ibook.book1D("Ele_pT_MTD_3_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_3_EE_ = ibook.book1D("Ele_eta_MTD_3_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_3_EE_ = ibook.book1D("Ele_phi_MTD_3_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_4_EE_ = ibook.book1D("Ele_pT_MTD_4_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_4_EE_ = ibook.book1D("Ele_eta_MTD_4_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_4_EE_ = ibook.book1D("Ele_phi_MTD_4_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_5_EE_ = ibook.book1D("Ele_pT_MTD_5_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_5_EE_ = ibook.book1D("Ele_eta_MTD_5_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_5_EE_ = ibook.book1D("Ele_phi_MTD_5_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_6_EE_ = ibook.book1D("Ele_pT_MTD_6_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_6_EE_ = ibook.book1D("Ele_eta_MTD_6_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_6_EE_ = ibook.book1D("Ele_phi_MTD_6_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_7_EE_ = ibook.book1D("Ele_pT_MTD_7_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_7_EE_ = ibook.book1D("Ele_eta_MTD_7_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_7_EE_ = ibook.book1D("Ele_phi_MTD_7_EE", "Electron phi MTD", 128, -3.2, 3.2);
  

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void MtdEleIsoValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/Tracks");
  desc.add<edm::InputTag>("inputTagG", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("inputTagT", edm::InputTag("trackExtenderWithMTD"));
  desc.add<edm::InputTag>("inputTag_4D", edm::InputTag("offlinePrimaryVertices4D")); //  offlinePrimaryVertices4D
  desc.add<edm::InputTag>("inputTag_3D", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("inputTagH", edm::InputTag("generatorSmeared"));

  desc.add<edm::InputTag>("inputEle_EB", edm::InputTag("gedGsfElectrons")); // Adding the electron collection name ! (At least I think I'm doing that) // barrel ecal and track driven electrons
  //desc.add<edm::InputTag>("inputEle", edm::InputTag("ecalDrivenGsfElectrons")); // barrel + endcap, but without track seeded electrons
  desc.add<edm::InputTag>("inputEle_EE", edm::InputTag("ecalDrivenGsfElectronsHGC")); // only endcap electrons
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
  desc.addUntracked<bool>("optionTrackMatchToPV", false);
  desc.addUntracked<bool>("option_dtToPV", true);
  desc.addUntracked<bool>("option_dtToTrack", false);
  desc.addUntracked<bool>("option_vtx_3D", false);
  desc.addUntracked<bool>("option_vtx_4D", true);
  
  descriptions.add("mtdEleIsoValid", desc);
   // Added option to calculate ele iso
}


DEFINE_FWK_MODULE(MtdEleIsoValidation); 

//*/

