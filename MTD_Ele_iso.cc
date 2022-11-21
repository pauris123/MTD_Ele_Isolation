//stuff

#include <tuple>
#include "MTD_Ele_iso.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/Track.h"
// No idea why in the "example" I was given, they define another constructor like here below. Works just fine without it...
/*
MTD_Ele_iso::MTD_Ele_iso(float min_dR_cut,
                         float max_dR_cut,
                         float min_pt_cut,
                         float max_dz_cut,
                         float max_dt_cut,
                         float min_strip_cut,
                         float min_track_mtd_mva_cut,
                         edm::Handle <reco::TrackCollection> &trackCollectionH,
                         edm::ValueMap<float> &mtdt0,
                         edm::ValueMap<float> &mtdSigmat0,
                         edm::ValueMap<float> &mtdTrkQualMVA,
                         reco::Vertex &vertex)
      : OutConeSize_(max_dR_cut),
        InConeSize_(min_dR_cut),
        Ptcut_(min_pt_cut),
        dZcut_(max_dz_cut),
        dtcut_(max_dt_cut),
        StripCut_(min_strip_cut),
        MtdMvaCut_(min_track_mtd_mva_cut),
        Track_coll_(trackCollectionH),
        mtdt0_(mtdt0),
        mtdSigmat0_(mtdSigmat0),
        mtdTrkQualMVA_(mtdTrkQualMVA),
        vertex_(vertex)
        {};

*/
MTD_Ele_iso::~MTD_Ele_iso() {} //{} destructor with return???

std::tuple<int,float,float> MTD_Ele_iso::ele_iso(const reco::GsfTrack* eleTrack) const {

    int N_tracks = 0;
    double pT_sum = 0;
    double rel_pT_sum = 0;

    // Checking for primary vertex

    
    reco::Vertex Vtx_chosen; // defining event vertex (This needs to be discussed)
    std::vector <reco::Vertex> vertices = *vertexColH_; 
    int nGoodVertex = 0;

    for(int iVtx = 0; iVtx < (int) vertices.size(); iVtx++){
      const reco::Vertex &vertex = vertices.at(iVtx);
      
      bool isGoodVertex = (
          !vertex.isFake() &&
          vertex.ndof() >= 4 //&&
          //fabs(vertices.z()) <= 24.0 &&
          //fabs(vertices.position().rho()) <= 2.0
      );
      
      nGoodVertex += (int) isGoodVertex;
      
      if(nGoodVertex)
      {
        Vtx_chosen = vertex;
        break;
      }
    }


    math::XYZVector EleSigTrackMomentumAtVtx = eleTrack->momentum(); // not sure if correct, but works
    double EleSigTrackEtaAtVtx = eleTrack->eta();


    int ele_SigTrkIdx = -1; // used for IDing the matched electron track from general track collection
    int eletrkIdx = -1; 

    for (const auto& trackGen : *Track_coll_) { // looping over tracks to find the match to electron track
        
        eletrkIdx++;
        //double dr = ROOT::Math::VectorUtil::DeltaR(trackGen->momentum(), EleSigTrackMomentumAtVtx);
        double dr = reco::deltaR(trackGen.momentum(), EleSigTrackMomentumAtVtx);

        int check = 0; 

        if (dr < InConeSize_) {
          ele_SigTrkIdx = eletrkIdx;
          check++; // should do a check if we match 2 or more tracks to one electron track
        }
    }

    reco::TrackRef ele_sigTrkRef; // matched electron track ref variable
    double ele_sigTrkTime = -1;
    double ele_sigTrkTimeErr = -1;
    double ele_sigTrkMtdMva = -1;

    if (ele_SigTrkIdx >= 0) { // if we found a match, we add MTD timing information for this matched track
      

        //meEle_test_->Fill(2.0); // checking for how many electrons we found a match

        const reco::TrackRef ele_sigTrkRef(Track_coll_, ele_SigTrkIdx);
        //const reco::TrackRef ele_sigTrkRef(iEvent.getHandle(GenRecTrackToken_), ele_SigTrkIdx);
        //sigTrkRef = reco::TrackRef(trackCollectionH_, sigTrkIdx); // try this 
    
        ele_sigTrkTime = mtdt0_[ele_sigTrkRef]; // matched electron track mtd time // probably will need iEvent.get()
        ele_sigTrkTimeErr = mtdSigmat0_[ele_sigTrkRef]; // matched electron track mtd time error
        ele_sigTrkMtdMva = mtdTrkQualMVA_[ele_sigTrkRef]; // matched electron track mtd MVA score
    
        ele_sigTrkTimeErr = (ele_sigTrkMtdMva > MtdMvaCut_)? ele_sigTrkTimeErr: -1; // track MTD MVA cut/check
    }

    int general_index = 0;
    
    for (const auto& trackGen : *Track_coll_) { // looping over all the tracks to get charged rel_iso for electrons
      
        // for testing
        //meEle_test_->Fill(3.0); // checking for how many tracks we scan thorugh over all loops

        const reco::TrackRef trackref_general(Track_coll_, general_index);
        //const reco::TrackRef trackref_general(iEvent.getHandle(GenRecTrackToken_), general_index);
        //const reco::TrackRef trkRef(trackCollectionH_, trkIdx); //try this
        general_index++;

        if (trackGen.pt() < Ptcut_ ){ // track pT cut, might ignore at the first iterations of the code
          continue;
        }

        if (fabs(trackGen.vz() - eleTrack->vz()) > dZcut_ ){ // should do some dz or some other Z axis cut/check
          continue;
        }

        //meEle_test_->Fill(4.0); // how many tracks survive dz and pt cuts
        ///
        // Spot for some other track cuts we'd like to use for electron isolation //
        ///

        double dt_vtx = 0; // dt regular track vs vtx
        //double dt_vtx_signif = 0; // significance?

        double dt_sigTrk = 0; // dt regular track vs signal track (absolute value)
        //double dt_sigTrk_signif = 0;

        double TrkMTDTime = mtdt0_[trackref_general]; // MTD timing info for particular track
        double TrkMTDTimeErr = mtdSigmat0_[trackref_general];
        double TrkMTDMva = mtdTrkQualMVA_[trackref_general];

        TrkMTDTimeErr = (TrkMTDMva > MtdMvaCut_)? TrkMTDTimeErr: -1; // track MTD MVA cut/check
        
        // check if vertex is good or not? (No, the check is done in MTDValidation code)
        

        //if(isGoodVtx && TrkMTDTimeErr > 0 && primRecoVtx.tError() > 0){
        if(TrkMTDTimeErr > 0 && Vtx_chosen.tError() > 0){
           // for testing
          //meEle_test_->Fill(5.0); // how many pass dt vtx "if"


          //dt_vtx = TrkMTDTime - primRecoVtx.t();
          dt_vtx = TrkMTDTime - Vtx_chosen.t();
          //dt_vtx_signif = dt_vtx/std::sqrt(TrkMTDTimeErr*TrkMTDTimeErr + primRecoVtx.tError()*primRecoVtx.tError());

          if(dt_vtx > dtcut_){ // if dt larger then dtCut, then we skip this track
            continue;
          }

        }

        //meEle_test_->Fill(6.0); // How many pass vertex "if" and dt vtx cut

        if(TrkMTDTimeErr > 0 && ele_sigTrkTimeErr > 0){
          
          //meEle_test_->Fill(7.0); // how many pass dt signal"if", and previous if's

          dt_sigTrk = fabs(TrkMTDTime - ele_sigTrkTime);
          //dt_sigTrk_signif = dt_sigTrk/std::sqrt(TrkMTDTimeErr*TrkMTDTimeErr + ele_sigTrkTimeErr*ele_sigTrkTimeErr); 

          if(dt_sigTrk > dtcut_){ // if dt larger then dtCut, then we skip this track
            continue;
          }

        }

        //meEle_test_->Fill(8.0); // How many pass signal"if" and dt vtx cut, and previous if's

        double dr_check = reco::deltaR(trackGen.momentum(), EleSigTrackMomentumAtVtx);
        double deta = fabs(trackGen.eta() - EleSigTrackEtaAtVtx);

        if(dr_check > InConeSize_ && dr_check < OutConeSize_ && deta > StripCut_){
          ++N_tracks;
          pT_sum += trackGen.pt(); 
        }


    } // looping over tracks ///

    rel_pT_sum = pT_sum/eleTrack->pt(); //ele.pt();

std::tuple<int,float,float> iso_output;
iso_output = std::make_tuple(N_tracks,pT_sum,rel_pT_sum);
//iso_output.first = N_tracks;
//iso_output.second = pT_sum;
//iso_output.third = rel_pT_sum;

return iso_output;

}





