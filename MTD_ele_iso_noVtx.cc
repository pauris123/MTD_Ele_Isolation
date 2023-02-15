/*

#include <tuple>
#include "MTD_ele_iso.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/Track.h"



MTD_Ele_iso::~MTD_Ele_iso() {} //{} destructor with return???

std::tuple<int,float,float,int,float,float,int,float,float,int,float,float,int,float,float,int,float,float,int,float,float,int,float,float> MTD_Ele_iso::ele_iso(const reco::GsfTrack* eleTrack) const {


    int N_tracks = -1; // No MTD case
    double pT_sum = -1;
    double rel_pT_sum = -1;

    int N_tracks_MTD_1 = -1; // MTD case //so we can seperate electrons that are not matched to primary vertex
    double pT_sum_MTD_1 = -1;
    double rel_pT_sum_MTD_1 = -1;

    int N_tracks_MTD_2 = -1; // MTD case //so we can seperate electrons that are not matched to primary vertex
    double pT_sum_MTD_2 = -1;
    double rel_pT_sum_MTD_2 = -1;

    int N_tracks_MTD_3 = -1; // MTD case //so we can seperate electrons that are not matched to primary vertex
    double pT_sum_MTD_3 = -1;
    double rel_pT_sum_MTD_3 = -1;

    int N_tracks_MTD_4 = -1; // MTD case //so we can seperate electrons that are not matched to primary vertex
    double pT_sum_MTD_4 = -1;
    double rel_pT_sum_MTD_4 = -1;

    int N_tracks_MTD_5 = -1; // MTD case //so we can seperate electrons that are not matched to primary vertex
    double pT_sum_MTD_5 = -1;
    double rel_pT_sum_MTD_5 = -1;

    int N_tracks_MTD_6 = -1; // MTD case //so we can seperate electrons that are not matched to primary vertex
    double pT_sum_MTD_6 = -1;
    double rel_pT_sum_MTD_6 = -1;

    int N_tracks_MTD_7 = -1; // MTD case //so we can seperate electrons that are not matched to primary vertex
    double pT_sum_MTD_7 = -1;
    double rel_pT_sum_MTD_7 = -1;

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

    // Checking if electron comes from primary vertex

    //float ele_track_source = 0;
    //ele_track_source = fabs(eleTrack->vertex().z()-Vtx_chosen.z()); // eleTrack->vertex()

    // checking only elecctrons who's tracks are matched with primary vtx < 0.1 cut
      
    for (const auto& trackGen : *Track_coll_) { // looping over tracks to find the match to electron track
          
        eletrkIdx++;
        //double dr = ROOT::Math::VectorUtil::DeltaR(trackGen->momentum(), EleSigTrackMomentumAtVtx);
        double dr = reco::deltaR(trackGen.momentum(), EleSigTrackMomentumAtVtx); // should do a pT cut aswell

        int check = 0; 

        if (dr < InConeSize_){
          ele_SigTrkIdx = eletrkIdx;
          check++; // should do a check if we match 2 or more tracks to one electron track
        }
    }

    reco::TrackRef ele_sigTrkRef; // matched electron track ref variable
    double ele_sigTrkTime = -1; // comment in when we use signal track vs track timing
    double ele_sigTrkTimeErr = -1;
    double ele_sigTrkMtdMva = -1;


    if (ele_SigTrkIdx >= 0) { // if we found a match, we add MTD timing information for this matched track and check if track is in vertex fit track collection
        
        const reco::TrackRef ele_sigTrkRef(Track_coll_, ele_SigTrkIdx);
        //const reco::TrackRef ele_sigTrkRef(iEvent.getHandle(GenRecTrackToken_), ele_SigTrkIdx);
        //sigTrkRef = reco::TrackRef(trackCollectionH_, sigTrkIdx); // try this 
      
        ele_sigTrkTime = mtdt0_[ele_sigTrkRef]; // matched electron track mtd time // probably will need iEvent.get()
        ele_sigTrkTimeErr = mtdSigmat0_[ele_sigTrkRef]; // matched electron track mtd time error
        ele_sigTrkMtdMva = mtdTrkQualMVA_[ele_sigTrkRef]; // matched electron track mtd MVA score
      
        ele_sigTrkTimeErr = (ele_sigTrkMtdMva > MtdMvaCut_)? ele_sigTrkTimeErr: -1; // track MTD MVA cut/check

        int general_index = 0;
      
        for (const auto& trackGen : *Track_coll_) { // looping over all the tracks to get charged rel_iso for electrons
        
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

            double dr_check = reco::deltaR(trackGen.momentum(), EleSigTrackMomentumAtVtx);
            double deta = fabs(trackGen.eta() - EleSigTrackEtaAtVtx);

            ele_sigTrkTimeErr = (ele_sigTrkMtdMva > MtdMvaCut_)? ele_sigTrkTimeErr: -1;

            if(dr_check > InConeSize_ && dr_check < OutConeSize_ && deta > StripCut_){
            
              if(N_tracks == -1 && pT_sum == -1){ // 
               ++N_tracks; // safety check for electrons that are not matched to primary vertex // Now all these 3 parameters should be 0 and are ready to be counted
               ++pT_sum;
               ++rel_pT_sum;
              }
            
              ++N_tracks;
              pT_sum += trackGen.pt();
            } else {
              continue;
            }

            // Here the noMTD case ends

            double dt_vtx = 0; // dt regular track vs vtx
            //double dt_vtx_signif = 0; // significance?

            double dt_sigTrk = 0; // dt regular track vs signal track (absolute value) // comment in when using signal track vs track timing
            //double dt_sigTrk_signif = 0;

            double TrkMTDTime = mtdt0_[trackref_general]; // MTD timing info for particular track
            double TrkMTDTimeErr = mtdSigmat0_[trackref_general];
            double TrkMTDMva = mtdTrkQualMVA_[trackref_general];

            TrkMTDTimeErr = (TrkMTDMva > MtdMvaCut_)? TrkMTDTimeErr: -1; // track MTD MVA cut/check
          
            
            //// Track vs track sig cut part

            if(TrkMTDTimeErr > 0 && ele_sigTrkTimeErr > 0){

              dt_sigTrk = fabs(TrkMTDTime - ele_sigTrkTime);
              //dt_sigTrk_signif = dt_sigTrk/std::sqrt(TrkMTDTimeErr*TrkMTDTimeErr + ele_sigTrkTimeErr*ele_sigTrkTimeErr); 

              if(dt_sigTrk > dtcut_track_[0]){ // if dt larger then dtCut, then we skip this track
                continue;
              }

            }

            //// Track vs vtx cut part // can change to track vs track quite easily

            if(TrkMTDTimeErr > 0 && Vtx_chosen.tError() > 0){

              dt_vtx = TrkMTDTime - Vtx_chosen.t();
              //dt_vtx_signif = dt_vtx/std::sqrt(TrkMTDTimeErr*TrkMTDTimeErr + Vtx_chosen.tError()*Vtx_chosen.tError());


              if(dt_vtx < dtcut_track_[0]){ // if dt larger then dtCut, then we skip this track
                if(N_tracks_MTD_1 == -1 && pT_sum_MTD_1 == -1){  
                ++N_tracks_MTD_1; // safety check for electrons that are not matched to primary vertex // Now all these 3 parameters should be 0 and are ready to be counted
                ++pT_sum_MTD_1;
                ++rel_pT_sum_MTD_1;
                }
                ++N_tracks_MTD_1;
                pT_sum_MTD_1 += trackGen.pt();
              }else{
                continue;
              }

              if(dt_vtx < dtcut_track_[1]){ 
                if(N_tracks_MTD_2 == -1 && pT_sum_MTD_2 == -1){  
                ++N_tracks_MTD_2; 
                ++pT_sum_MTD_2;
                ++rel_pT_sum_MTD_2;
                }
                ++N_tracks_MTD_2;
                pT_sum_MTD_2 += trackGen.pt();
              }else{
                continue;
              }

              if(dt_vtx < dtcut_track_[2]){ 
                if(N_tracks_MTD_3 == -1 && pT_sum_MTD_3 == -1){  
                ++N_tracks_MTD_3; 
                ++pT_sum_MTD_3;
                ++rel_pT_sum_MTD_3;
                }
                ++N_tracks_MTD_3;
                pT_sum_MTD_3 += trackGen.pt();
              }else{
                continue;
              }

              if(dt_vtx < dtcut_track_[3]){ 
                if(N_tracks_MTD_4 == -1 && pT_sum_MTD_4 == -1){  
                ++N_tracks_MTD_4; 
                ++pT_sum_MTD_4;
                ++rel_pT_sum_MTD_4;
                }
                ++N_tracks_MTD_4;
                pT_sum_MTD_4 += trackGen.pt();
              }else{
                continue;
              }

              if(dt_vtx < dtcut_track_[4]){ 
                if(N_tracks_MTD_5 == -1 && pT_sum_MTD_5 == -1){  
                ++N_tracks_MTD_5; 
                ++pT_sum_MTD_5;
                ++rel_pT_sum_MTD_5;
                }
                ++N_tracks_MTD_5;
                pT_sum_MTD_5 += trackGen.pt();
              }else{
                continue;
              }

              if(dt_vtx < dtcut_track_[5]){ 
                if(N_tracks_MTD_6 == -1 && pT_sum_MTD_6 == -1){  
                ++N_tracks_MTD_6; 
                ++pT_sum_MTD_6;
                ++rel_pT_sum_MTD_6;
                }
                ++N_tracks_MTD_6;
                pT_sum_MTD_6 += trackGen.pt();
              }else{
                continue;
              }

              if(dt_vtx < dtcut_track_[6]){ 
                if(N_tracks_MTD_7 == -1 && pT_sum_MTD_7 == -1){  
                ++N_tracks_MTD_7; 
                ++pT_sum_MTD_7;
                ++rel_pT_sum_MTD_7;
                }
                ++N_tracks_MTD_7;
                pT_sum_MTD_7 += trackGen.pt();
              }else{
                continue;
              }     
                       
            } 
        } // looping over tracks ///

        if(N_tracks == -1 && pT_sum == -1){
            ++N_tracks; // if we have 0 tracks in electron iso cone, we get 0 track result, not -1 tracks, same for all MTD cases
            ++pT_sum;
            ++rel_pT_sum;
        }
        if(N_tracks_MTD_1 == -1 && pT_sum_MTD_1 == -1){
            ++N_tracks_MTD_1; // if regular noMTD case has 0 tracks, than ofcourse MTD one will have 0 aswell
            ++pT_sum_MTD_1;
            ++rel_pT_sum_MTD_1;
        }
        if(N_tracks_MTD_2 == -1 && pT_sum_MTD_2 == -1){
            ++N_tracks_MTD_2; 
            ++pT_sum_MTD_2;
            ++rel_pT_sum_MTD_2;
        }
        if(N_tracks_MTD_3 == -1 && pT_sum_MTD_3 == -1){
            ++N_tracks_MTD_3; 
            ++pT_sum_MTD_3;
            ++rel_pT_sum_MTD_3;
        }
        if(N_tracks_MTD_4 == -1 && pT_sum_MTD_4 == -1){
            ++N_tracks_MTD_4; 
            ++pT_sum_MTD_4;
            ++rel_pT_sum_MTD_4;
        }
        if(N_tracks_MTD_5 == -1 && pT_sum_MTD_5 == -1){
            ++N_tracks_MTD_5; 
            ++pT_sum_MTD_5;
            ++rel_pT_sum_MTD_5;
        }
        if(N_tracks_MTD_6 == -1 && pT_sum_MTD_6 == -1){
            ++N_tracks_MTD_6; 
            ++pT_sum_MTD_6;
            ++rel_pT_sum_MTD_6;
        }
        if(N_tracks_MTD_7 == -1 && pT_sum_MTD_7 == -1){
            ++N_tracks_MTD_7; 
            ++pT_sum_MTD_7;
            ++rel_pT_sum_MTD_7;
        }

        // calculating rel_ch_iso
          
        rel_pT_sum = pT_sum/eleTrack->pt(); //ele.pt();
        rel_pT_sum_MTD_1 = pT_sum_MTD_1/eleTrack->pt(); //ele.pt();
        rel_pT_sum_MTD_2 = pT_sum_MTD_2/eleTrack->pt();
        rel_pT_sum_MTD_3 = pT_sum_MTD_3/eleTrack->pt();
        rel_pT_sum_MTD_4 = pT_sum_MTD_4/eleTrack->pt();
        rel_pT_sum_MTD_5 = pT_sum_MTD_5/eleTrack->pt();
        rel_pT_sum_MTD_6 = pT_sum_MTD_6/eleTrack->pt();
        rel_pT_sum_MTD_7 = pT_sum_MTD_7/eleTrack->pt();
        
          
    }    

std::tuple<int,float,float,int,float,float,int,float,float,int,float,float,int,float,float,int,float,float,int,float,float,int,float,float> iso_output;
iso_output = std::make_tuple(N_tracks,pT_sum,rel_pT_sum , N_tracks_MTD_1,pT_sum_MTD_1,rel_pT_sum_MTD_1,
N_tracks_MTD_2,pT_sum_MTD_2,rel_pT_sum_MTD_2 , N_tracks_MTD_3,pT_sum_MTD_3,rel_pT_sum_MTD_3 , N_tracks_MTD_4,pT_sum_MTD_4,rel_pT_sum_MTD_4,
N_tracks_MTD_5,pT_sum_MTD_5,rel_pT_sum_MTD_5 , N_tracks_MTD_6,pT_sum_MTD_6,rel_pT_sum_MTD_6 , N_tracks_MTD_7,pT_sum_MTD_7,rel_pT_sum_MTD_7);
//iso_output.first = N_tracks;
//iso_output.second = pT_sum;
//iso_output.third = rel_pT_sum;

return iso_output;

}



*/


