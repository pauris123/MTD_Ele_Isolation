#include <tuple>
#include <string>
#include <vector>

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTrackSelector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"

//reco::Vertex &vertex
//vertex_(vertex),
//reco::Vertex vertex_;

class MTD_Ele_iso {
public:
    MTD_Ele_iso(float min_dR_cut,
                float max_dR_cut,
                float min_pt_cut,
                float max_dz_cut,
                float max_dz_vtx_cut,
                float max_dxy_vtx_cut,
                std::vector<float> max_dt_cut_vtx,
                std::vector<float> max_dt_cut_track,
                float min_strip_cut,
                float min_track_mtd_mva_cut,
                edm::Handle <reco::TrackCollection> trackCollectionH, // was &trackCollectionH, couldn't get it working
                edm::ValueMap<float> mtdt0, // same as above
                edm::ValueMap<float> mtdSigmat0, // same as above
                edm::ValueMap<float> mtdTrkQualMVA, // same as above
                edm::Handle <reco::VertexCollection> vertexCol)
                //edm::Handle <reco::GenParticleCollection> GenPartCol) // same as above
      : OutConeSize_(max_dR_cut),
        InConeSize_(min_dR_cut),
        Ptcut_(min_pt_cut),
        dZcut_(max_dz_cut),
        dZvtxCut_(max_dz_vtx_cut),
        dxyVtxCut_(max_dxy_vtx_cut),
        dtcut_vtx_(max_dt_cut_vtx),
        dtcut_track_(max_dt_cut_track),
        StripCut_(min_strip_cut),
        MtdMvaCut_(min_track_mtd_mva_cut),
        Track_coll_(trackCollectionH),
        mtdt0_(mtdt0),
        mtdSigmat0_(mtdSigmat0),
        mtdTrkQualMVA_(mtdTrkQualMVA),
        vertexColH_(vertexCol)
        //GenPartColH_(GenPartCol)
        {}; // ? no idea how to end the constructor
    
    
    ~MTD_Ele_iso();

    std::tuple<int,float,float,int,float,float,int,float,float,int,float,float,int,float,float,int,float,float,int,float,float,int,float,float,int> ele_iso(const reco::GsfTrack*) const ; // const at the end???

private:

    float OutConeSize_;
    float InConeSize_;
    float Ptcut_;
    float dZcut_;
    float dZvtxCut_;
    float dxyVtxCut_;
    std::vector<float> dtcut_vtx_;
    std::vector<float> dtcut_track_;
    float StripCut_;
    float MtdMvaCut_;
    edm::Handle <reco::TrackCollection> Track_coll_;
    edm::ValueMap<float> mtdt0_;
    edm::ValueMap<float> mtdSigmat0_;
    edm::ValueMap<float> mtdTrkQualMVA_;
    edm::Handle <reco::VertexCollection> vertexColH_;
    //edm::Handle <reco::GenParticleCollection> GenPartColH_;


};
