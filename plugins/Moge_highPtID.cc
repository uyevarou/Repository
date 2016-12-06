// -*- C++ -*-
//
// Package:    DelauIdeal/Moge_highPtID
// Class:      Moge_highPtID
// 
/**\class Moge_highPtID Moge_highPtID.cc DelauIdeal/Moge_highPtID/plugins/Moge_highPtID.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Uladzislava Yevarouskaya
//         Created:  Tue, 15 Nov 2016 10:59:21 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <string>
#include <set>
#include <stdio.h>
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFile.h"
#include <math.h>
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "Math/LorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include <algorithm>
#include <iostream>
#include "TROOT.h"
#include "TSystem.h"
#include "TH1I.h"
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"    
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/Common/interface/RefCore.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"

#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TH1.h"
//#include "MSETools.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;
//
// class declaration
//

class Moge_highPtID : public edm::EDAnalyzer {
   public:
      explicit Moge_highPtID(const edm::ParameterSet&);
      ~Moge_highPtID();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      bool isPassHLT();
      bool RecoHLTMuonMatching(float RecoEta,float RecoPhi);
      float delR(float eta1,float phi1,float eta2,float phi2);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      void clearTree();
      // ----------member data ---------------------------
	TTree* mytree;

	      // Tokens
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::TauCollection> tauToken_;
      edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      edm::EDGetTokenT<pat::PackedGenParticleCollection> genPackedCollToken_;
      edm::EDGetTokenT<reco::GenParticleCollection> genPrunedCollToken_;
      edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
      edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
      edm::EDGetTokenT<vector<PileupSummaryInfo> > puCollection_;
      edm::InputTag mSrcRho;
      // structures
      struct PATRONUS{
       TLorentzVector p4;
       float Spt;
       float Seta;
       float Sphi;
       int Scharge;
       reco::TransientTrack Strack;
       float dxy;
       float ptError;
       int globalMuonHits;
       int pixelHits;
       int trackerLayers;
       int matchedMuStations;
       float Spx;
       float Spy;
       float Spz;
       float Sp;
       };
	static bool lepSortingRule(PATRONUS x, PATRONUS y) {return x.p4.Pt() > y.p4.Pt();}
      struct GENPARTICLE {
         // ---- momentum 4-vector ---------------------------------------
         TLorentzVector p4;
         // ---- pdgid ---------------------------------------------------
         int pdgId;
         // ---- motherid ---------------------------------------------------
         int motherId;
         int mmotherId;
         // ---- Charge ----------------------------
         int charge;
      };
      static bool lepSortingRuleGen(GENPARTICLE x, GENPARTICLE y) {return x.p4.Pt() > y.p4.Pt();}
  /*|-------------------------------------------------------|
    |                Histogramms Part                       |
    |-------------------------------------------------------|*/
      TH1D  *hWEvents_;
      TH1I *h_NbVtx;
      TH1I *h_NbGoodvtx;
      TH1F *hPU_;
      TH1F *hPUTrue_; 
      TH1F *hPUTrueTrue_;
      TH1F *DiMuon_Vrt;
      TH1F *DiMuon_Mass_fromVrt;
      TH1F *DYMass_hist;
      TH1F *DYZpkM_hist;
      TH1F *DYPt_hist;
      TH1F *DYPt_lead;
      TH1F *DYPt_sub;
      TH1F *DYEta_hist;
      TH1F *DYEta_lead;
      TH1F *DYEta_sub;
      TH1F *DYPhi_hist;
      TH1F *DYPhi_lead;
      TH1F *DYPhi_sub;
      TH1F *DYGenMass;
      TH1F *DYGenMass_nom;
      TH1F *DYGenMass_den;

  /*|-------------------------------------------------------|
    |                Run Info Part                           |
    |-------------------------------------------------------|*/     
      int Run;
      int Event;
      int lumi;
      int bunch;
      int isRealData;
  /*|-------------------------------------------------------|
    |                Trigger Part                           |
    |-------------------------------------------------------|*/
  	std::vector<int>    HLT_nb;
  	std::vector<string> HLT_name;
  	std::vector<bool>   HLT_isaccept;
  	std::vector<int>    HLTObj_nbObj;
  	std::vector<int>    HLTObj_nbObjBefore;
  	std::vector<float>  HLTObj_pt;
  	std::vector<float>  HLTObj_eta;
  	std::vector<float>  HLTObj_phi;
  	std::vector<string> HLTObj_collection;
  	std::vector<bool>   HLT_match;
  /*|-------------------------------------------------------|
    |                Gen Part                           |
    |-------------------------------------------------------|*/
      float  mcWeight_;
      float GenMass;
      float GenMass_den; 
      float GenMass_nom; 

  /*|-------------------------------------------------------|
    |                Reco Part                           |
    |-------------------------------------------------------|*/
      float  mmM_Sel;
      float  mmM_Vrt_Good;
      float  mmM_Vrt_Sel; 
      double Angle_Good;
      double Angle_Sel;
      float  vtxNormChi_Good;
      float  minVrt_Sel;
      std::vector<float>  muPt;
      std::vector<float>  muPt_Good;
      std::vector<float>  muPt_Sel;
      std::vector<float>  muEta;
      std::vector<float>  muEta_Good;
      std::vector<float>  muEta_Sel;
      std::vector<float>  muPhi;
      std::vector<float>  muPhi_Good;
      std::vector<float>  muPhi_Sel;
      std::vector<int>    muChrg;
      std::vector<int>    muChrg_Good;
      std::vector<int>    muChrg_Sel;

  /*|-------------------------------------------------------|
    |                PileUp part                            |
    |-------------------------------------------------------|*/
	Int_t                 nPUInfo_;
	std::vector<int>      nPU_;
	std::vector<int>      puBX_;
	std::vector<float>    puTrue_;
  /*|-------------------------------------------------------|
    |                PrimaryVertex Part                     |
    |-------------------------------------------------------|*/
	std::vector<int>    nbPv;
  	std::vector<int>    Nbdof;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Moge_highPtID::Moge_highPtID(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
   tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
   photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
   jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
   fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
   metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
   triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
   triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
   triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),   
   genPackedCollToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGenParticles"))),   
   genPrunedCollToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
   genJetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets")))
{
   mSrcRho            = iConfig.getParameter<edm::InputTag>             ("srcRho");
   puCollection_      = consumes<vector<PileupSummaryInfo> >    (iConfig.getParameter<InputTag>("pileupCollection"));
   conversionsToken_  = mayConsume< reco::ConversionCollection >(iConfig.getParameter<edm::InputTag>("conversions"));
   beamSpotToken_     = consumes<reco::BeamSpot> (iConfig.getParameter <edm::InputTag>("beamSpot"));

}


Moge_highPtID::~Moge_highPtID()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Moge_highPtID::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  clearTree();

  Run   = iEvent.id().run();
  Event = iEvent.id().event();
  lumi  = iEvent.luminosityBlock();
  bunch = iEvent.bunchCrossing();
  isRealData = iEvent.isRealData() ? 1:0;
  /*|-------------------------------------------------------|
    |                PrimaryVertex Part                     |
    |-------------------------------------------------------|*/
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  if (vertices->empty()) return;
  VertexCollection::const_iterator PV = vertices->end();
  for (VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++PV) {
      if ( !(vtx->isFake())
          && vtx->ndof() >= 4. && vtx->position().Rho() < 2.0
          && fabs(vtx->position().Z()) < 24.0) {
          PV = vtx;
          break;
      }
  }
  if (PV==vertices->end()) return; 
  // count how many good vertices we have
  int nGoodVtxs = 0;
  int nVtxs = 0;
  for (VertexCollection::const_iterator vtx = vertices->begin();vtx != vertices->end(); ++vtx) {
      nVtxs++;
      if ( !(vtx->isFake()) && vtx->ndof() >= 4. 
            && vtx->position().Rho() <= 2.0 
            && fabs(vtx->position().Z()) <= 24.) nGoodVtxs++; 
      h_NbVtx->Fill(nVtxs);
      h_NbGoodvtx->Fill(nGoodVtxs);
      nbPv.push_back(nGoodVtxs);
      Nbdof.push_back(vtx->ndof());
   }
  /*|-------------------------------------------------------|
    |                PileUp part                            |
    |-------------------------------------------------------|*/
  edm::Handle<vector<PileupSummaryInfo> > genPileupHandle;
  iEvent.getByToken(puCollection_, genPileupHandle);
  if (genPileupHandle.isValid()) {
    for (vector<PileupSummaryInfo>::const_iterator pu = genPileupHandle->begin(); pu != genPileupHandle->end(); ++pu) {
      if (pu->getBunchCrossing() == 0) {
        hPU_->Fill(pu->getPU_NumInteractions());
        hPUTrue_->Fill(pu->getTrueNumInteractions());
      }
      nPU_   .push_back(pu->getPU_NumInteractions());
      puTrue_.push_back(pu->getTrueNumInteractions());
      hPUTrueTrue_->Fill(pu->getTrueNumInteractions());
      puBX_  .push_back(pu->getBunchCrossing());

      nPUInfo_++;
    }
  } 
  /*|-------------------------------------------------------|
    |                Trigger Part                           |
    |-------------------------------------------------------|*/
  int NbTriggers = 0;
  int NbTriggerObj = 0;
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::TriggerObjectStandAloneCollection> trigobj_handle;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  std::string full_name;
  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerObjects_, trigobj_handle);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    std::string const& name = names.triggerName(i);
    full_name = name;
    NbTriggers++;
    HLT_nb.push_back(NbTriggers);
    HLT_name.push_back(names.triggerName(i));   
    HLT_isaccept.push_back(triggerBits->accept(i));
  }
  for (unsigned i = 0; i < trigobj_handle->size(); ++i) {
    pat::TriggerObjectStandAlone src = trigobj_handle->at(i);
    src.unpackPathNames(names);
    std::vector<std::string> const& pathnames = src.pathNames();
    for (unsigned j = 0; j < pathnames.size(); ++j) {
      if( src.pt() < 10.0 ) continue;
      if( fabs(src.eta()) > 2.5 ) continue;
      NbTriggerObj++;
      HLTObj_nbObj.push_back(NbTriggerObj);
      HLTObj_pt.push_back(src.pt()); 
      HLTObj_eta.push_back(src.eta());
      HLTObj_phi.push_back(src.phi());
      HLTObj_collection.push_back(pathnames[j]);
      //std::cout << "triggernames" <<pathnames[j]<<std::endl; 
    }
  }
  /*|-------------------------------------------------------|
    |                Gen Part                               |
    |-------------------------------------------------------|*/
  vector<GENPARTICLE> myGenLeptons;
  edm::Handle<pat::PackedGenParticleCollection> genPackedColl;
  edm::Handle<reco::GenParticleCollection> genPrunedColl;
  if (!isRealData) {
     iEvent.getByToken(genPackedCollToken_, genPackedColl);
     iEvent.getByToken(genPrunedCollToken_, genPrunedColl);

     Handle<GenEventInfoProduct> geninfo;
     iEvent.getByLabel("generator",geninfo);
     mcWeight_ = geninfo->weight();
     if (mcWeight_>0){
        hWEvents_->Fill(0.5,mcWeight_);
        hWEvents_->Fill(1.5,mcWeight_*mcWeight_);
     }
     if (mcWeight_<0) {
         hWEvents_->Fill(2.5,-1*mcWeight_);
         hWEvents_->Fill(3.5,mcWeight_*mcWeight_);
      }
      hWEvents_->Fill(4.5,mcWeight_);
      hWEvents_->Fill(5.5,1);
  
     // Find the generated Z

     for (const reco::GenParticle &i_gen : *genPrunedColl) {           
        // ----------- Muons and Electrons -------------
        if ( abs(i_gen.pdgId()) != 13 && abs(i_gen.pdgId()) != 11 ) continue; 
        // require stable particles  
        if ( i_gen.status() != 1 ) continue; 
        // set lepton 4-vector and ID
        GENPARTICLE aGenLepton;
        TLorentzVector lepP4GEN(i_gen.p4().Px(),i_gen.p4().Py(),i_gen.p4().Pz(),i_gen.p4().E());
        aGenLepton.p4 = lepP4GEN;
        aGenLepton.pdgId = i_gen.pdgId();
        // Mother particles are of type reco::GenParticle
        const reco::Candidate * motherCand = i_gen.mother(0) ;
        aGenLepton.motherId = motherCand->pdgId();
           if (fabs(aGenLepton.pdgId)==13){
               if (aGenLepton.motherId == aGenLepton.pdgId) {
	           aGenLepton.motherId = motherCand->mother()->pdgId();
	           if (aGenLepton.motherId == aGenLepton.pdgId) {
                      aGenLepton.motherId = motherCand->mother()->mother()->pdgId();
	              if (aGenLepton.motherId == aGenLepton.pdgId) {
	     	         aGenLepton.motherId = motherCand->mother()->mother()->mother()->pdgId();
	     	       }
	            }
	       }
           }
           aGenLepton.charge = i_gen.charge();
           // Keep leptons only within kinematic and geometric acceptance
           // if ( (aGenLepton.p4.Pt() < 22) && (fabs(aGenLepton.p4.Eta()) > 2.4) ) continue; 
           myGenLeptons.push_back(aGenLepton);
     } // end PackedGenParticles loop

    sort(myGenLeptons.begin(), myGenLeptons.end(), lepSortingRuleGen);
    TLorentzVector llP4GEN(0,0,0,0);
    if (myGenLeptons.size() > 1) llP4GEN = myGenLeptons[0].p4 + myGenLeptons[1].p4;
    GenMass = llP4GEN.M();
    DYGenMass -> Fill(GenMass); 
    if (myGenLeptons.size() > 1
        && myGenLeptons[0].pdgId * myGenLeptons[1].pdgId == -169
        && myGenLeptons[0].charge*myGenLeptons[1].charge<0
        && ((myGenLeptons[0].motherId == 23 && myGenLeptons[1].motherId == 23) || (myGenLeptons[0].motherId == 22 && myGenLeptons[1].motherId == 22))){
        TLorentzVector llP4GEN_den(0,0,0,0);
        llP4GEN_den = myGenLeptons[0].p4 + myGenLeptons[1].p4;
        GenMass_den = llP4GEN_den.M();
        DYGenMass_den -> Fill(GenMass_den); 
        if (fabs(myGenLeptons[0].p4.Eta() < 2.4) && fabs(myGenLeptons[1].p4.Eta() < 2.4)
           && myGenLeptons[0].p4.Pt() > 22 && myGenLeptons[1].p4.Pt() > 10){
           TLorentzVector llP4GEN_nom(0,0,0,0);
           llP4GEN_nom = myGenLeptons[0].p4 + myGenLeptons[1].p4;
           GenMass_nom = llP4GEN_nom.M();
           DYGenMass_nom -> Fill(GenMass_nom); 
        }
     }   
  }

  /*|-------------------------------------------------------|
    |                Muon Part                              |
    |-------------------------------------------------------|*/
// if you want to change on High Pt Id, use tuneP algoritm!!
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId and read for Run2

  int gdInx_mu1 = -999;
  int gdInx_mu2 = -999;
  vector<PATRONUS> GoodRecoMu;
  const reco::Vertex &PV_na = vertices->front();

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  for (const pat::Muon &mu : *muons) {
    PATRONUS muON;
    TLorentzVector muONP4(mu.p4().Px(), mu.p4().Py(), mu.p4().Pz(), mu.p4().E() );
    muON.p4      = muONP4;
    const reco::TrackRef& tunePTrack  = mu.tunePMuonBestTrack();
    reco::TrackRef glbTrack      = mu.globalTrack();
    reco::TrackRef trackerTrack  = mu.innerTrack();
    reco::TrackRef muonTrack     = mu.outerTrack();
    reco::TrackRef MuonBestTrack = mu.muonBestTrack();

    muON.matchedMuStations = mu.numberOfMatchedStations();

    if (tunePTrack.isNonnull()){
    muON.Spt = tunePTrack->pt();
    muON.Seta = tunePTrack->eta();
    muON.Sphi = tunePTrack->phi();
    muON.Scharge = tunePTrack->charge();
    muON.Spx = tunePTrack->px();
    muON.Spy = tunePTrack->py();
    muON.Spz = tunePTrack->pz();
    muON.Sp = tunePTrack->p();
    //muON.dxy = fabs(tunePTrack->dxy(PV_na.position()));
    //muON.ptError  = tunePTrack->ptError()/tunePTrack->pt();
    muPt.push_back(muON.Spt);
    muEta.push_back(muON.Seta);
    muPhi.push_back(muON.Sphi);
    muChrg.push_back(muON.Scharge);
    }
    if (mu.globalTrack().isNonnull()) {
      const reco::HitPattern & glbhit = glbTrack->hitPattern();
      muON.globalMuonHits  = glbhit.numberOfValidMuonHits();

      if (trackerTrack.isNonnull()){
         const reco::HitPattern & inhit = trackerTrack->hitPattern();
         muON.pixelHits     = inhit.numberOfValidPixelHits();
         muON.trackerLayers = inhit.trackerLayersWithMeasurement();
      }

    }
    else {
        if (trackerTrack.isNonnull()){
            const reco::HitPattern & inhit = trackerTrack->hitPattern();
            muON.pixelHits    = inhit.numberOfValidPixelHits();
            muON.trackerLayers = inhit.trackerLayersWithMeasurement();
            }
     }

    muON.dxy = fabs(mu.muonBestTrack()->dxy(PV_na.position()));
    muON.ptError  = mu.muonBestTrack()->ptError()/mu.muonBestTrack()->pt();

    if( !(muON.Spt > 10 && abs(muON.Seta) < 2.4)) continue;
    if( !mu.isGlobalMuon()) continue;
    if( !(muON.globalMuonHits > 0)) continue;
    if( !(mu.numberOfMatchedStations() > 1)) continue;
    if( !(muON.ptError < 0.3)) continue;
    if( !(muON.dxy < 0.2)) continue;
    if( !(muON.pixelHits > 0)) continue;
    if( !(muON.trackerLayers > 5)) continue;  
    float muIsoTrk = mu.isolationR03().sumPt/mu.pt();
    if (muIsoTrk >= 0.1) continue;
    bool fireHLT = isPassHLT();
    if(fireHLT == 0) continue;

    edm::ESHandle<TransientTrackBuilder> ttkb;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb);
    reco::TransientTrack ttv;
    ttv = ttkb->build(tunePTrack);
    muON.Strack  = ttv;

    GoodRecoMu.push_back(muON);
  }
  //std::cout << "size =" << GoodRecoMu.size() <<std::endl;
  sort(GoodRecoMu.begin(), GoodRecoMu.end(), lepSortingRule);
  
  for (unsigned j = 0; j < GoodRecoMu.size(); j++) {
      
      muPt_Good.push_back(GoodRecoMu[j].Spt);
      muEta_Good.push_back(GoodRecoMu[j].Seta);
      muPhi_Good.push_back(GoodRecoMu[j].Sphi);
      muChrg_Good.push_back(GoodRecoMu[j].Scharge);
      bool HLT_mu = RecoHLTMuonMatching(GoodRecoMu[j].Seta,GoodRecoMu[j].Sphi);
      HLT_match.push_back(HLT_mu);

  }
  
  if (GoodRecoMu.size() > 1){

   for(unsigned led = 0; led < GoodRecoMu.size(); led++) {
      if (!(GoodRecoMu[led].Spt > 22)) continue;

      for(unsigned sub = 0; sub < GoodRecoMu.size(); sub++) {
         if (led == sub) continue;
         // Angle
         TVector3 Mu1(GoodRecoMu[led].Spx,GoodRecoMu[led].Spy, GoodRecoMu[led].Spz);
         TVector3 Mu2(GoodRecoMu[sub].Spx, GoodRecoMu[sub].Spy, GoodRecoMu[sub].Spz);
         double cos_angle = Mu1.Dot(Mu2) / GoodRecoMu[led].Sp / GoodRecoMu[sub].Sp;
         Angle_Good = cos_angle;
         //trigger
         bool HLT_mu1 = RecoHLTMuonMatching(GoodRecoMu[led].Seta,GoodRecoMu[led].Sphi);
         bool HLT_mu2 = RecoHLTMuonMatching(GoodRecoMu[sub].Seta,GoodRecoMu[sub].Sphi);
         //vertex
         //for HightPt Muons tunePBestTrack
         std::vector<reco::TransientTrack> ttV;
         ttV.clear();
         ttV.push_back(GoodRecoMu[led].Strack);
         ttV.push_back(GoodRecoMu[sub].Strack);
         KalmanVertexFitter kvf(true);
         CachingVertex<5> vtx = kvf.vertex(ttV);
         float vtxNormChi = vtx.totalChiSquared()/vtx.degreesOfFreedom();
         //std::cout<<"vtxNormChi =" << vtxNormChi << std::endl;
         InvariantMassFromVertex imfv;
         static const double muon_mass = 0.1056583;
         vtxNormChi_Good = vtxNormChi;
         Measurement1D massOK = imfv.invariantMass(vtx, muon_mass);
         mmM_Vrt_Good = massOK.value();

         if (cos_angle < -0.999987) continue;
         if(HLT_mu1 == 0 || HLT_mu2 == 0) continue;
         if (!(GoodRecoMu[sub].Spt > 10)) continue;
         if (!(GoodRecoMu[led].Scharge * GoodRecoMu[sub].Scharge < 0)) continue;
         if (!(vtxNormChi < 20.0)) continue;
         if (vtxNormChi < minVrt_Sel) {
             minVrt_Sel = vtxNormChi;
             gdInx_mu1 = led;
             gdInx_mu2 = sub;
             Measurement1D massO = imfv.invariantMass(vtx, muon_mass);
             mmM_Vrt_Sel = massO.value();
             Angle_Sel = cos_angle;
          }  
       }
     }
  }
   if (minVrt_Sel< 20) DiMuon_Vrt -> Fill(minVrt_Sel);
   DiMuon_Mass_fromVrt -> Fill(mmM_Vrt_Sel);

   if (gdInx_mu1 > -1 && gdInx_mu2 > -1){
   muPt_Sel.push_back(GoodRecoMu[gdInx_mu1].Spt);
   muPt_Sel.push_back(GoodRecoMu[gdInx_mu2].Spt);
   DYPt_hist     -> Fill(GoodRecoMu[gdInx_mu1].Spt);
   DYPt_hist     -> Fill(GoodRecoMu[gdInx_mu2].Spt);
   DYPt_lead     -> Fill(GoodRecoMu[gdInx_mu1].Spt);
   DYPt_sub      -> Fill(GoodRecoMu[gdInx_mu2].Spt);

   muEta_Sel.push_back(GoodRecoMu[gdInx_mu1].Seta);
   muEta_Sel.push_back(GoodRecoMu[gdInx_mu2].Seta);
   DYEta_hist    -> Fill(GoodRecoMu[gdInx_mu1].Seta);
   DYEta_hist    -> Fill(GoodRecoMu[gdInx_mu2].Seta);
   DYEta_lead    -> Fill(GoodRecoMu[gdInx_mu1].Seta);
   DYEta_sub     -> Fill(GoodRecoMu[gdInx_mu2].Seta);

   muPhi_Sel.push_back(GoodRecoMu[gdInx_mu1].Sphi);
   muPhi_Sel.push_back(GoodRecoMu[gdInx_mu2].Sphi);
   DYPhi_hist    -> Fill(GoodRecoMu[gdInx_mu1].Sphi);
   DYPhi_hist    -> Fill(GoodRecoMu[gdInx_mu2].Sphi);
   DYPhi_lead    -> Fill(GoodRecoMu[gdInx_mu1].Sphi);
   DYPhi_sub     -> Fill(GoodRecoMu[gdInx_mu2].Sphi);

   muChrg_Sel.push_back(GoodRecoMu[gdInx_mu1].Scharge);
   muChrg_Sel.push_back(GoodRecoMu[gdInx_mu2].Scharge);

   TLorentzVector mu1, mu2;
   mu1.SetPtEtaPhiM(GoodRecoMu[gdInx_mu1].Spt,GoodRecoMu[gdInx_mu1].Seta,GoodRecoMu[gdInx_mu1].Sphi, 0.1056583);
   mu2.SetPtEtaPhiM(GoodRecoMu[gdInx_mu2].Spt,GoodRecoMu[gdInx_mu2].Seta,GoodRecoMu[gdInx_mu2].Sphi, 0.1056583);
   mmM_Sel = (mu1+mu2).M();
   DYMass_hist -> Fill(mmM_Sel);
   if (mmM_Sel > 50 && mmM_Sel < 200) DYZpkM_hist -> Fill(mmM_Sel);

  }
  mytree->Fill();
}
// ------------ method called once each job just before starting event loop  ------------
void 
Moge_highPtID::beginJob()
{
    Float_t Array[46] = {15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,
    64.,68.,72.,76.,81.,86.,91.,96.,101.,106.,110.,115.,
    120.,126.,133.,141.,150.,160.,171.,185.,
    200.,220.,243.,273.,320.,380.,440.,510.,
    600.,700.,830.,1000.,1200.,1500.,2000.,3000.};
    edm::Service<TFileService> fs;
    mytree = fs->make<TTree>("tree","tr");
//------------------ Weight------------------------------------------------------
    hWEvents_   = fs->make<TH1D>("WEvents", "Weighted Events",6,0,6);
    hWEvents_->Sumw2();
//-------------------Pile up-----------------------------------------------
    hPU_            = fs->make<TH1F>("hPU",    "number of pileup",       200,  0, 200);
    hPUTrue_        = fs->make<TH1F>("hPUTrue"," nb with BeamCross",        50, 0, 50);
    hPUTrueTrue_    = fs->make<TH1F>("hPUTrueTrue","number of true pilepu", 50, 0, 50);
//------------------ Vertex------------------------------------------------------
    h_NbVtx       = fs->make<TH1I>("NVtx","NVtx",50,0,50);
    h_NbGoodvtx   = fs->make<TH1I>("NGoodVtx","NGoodVtx",50,0,50);
//------------------Histogramms-------------------------------------------
// next histogramms are covered just selected events
//------------------DiMuon Vertex --------------------------------------------
    DiMuon_Vrt    = fs->make<TH1F>("DiMuonVrt","",30,0,30);
//------------------ Mass------------------------------------------------------
    DiMuon_Mass_fromVrt = fs->make<TH1F>("DiMass_frVrt","",45,Array);
    DYMass_hist   = fs->make<TH1F>("DYMass","DYMass", 45, Array);
    DYZpkM_hist   = fs->make<TH1F>("DYMass_Zpk","DYMass_Zpk", 60, 60, 120);
//--------------------Pt-------------------------------------------------------------
    DYPt_hist     = fs->make<TH1F>("DYPt","DYPt", 100, 0., 250.);
    DYPt_lead     = fs->make<TH1F>("DYLeadPt","DYLeadPt", 100, 0., 250.);
    DYPt_sub      = fs->make<TH1F>("DYSubPt","DYSubPt", 100, 0., 250.);
//--------------------Eta & Phi--------------------------------------------------------------
    DYEta_hist    = fs->make<TH1F>("DYEta","DYEta", 100, -3, 3.);
    DYEta_lead    = fs->make<TH1F>("DYLeadEta","DYLeadEta", 100, -3, 3.);
    DYEta_sub     = fs->make<TH1F>("DYSubEta","DYSubEta", 100, -3, 3.);
    DYPhi_hist    = fs->make<TH1F>("DYPhi","DYPhi", 100, -4., 4.);
    DYPhi_lead    = fs->make<TH1F>("DYLeadPhi","DYLeadPhi", 100, -4., 4.);
    DYPhi_sub     = fs->make<TH1F>("DYSubPhi","DYSubPhi", 100, -4., 4.);
//--------------------GenMass--------------------------------------------------------------
    DYGenMass       = fs->make<TH1F>("DYGenMass","DYGenMass", 45, Array);
    DYGenMass_nom   = fs->make<TH1F>("DYGenMass_nom","DYGenMass_nom", 45, Array);
    DYGenMass_den   = fs->make<TH1F>("DYGenMass_den","DYGenMass_den", 45, Array);
//--------------------Tree Branches--------------------------------------------------------
//--------------------Gen Part--------------------------------------------------------
   mytree->Branch("GenMass"          ,&GenMass            ,"GenMass/F");
   mytree->Branch("GenMass_den"      ,&GenMass_den        ,"GenMass_den/F");
   mytree->Branch("GenMass_nom"      ,&GenMass_nom        ,"GenMass_nom/F");
//------------------ Weight------------------------------------------------------
   mytree->Branch("mcWeight",&mcWeight_,"mcWeight/F");
//-------------------Pile up-----------------------------------------------
   mytree->Branch("puTrue"      ,&puTrue_);
//------------------ Mass------------------------------------------------------
   mytree->Branch("mmM_Sel"          ,&mmM_Sel            ,"mmM_Sel/F");
   mytree->Branch("mmM_Vrt_Good"     ,&mmM_Vrt_Good       ,"mmM_Vrt_Good/F");
   mytree->Branch("mmM_Vrt_Sel"      ,&mmM_Vrt_Sel        ,"mmM_Vrt_Sel/F");
//--------------------Eta & Phi--------------------------------------------------------------
   mytree->Branch("muPt"             ,&muPt);
   mytree->Branch("muPt_Good"        ,&muPt_Good);
   mytree->Branch("muPt_Sel"         ,&muPt_Sel);
   mytree->Branch("muEta"            ,&muEta);
   mytree->Branch("muEta_Good"       ,&muEta_Good);
   mytree->Branch("muEta_Sel"        ,&muEta_Sel);
   mytree->Branch("muPhi"            ,&muPhi);
   mytree->Branch("muPhi_Good"       ,&muPhi_Good);
   mytree->Branch("muPhi_Sel"        ,&muPhi_Sel);
//--------------------Charge-----------------------------------------------
   mytree->Branch("muChrg"           ,&muChrg);
   mytree->Branch("muChrg_Good"      ,&muChrg_Good);
   mytree->Branch("muChrg_Sel"       ,&muChrg_Sel);

//--------------------Addition--------------------------------------------------------------
   mytree->Branch("Angle_Good"   ,&Angle_Good,"Angle_Good/D");
   mytree->Branch("Angle_Sel"   ,&Angle_Sel,"Angle_Sel/D");
   mytree->Branch("vtxNormChi_Good" ,&vtxNormChi_Good,"vtxNormChi_Good/F");
   mytree->Branch("minVrt_Sel" ,&minVrt_Sel,"minVrt_Sel/F");
//--------------------Trigger---------------------------------------------------
   mytree->Branch("HLT_match"  ,&HLT_match);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Moge_highPtID::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

void 
Moge_highPtID::beginRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when ending the processing of a run  ------------

void 
Moge_highPtID::endRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when starting to processes a luminosity block  ------------

void 
Moge_highPtID::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method called when ending the processing of a luminosity block  ------------

void 
Moge_highPtID::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Moge_highPtID::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void Moge_highPtID::clearTree() {
  /*|-------------------------------------------------------|
    |                Run Info Part                           |
    |-------------------------------------------------------|*/     
	Run = -999;
	Event = -999;
	lumi = -999;
	bunch = -999;
	isRealData = -999;
      //vectors
  /*|-------------------------------------------------------|
    |                Trigger Part                           |
    |-------------------------------------------------------|*/
	HLT_nb.clear();
	HLT_name.clear();
	HLT_isaccept.clear();
	HLTObj_nbObj.clear();
	HLTObj_nbObjBefore.clear();
	HLTObj_pt.clear();
	HLTObj_eta.clear();
	HLTObj_phi.clear();
	HLTObj_collection.clear();
        HLT_match.clear();
  /*|-------------------------------------------------------|
    |                Gen Part                           |
    |-------------------------------------------------------|*/
	mcWeight_   = -999;
        GenMass     = -999;
        GenMass_nom = -999;
        GenMass_den = -999;
  /*|-------------------------------------------------------|
    |                Reco Part                           |
    |-------------------------------------------------------|*/
	mmM_Sel         = -999;
	mmM_Vrt_Good    = -999; 
	mmM_Vrt_Sel     = -999; 
	Angle_Good      = -999;
	Angle_Sel       = -999;
	vtxNormChi_Good = 999;
	minVrt_Sel      = 21;
	muPt.clear();
	muPt_Good.clear();
	muPt_Sel.clear();
	muEta.clear();
	muEta_Good.clear();
	muEta_Sel.clear();
	muPhi.clear();
	muPhi_Good.clear();
	muPhi_Sel.clear();
	muChrg.clear();
	muChrg_Good.clear();
	muChrg_Sel.clear();

  /*|-------------------------------------------------------|
    |                PileUp part                            |
    |-------------------------------------------------------|*/
	nPUInfo_ = -999;
	nPU_.clear();
	puBX_.clear();
	puTrue_.clear();
  /*|-------------------------------------------------------|
    |                PrimaryVertex Part                     |
    |-------------------------------------------------------|*/
	nbPv.clear();
	Nbdof.clear();
}
//define this as a plug-in
DEFINE_FWK_MODULE(Moge_highPtID);

bool Moge_highPtID::isPassHLT(){ 
  int nbMatch = 0;
  //std::cout<<" HLT_nb.size() = "<< HLT_nb.size() << endl;
  for(unsigned i=0; i < HLT_nb.size(); i++){
    if( (HLT_name[i] == "HLT_IsoMu20_v1" || 
         HLT_name[i] == "HLT_IsoTkMu20_v1") && HLT_isaccept[i] == 1 ) {
	  //std::cout<<"triggerName = "<< HLT_name[i] << endl;
      //cout<<"triggerName = "<<name<<endl;
      nbMatch++;
    }
  }
//std::cout<<" Matchpass = "<< nbMatch << endl;
  if(nbMatch>0) {
    return true;
  }
  else return false;
}
bool Moge_highPtID::RecoHLTMuonMatching(float RecoEta,float RecoPhi){
  int nbMatch  = 0;
  float deltaR = -10000.0;
   //std::cout<<" HLTObj_nbObj = "<< HLTObj_nbObj.size() << endl;
	
  for(unsigned i=0; i<HLTObj_nbObj.size(); i++){
    if( HLTObj_collection[i] == "HLT_IsoMu20_v1" || 
		HLTObj_collection[i] == "HLT_IsoTkMu20_v1"){

		deltaR = delR(HLTObj_eta[i],HLTObj_phi[i],RecoEta,RecoPhi);
		if(fabs(deltaR) > 0.1) continue;
                nbMatch++;
    }
  }
//std::cout<<" Matchobj = "<< nbMatch << endl;
  if(nbMatch>0){
  return true;
  }
  else return false;
}
float Moge_highPtID::delR(float eta1,float phi1,float eta2,float phi2){
	float mpi=3.14;
	float dp=std::abs(phi1-phi2);
	if (dp>mpi) dp-=float(2*mpi);
	return sqrt((eta1-eta2)*(eta1-eta2) + dp*dp);
	}
