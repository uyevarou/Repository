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

#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "Math/VectorUtil.h"

#include "TH1.h"
//#include "MSETools.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;
//
// class declaration
//

class TauEl_Both_b : public edm::EDAnalyzer {
   public:
      explicit TauEl_Both_b(const edm::ParameterSet&);
      ~TauEl_Both_b();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      bool isPassHLT();
      bool RecoHLTMuonMatching(float RecoEta,float RecoPhi);
      float delR(float eta1,float phi1,float eta2,float phi2);
      float TMass_F(float pt3lep, float px3lep, float py3lep, float met, float metPhi);

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
      edm::EDGetToken electronsToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<reco::GsfElectronCoreCollection> gsfelectronToken_;
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

      // ID decisions objects
      edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;

      // One example of full information about the cut flow
      edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleMediumIdFullInfoMapToken_;

    // Verbose output for ID
    bool verboseIdFlag_;

      // structures
      struct PATRONUS{
   TLorentzVector p4;
       float Spt;
       float Seta;
       float Sphi;
       int Scharge;
       int Jcharge;
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
       float IsoTrk;
       float IsoPf;
       float Sdz;
       float SEn;
       bool isMedium;
       bool STauIdIso;
       bool STauIdMuRj;
       bool STauIdMVA;
       float Smass;
       };
	static bool lepSortingRule(PATRONUS x, PATRONUS y) {return x.p4.Pt() > y.p4.Pt();}
	static bool lepSortingRule_D(PATRONUS x, PATRONUS y) {return x.Spt > y.Spt;}

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
         int status;
       float rap;
       float en;
       float pt;
       float eta;
       float phi;
       float px;
       float py;
       float pz;
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

      TH1F *DiEle_Vrt;

      TH1F *DiEle_Mass_fromVrt;
      TH1F *El_DYMass_hist;
      TH1F *El_DYZpkM_hist;

      TH1F *El_DYPt_hist;
      TH1F *El_DYPt_lead; 
      TH1F *El_DYPt_sub;

      TH1F *El_DYEta_hist;
      TH1F *El_DYEta_lead;
      TH1F *El_DYEta_sub;
      TH1F *El_DYPhi_hist;
      TH1F *El_DYPhi_lead;
      TH1F *El_DYPhi_sub;

  /*|-------------------------------------------------------|
    |                Run Info Part                          |
    |-------------------------------------------------------|*/     
      int Run;
      int Event;
      int lumi;
      int bunch;
      int isRealData;
      bool isMC_;
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
      std::vector<int>    genPdgId;
      std::vector<int>    genStatus;
      std::vector<int>    genMotherId;
      std::vector<float>  genPt;
      std::vector<float>  genEta;
      std::vector<float>  genPhi;
      std::vector<int>    genCharge;
      std::vector<float>  genPx;
      std::vector<float>  genPy;
      std::vector<float>  genPz;
      std::vector<float>  genP;

      std::vector<bool>   fromHProcessFinalState_;
      std::vector<bool>   fromHProcessDecayed_;
      std::vector<bool>   isHProcess_;
      std::vector<float>  genPreRap;
      std::vector<float>  genPreEn;
      std::vector<float>  genPtP4;
      std::vector<float>  genEtaP4;
      std::vector<float>  genPhiP4;
      std::vector<float>  genPxP4;
      std::vector<float>  genPyP4;
      std::vector<float>  genPzP4;
      std::vector<float>  genPP4;

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
      int gdInx_mu1;
      int gdInx_mu2;
      std::vector<float>  muPt_Good;
      std::vector<float>  muPt_Sel;
      std::vector<float>  muEta_Good;
      std::vector<float>  muEta_Sel;
      std::vector<float>  muPhi_Good;
      std::vector<float>  muPhi_Sel;
      std::vector<int>    muChrg_Good;
      std::vector<int>    muChrg_Sel;
   //----------------Initial tree----------------------------
     std::vector<bool> Mu_isTightMuon;
     std::vector<bool> Mu_isMediumMuon;
     std::vector<bool> Mu_isGlobalMuon;
     std::vector<bool> Mu_isHighPtMuon;
     std::vector<bool> Mu_isTrackerMuon;
     std::vector<float>  muPtP4;
     std::vector<float>  muEtaP4;
     std::vector<float>  muPhiP4;
     std::vector<float>  muPxP4;
     std::vector<float>  muPyP4;
     std::vector<float>  muPzP4;
     std::vector<float>  muPP4;
     std::vector<int>    muChargP4;
     std::vector<float>  muPtTuneP;
     std::vector<float>  muEtaTuneP;
     std::vector<float>  muPhiTuneP;
     std::vector<float>  muPxTuneP;
     std::vector<float>  muPyTuneP;
     std::vector<float>  muPzTuneP;
     std::vector<float>  muPTuneP;
     std::vector<int>    muChargTuneP;
     std::vector<int>    muglobalMuonHits;
     std::vector<int>    munbMatchSt;
     std::vector<int>    mupixelHits;
     std::vector<int>    mutrackerLayers;
     std::vector<float>  mudxy;
     std::vector<float>  mupterror;
     std::vector<float>  muisoTrk;
     std::vector<float>  muisoPf;
     std::vector<float>  mudz;
//dimuon vectors
     std::vector<float>  init_vtxNormChi;
     std::vector<float>  init_vrtMass;

  /*|-------------------------------------------------------|
    |                PileUp part                            |
    |-------------------------------------------------------|*/
	Int_t                 nPUInfo_;
	std::vector<int>      nPU_;
	std::vector<int>      puBX_;
	std::vector<float>    puTrue_;
        int nPV;
  /*|-------------------------------------------------------|
    |                PrimaryVertex Part                     |
    |-------------------------------------------------------|*/
	std::vector<int>    nbPv;
  	std::vector<int>    Nbdof;
        std::vector<float>  nbPv_del;
        Int_t nbPv_dl; // number of reconsrtucted primary vertices 
        //nbPv_dl;
  /*|-------------------------------------------------------|
    |                Electron Part                     |
    |-------------------------------------------------------|*/
	std::vector<float>    elPt_Good;
	std::vector<float>    elEta_Good;
	std::vector<float>    elPhi_Good;
	std::vector<int>      elChrg_Good;
	std::vector<float>    elP_Good;
	std::vector<float>    elEnSC_Good;
	std::vector<float>    elEtaSC_Good;
	std::vector<float>    elPhiSC_Good;

        float eeM_Sel;
	float El_mmM_Vrt_Good; 
	float El_mmM_Vrt_Sel; 
	float El_Angle_Good;
	float El_Angle_Sel;
	float El_vtxNormChi_Good;
	float El_minVrt_Sel;
	int gdInx_el1;
	int gdInx_el2;

	std::vector<float>    elPt_Sel;
	std::vector<float>    elEta_Sel;
	std::vector<float>    elPhi_Sel;
	std::vector<int>      elChrg_Sel;
   //--------Init Tree--------------------------------------

	std::vector<bool>    elisMedium;

	std::vector<float>    elPtP4;
	std::vector<float>    elEtaP4;
	std::vector<float>    elPhiP4;
	std::vector<float>    elPxP4;
	std::vector<float>    elPyP4;
	std::vector<float>    elPzP4;
	std::vector<float>    elPP4;
	std::vector<int>      elChargP4;
        
	std::vector<float>    elEnSC;
	std::vector<float>    elEtaSC;
	std::vector<float>    elPhiSC;

	std::vector<float>    El_init_vtxNormChi;
	std::vector<float>    El_init_vrtMass;
 /*|-------------------------------------------------------|
   |                Tau Part                               |
   |-------------------------------------------------------|*/
	std::vector<float>    tauPt;
	std::vector<float>    tauEta;
	std::vector<float>    tauPhi;
	std::vector<int>      tauCharge;
	std::vector<float>    tauMass;
	std::vector<float>    taudxy;
	std::vector<bool>     TauIdIso;
	std::vector<bool>     TauIdMuRj;
	std::vector<bool>     TauIdMVA;
	std::vector<float>    MassDYtaumuSS;
	std::vector<float>    MassDYtaumuOS;
 /*|-------------------------------------------------------|
   |                MET Part                               |
   |-------------------------------------------------------|*/
	float pfMET_;
	float pfMETPhi_;
	float MuMetTranverseMass;
	float pfMETSumEt_;
 /*|-------------------------------------------------------|
   |                Photon Part                            |
   |-------------------------------------------------------|*/
	std::vector<float>    ptPhoton_;
	std::vector<float>    etaPhoton_;
	std::vector<float>    phiPhoton_;
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
TauEl_Both_b::TauEl_Both_b(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
   //electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
   //gsfelectronToken_(consumes<reco::GsfElectronCoreCollection>(iConfig.getParameter<edm::InputTag>("gsfelectrons"))),
   tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
   photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
   jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
   fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
   metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
   triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
   triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
   triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),   
   genPackedCollToken_(mayConsume<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGenParticles"))),   
   genPrunedCollToken_(mayConsume<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
   genJetToken_(mayConsume<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
   eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
   eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
   eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
   eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
   eleHEEPIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap")))
{
   mSrcRho            = iConfig.getParameter<edm::InputTag>             ("srcRho");
   puCollection_      = consumes<vector<PileupSummaryInfo> >    (iConfig.getParameter<InputTag>("pileupCollection"));
   conversionsToken_  = mayConsume< reco::ConversionCollection >(iConfig.getParameter<edm::InputTag>("conversions"));
   beamSpotToken_     = consumes<reco::BeamSpot> (iConfig.getParameter <edm::InputTag>("beamSpot"));
   electronsToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
("electrons"));
   isMC_ = iConfig.getUntrackedParameter<bool>("isMC");

}


TauEl_Both_b::~TauEl_Both_b()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TauEl_Both_b::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  clearTree();

  Run   = iEvent.id().run();
  Event = iEvent.id().event();
  lumi  = iEvent.luminosityBlock();
  bunch = iEvent.bunchCrossing();
  isRealData = iEvent.isRealData() ? 1:0;
/*|-------------------------------------------------------|
  |                PileUp part                            |
  |-------------------------------------------------------|*/
  if (!isRealData){
  edm::Handle<vector<PileupSummaryInfo> > genPileupHandle;
  iEvent.getByToken(puCollection_, genPileupHandle);
  if (genPileupHandle.isValid()) {
    for (vector<PileupSummaryInfo>::const_iterator pu = genPileupHandle->begin(); pu != genPileupHandle->end(); ++pu) {
      if (pu->getBunchCrossing() == 0) {
        hPU_->Fill(pu->getPU_NumInteractions());
        hPUTrue_->Fill(pu->getTrueNumInteractions());
      
        nPU_   .push_back(pu->getPU_NumInteractions());
        puTrue_.push_back(pu->getTrueNumInteractions());
        //hPUTrueTrue_->Fill(pu->getTrueNumInteractions());
        puBX_  .push_back(pu->getBunchCrossing());

        nPUInfo_++;
      }
    }
  } 
 } 
/*|-------------------------------------------------------|
  |                PrimaryVertex Part                     |
  |-------------------------------------------------------|*/
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  if (vertices->empty()) return;
  nbPv_dl = vertices->size();
  int nGoodVtxsId = 0;

  VertexCollection::const_iterator nGoodVtxs = vertices->end();
  for (VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++nGoodVtxsId) {
    nGoodVtxs = vtx;
    break;
  }

  if (nGoodVtxs == vertices->end())
  return; // skip event if there are no good PVs


  h_NbVtx -> Fill(nbPv_dl);
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
      //if( src.pt() < 10.0 ) continue;
      //if( fabs(src.eta()) > 2.5 ) continue;
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
  if (isMC_){
  vector<GENPARTICLE> myGenLeptons;
  edm::Handle<pat::PackedGenParticleCollection> genPackedColl;
  edm::Handle<reco::GenParticleCollection> genPrunedColl;
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
  

   for (const reco::GenParticle &i_gen : *genPrunedColl) {           
      // ----------- Muons and Electrons -------------
      // set lepton 4-vector and ID
      fromHProcessFinalState_.push_back(i_gen.fromHardProcessFinalState());
      fromHProcessDecayed_.push_back(i_gen.fromHardProcessDecayed());
      isHProcess_.push_back(i_gen.isHardProcess());

      GENPARTICLE aGenLepton;
      TLorentzVector lepP4GEN(i_gen.p4().Px(),i_gen.p4().Py(),i_gen.p4().Pz(),i_gen.p4().E());
      aGenLepton.p4 = lepP4GEN;
      aGenLepton.pdgId = i_gen.pdgId();
      aGenLepton.charge = i_gen.charge();
      aGenLepton.status= i_gen.status();
      aGenLepton.rap = i_gen.rapidity();
      aGenLepton.en  = i_gen.energy();
      aGenLepton.pt  = i_gen.pt();
      aGenLepton.eta = i_gen.eta();
      aGenLepton.phi = i_gen.phi();
      aGenLepton.px  = i_gen.px();
      aGenLepton.py  = i_gen.py();
      aGenLepton.pz  = i_gen.pz();

 if ( abs(i_gen.pdgId()) != 13 && abs(i_gen.pdgId()) != 11 ) continue; 
 if ( i_gen.status() != 1 ) continue; 


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
      
           // Keep leptons only within kinematic and geometric acceptance
           // if ( (aGenLepton.p4.Pt() < 22) && (fabs(aGenLepton.p4.Eta()) > 2.4) ) continue; 
     myGenLeptons.push_back(aGenLepton);
     } // end PackedGenParticles loop

  for(unsigned l = 0; l < myGenLeptons.size(); l++) {
     genPdgId.push_back(myGenLeptons[l].pdgId);
     genStatus.push_back(myGenLeptons[l].status);
     genMotherId.push_back(myGenLeptons[l].motherId);
     genPtP4.push_back(myGenLeptons[l].p4.Pt());
     genEtaP4.push_back(myGenLeptons[l].p4.Eta());
     genPhiP4.push_back(myGenLeptons[l].p4.Phi());
     genCharge.push_back(myGenLeptons[l].charge);
     genPxP4.push_back(myGenLeptons[l].p4.Px());
     genPyP4.push_back(myGenLeptons[l].p4.Py());
     genPzP4.push_back(myGenLeptons[l].p4.Pz());
     genPP4.push_back(myGenLeptons[l].p4.P());
     genPreRap.push_back(myGenLeptons[l].rap);
     genPreEn.push_back(myGenLeptons[l].en);
     genPx.push_back(myGenLeptons[l].px);
     genPy.push_back(myGenLeptons[l].py);
     genPz.push_back(myGenLeptons[l].pz);
     genPt.push_back(myGenLeptons[l].pt);
     genEta.push_back(myGenLeptons[l].eta);
     genPhi.push_back(myGenLeptons[l].phi);
  }

    sort(myGenLeptons.begin(), myGenLeptons.end(), lepSortingRuleGen);
    TLorentzVector llP4GEN(0,0,0,0);
    if (myGenLeptons.size() > 1) llP4GEN = myGenLeptons[0].p4 + myGenLeptons[1].p4;
    GenMass = llP4GEN.M();
    DYGenMass -> Fill(GenMass); 
}

/*
    if (myGenLeptons.size() > 1
        && myGenLeptons[0].pdgId * myGenLeptons[1].pdgId == -169
        && myGenLeptons[0].charge*myGenLeptons[1].charge<0)
        //&& ((myGenLeptons[0].motherId == 23 && myGenLeptons[1].motherId == 23) || (myGenLeptons[0].motherId == 22 && myGenLeptons[1].motherId == 22)))
        {
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
*/
  /*|-------------------------------------------------------|
    |                Electron Part                          |
    |-------------------------------------------------------|*/

    edm::Handle<edm::View<reco::GsfElectron> > electrons;
    //edm::Handle<reco::GsfElectronCoreCollection> electrons;
    iEvent.getByToken(electronsToken_, electrons);

    edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
    iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);

    vector<PATRONUS> GoodRecoEl;
    vector<PATRONUS> GoodLeptons;
    vector<PATRONUS> InitEl;
    //vector<int> passMediumId_;

    for (size_t i = 0; i < electrons->size(); ++i){
    const auto el = electrons->ptrAt(i);
    PATRONUS elON;
    TLorentzVector elONP4(el->p4().Px(), el->p4().Py(), el->p4().Pz(), el->p4().E() );
    elON.p4      = elONP4;
    elON.Jcharge = el->charge();
    elON.Sp      = el->p();

    //double R = sqrt(el->superCluster()->x()*el->superCluster()->x() + el->superCluster()->y()*el->superCluster()->y() +el->superCluster()->z()*el->superCluster()->z());
    //double Rt = sqrt(el->superCluster()->x()*el->superCluster()->x() + el->superCluster()->y()*el->superCluster()->y());

    elON.SEn = el->superCluster()->energy();
    //preEnSC_.push_back(el->superCluster()->preshowerEnergy());
    //rawEnSC_.push_back(el->superCluster()->rawEnergy());
    //elON.SEn =  (el->superCluster()->energy())*(Rt/R);
    elON.Seta = el->superCluster()->eta();
    elON.Sphi = el->superCluster()->phi();

   
    // Isolation
    //GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();

    edm::ESHandle<TransientTrackBuilder> ttkb;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb);
    reco::TransientTrack ttv;

    reco::GsfTrackRef theTrack = el->gsfTrack();

    ttv = ttkb->build(theTrack);
    elON.Strack  = ttv;
    bool isPassMedium = (*medium_id_decisions)[el];
    elON.isMedium = isPassMedium;
    
    InitEl.push_back(elON);

    if (el->pt() < 10) continue;
    if (fabs(el->eta()) >= 2.5) continue;
    if (isPassMedium == 0) continue;
    bool fireElHLT = isPassHLT();
    if(fireElHLT == 0) continue;

    GoodRecoEl.push_back(elON);
    GoodLeptons.push_back(elON);
    }

  sort(GoodRecoEl.begin(), GoodRecoEl.end(), lepSortingRule_D);

  
  for (unsigned j = 0; j < GoodRecoEl.size(); j++) {
      
      elPt_Good.push_back(GoodRecoEl[j].p4.Pt());
      elEta_Good.push_back(GoodRecoEl[j].p4.Eta());
      elPhi_Good.push_back(GoodRecoEl[j].p4.Phi());
      elChrg_Good.push_back(GoodRecoEl[j].Scharge);
      elP_Good.push_back(GoodRecoEl[j].p4.P());
      elEnSC_Good.push_back(GoodRecoEl[j].SEn);
      elEtaSC_Good.push_back(GoodRecoEl[j].Seta);
      elPhiSC_Good.push_back(GoodRecoEl[j].Sphi);

  }

  if (GoodRecoEl.size() > 1){

   for(unsigned led = 0; led < GoodRecoEl.size(); led++) {
      if (!(GoodRecoEl[led].p4.Pt() > 22)) continue;

      for(unsigned sub = 0; sub < GoodRecoEl.size(); sub++) {
         if (led == sub) continue;
         // Angle
         TVector3 El1(GoodRecoEl[led].p4.Px(),GoodRecoEl[led].p4.Py(), GoodRecoEl[led].p4.Pz());
         TVector3 El2(GoodRecoEl[sub].p4.Px(), GoodRecoEl[sub].p4.Py(), GoodRecoEl[sub].p4.Pz());
         double cos_angle = El1.Dot(El2) / GoodRecoEl[led].p4.P() / GoodRecoEl[sub].p4.P();
         El_Angle_Good = cos_angle;
         //cout<<"angle ="<<cos_angle<<endl;
         //trigger
         bool HLT_el1 = RecoHLTMuonMatching(GoodRecoEl[led].p4.Eta(),GoodRecoEl[led].p4.Phi());
         bool HLT_el2 = RecoHLTMuonMatching(GoodRecoEl[sub].p4.Eta(),GoodRecoEl[sub].p4.Phi());
         //vertex
         std::vector<reco::TransientTrack> ttV;
         ttV.clear();
         ttV.push_back(GoodRecoEl[led].Strack);
         ttV.push_back(GoodRecoEl[sub].Strack);
         KalmanVertexFitter kvf(true);
         CachingVertex<5> vtx = kvf.vertex(ttV);
         float vtxNormChi = vtx.totalChiSquared()/vtx.degreesOfFreedom();
         //std::cout<<"vtxNormChi =" << vtxNormChi << std::endl;
         InvariantMassFromVertex imfv;
         static const double ele_mass = 0.0005110;
         El_vtxNormChi_Good = vtxNormChi;
         Measurement1D massOK = imfv.invariantMass(vtx, ele_mass);
         El_mmM_Vrt_Good = massOK.value();

         if (cos_angle < -0.999987) continue;
         if (HLT_el1 == 0 || HLT_el2 == 0) continue;
         if (!(GoodRecoEl[sub].p4.Pt() > 10)) continue;
         if (!(GoodRecoEl[led].Jcharge * GoodRecoEl[sub].Jcharge < 0)) continue;
         if (!(vtxNormChi < 20.0)) continue;
         if (vtxNormChi < El_minVrt_Sel) {
             El_minVrt_Sel = vtxNormChi;
             gdInx_el1 = led;
             if (gdInx_el1 > 0) cout<<"ind ="<<gdInx_el1<<endl;
             gdInx_el2 = sub;
             Measurement1D massO = imfv.invariantMass(vtx, ele_mass);
             El_mmM_Vrt_Sel = massO.value();
             El_Angle_Sel = cos_angle;
          }  
       }
     }
  }
   DiEle_Vrt -> Fill(El_minVrt_Sel);
   DiEle_Mass_fromVrt -> Fill(El_mmM_Vrt_Sel);

   if (gdInx_el1 > -1 && gdInx_el2 > -1){
   elPt_Sel.push_back(GoodRecoEl[gdInx_el1].p4.Pt());
   elPt_Sel.push_back(GoodRecoEl[gdInx_el2].p4.Pt());
   El_DYPt_hist     -> Fill(GoodRecoEl[gdInx_el1].p4.Pt());
   El_DYPt_hist     -> Fill(GoodRecoEl[gdInx_el2].p4.Pt());
   El_DYPt_lead     -> Fill(GoodRecoEl[gdInx_el1].p4.Pt());
   El_DYPt_sub      -> Fill(GoodRecoEl[gdInx_el2].p4.Pt());

   elEta_Sel.push_back(GoodRecoEl[gdInx_el1].p4.Eta());
   elEta_Sel.push_back(GoodRecoEl[gdInx_el2].p4.Eta());
   El_DYEta_hist    -> Fill(GoodRecoEl[gdInx_el1].p4.Eta());
   El_DYEta_hist    -> Fill(GoodRecoEl[gdInx_el2].p4.Eta());
   El_DYEta_lead    -> Fill(GoodRecoEl[gdInx_el1].p4.Eta());
   El_DYEta_sub     -> Fill(GoodRecoEl[gdInx_el2].p4.Eta());

   elPhi_Sel.push_back(GoodRecoEl[gdInx_el1].p4.Phi());
   elPhi_Sel.push_back(GoodRecoEl[gdInx_el2].p4.Phi());
   El_DYPhi_hist    -> Fill(GoodRecoEl[gdInx_el1].p4.Phi());
   El_DYPhi_hist    -> Fill(GoodRecoEl[gdInx_el2].p4.Phi());
   El_DYPhi_lead    -> Fill(GoodRecoEl[gdInx_el1].p4.Phi());
   El_DYPhi_sub     -> Fill(GoodRecoEl[gdInx_el2].p4.Phi());

   elChrg_Sel.push_back(GoodRecoEl[gdInx_el1].Scharge);
   elChrg_Sel.push_back(GoodRecoEl[gdInx_el2].Scharge);

   TLorentzVector el1, el2;
   el1.SetPtEtaPhiM(GoodRecoEl[gdInx_el1].p4.Pt(),GoodRecoEl[gdInx_el1].p4.Eta(),GoodRecoEl[gdInx_el1].p4.Phi(), 0.0005110);
   el2.SetPtEtaPhiM(GoodRecoEl[gdInx_el2].p4.Pt(),GoodRecoEl[gdInx_el2].p4.Eta(),GoodRecoEl[gdInx_el2].p4.Phi(), 0.0005110);
   eeM_Sel = (el1+el2).M();
   El_DYMass_hist -> Fill(eeM_Sel);
   if (eeM_Sel > 60 && eeM_Sel < 120) El_DYZpkM_hist -> Fill(eeM_Sel);

  }
  
// filling initial tree
  for(unsigned l = 0; l < InitEl.size(); l++) {
     elPtP4.push_back(InitEl[l].p4.Pt());
     elEtaP4.push_back(InitEl[l].p4.Eta());
     elPhiP4.push_back(InitEl[l].p4.Phi());
     elPxP4.push_back(InitEl[l].p4.Px());
     elPyP4.push_back(InitEl[l].p4.Py());
     elPzP4.push_back(InitEl[l].p4.Pz());
     elPP4.push_back(InitEl[l].p4.P());
     elChargP4.push_back(InitEl[l].Jcharge);

     elEnSC.push_back(InitEl[l].SEn);
     elEtaSC.push_back(InitEl[l].Seta);
     elPhiSC.push_back(InitEl[l].Sphi);
     elisMedium.push_back(InitEl[l].isMedium);

  }

  int nbInEl = InitEl.size();
  sort(InitEl.begin(), InitEl.end(), lepSortingRule_D);

   for(int ld = 0; ld < nbInEl; ld++) {
      for(int sb = 0; sb < nbInEl; sb++) {
         if (ld == sb) continue;
         std::vector<reco::TransientTrack> ttV_D;
         ttV_D.clear();
         ttV_D.push_back(InitEl[ld].Strack);
         ttV_D.push_back(InitEl[sb].Strack);
         KalmanVertexFitter kvf(true);
         CachingVertex<5> vtx = kvf.vertex(ttV_D);
         float vtxNormChi = vtx.totalChiSquared()/vtx.degreesOfFreedom();
         InvariantMassFromVertex imfv;
         static const double ele_mass = 0.0005110;
         El_init_vtxNormChi.push_back(vtxNormChi);
         Measurement1D massK = imfv.invariantMass(vtx, ele_mass);
         El_init_vrtMass.push_back(massK.value());
         //cout<<ld<<sb<<"="<<vtxNormChi<<endl;
       }
     }




  /*|-------------------------------------------------------|
    |                Muon Part                              |
    |-------------------------------------------------------|*/

  vector<PATRONUS> GoodRecoMu;
  vector<PATRONUS> InitMuons;
  const reco::Vertex &PV_na = vertices->front();

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  for (const pat::Muon &mu : *muons) {

    Mu_isMediumMuon.push_back(mu.isMediumMuon());
    Mu_isTightMuon.push_back(mu.isTightMuon(vertices->at(0)));
    Mu_isHighPtMuon.push_back(mu.isHighPtMuon(vertices->at(0)));
    Mu_isGlobalMuon.push_back(mu.isGlobalMuon());
    Mu_isTrackerMuon.push_back(mu.isTrackerMuon());

    PATRONUS muON;
    TLorentzVector muONP4(mu.p4().Px(), mu.p4().Py(), mu.p4().Pz(), mu.p4().E() );
    muON.p4      = muONP4;
    muON.Jcharge = mu.charge();

    const reco::TrackRef& tunePTrack  = mu.tunePMuonBestTrack();
    reco::TrackRef glbTrack      = mu.globalTrack();
    reco::TrackRef trackerTrack  = mu.innerTrack();
    reco::TrackRef muonTrack     = mu.outerTrack();
    reco::TrackRef MuonBestTrack = mu.muonBestTrack();

    muON.matchedMuStations = mu.numberOfMatchedStations();
    muON.Sdz = fabs(mu.muonBestTrack()->dz(nGoodVtxs -> position()));

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

    muON.dxy = fabs(mu.muonBestTrack()->dxy(nGoodVtxs -> position()));
    muON.ptError  = mu.muonBestTrack()->ptError()/mu.muonBestTrack()->pt();
    muON.IsoTrk = mu.isolationR03().sumPt/mu.pt();

    edm::ESHandle<TransientTrackBuilder> ttkb;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb);
    reco::TransientTrack ttv;
    ttv = ttkb->build(tunePTrack);
    muON.Strack  = ttv;
    float pfiso = (mu.pfIsolationR04().sumChargedHadronPt + max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5*mu.pfIsolationR04().sumPUPt))/mu.pt();
    muON.IsoPf  = pfiso;  
    muON.isMedium = mu.isMediumMuon();      

    InitMuons.push_back(muON);

    if( !(muON.Spt > 10 && abs(muON.Seta) < 2.4)) continue;
    if( !mu.isGlobalMuon()) continue;
    if( !(muON.globalMuonHits > 0)) continue;
    if( !(muON.matchedMuStations > 1)) continue;
    if( !(muON.ptError < 0.3)) continue;
    if( !(muON.dxy < 0.2)) continue;
    if( !(muON.pixelHits > 0)) continue;
    if( !(muON.trackerLayers > 5)) continue;  
    if (muON.IsoTrk >= 0.1) continue;
    bool fireHLT = isPassHLT();
    if(fireHLT == 0) continue;

    GoodRecoMu.push_back(muON);
    GoodLeptons.push_back(muON);
  }
  sort(GoodRecoMu.begin(), GoodRecoMu.end(), lepSortingRule_D);
  
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
   DiMuon_Vrt -> Fill(minVrt_Sel);
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
   if (mmM_Sel > 60 && mmM_Sel < 120) DYZpkM_hist -> Fill(mmM_Sel);

  }
  
// filling initial tree
  for(unsigned l = 0; l < InitMuons.size(); l++) {
     muPtP4.push_back(InitMuons[l].p4.Pt());
     muEtaP4.push_back(InitMuons[l].p4.Eta());
     muPhiP4.push_back(InitMuons[l].p4.Phi());
     muPxP4.push_back(InitMuons[l].p4.Px());
     muPyP4.push_back(InitMuons[l].p4.Py());
     muPzP4.push_back(InitMuons[l].p4.Pz());
     muPP4.push_back(InitMuons[l].p4.P());
     muChargP4.push_back(InitMuons[l].Jcharge);

     muPtTuneP.push_back(InitMuons[l].Spt);
     muEtaTuneP.push_back(InitMuons[l].Seta);
     muPhiTuneP.push_back(InitMuons[l].Sphi);
     muPxTuneP.push_back(InitMuons[l].Spx);
     muPyTuneP.push_back(InitMuons[l].Spy);
     muPzTuneP.push_back(InitMuons[l].Spz);
     muPTuneP.push_back(InitMuons[l].Sp);
     muChargTuneP.push_back(InitMuons[l].Scharge);

     muglobalMuonHits.push_back(InitMuons[l].globalMuonHits);
     munbMatchSt.push_back(InitMuons[l].matchedMuStations);
     mupixelHits.push_back(InitMuons[l].pixelHits);    
     mutrackerLayers.push_back(InitMuons[l].trackerLayers); 
     mudxy.push_back(InitMuons[l].dxy);
     mupterror.push_back(InitMuons[l].ptError);
     muisoTrk.push_back(InitMuons[l].IsoTrk);
     muisoPf.push_back(InitMuons[l].IsoPf);
     mudz.push_back(InitMuons[l].Sdz);


  }

  int nbInMu = InitMuons.size();
  sort(InitMuons.begin(), InitMuons.end(), lepSortingRule_D);
   for(int ld = 0; ld < nbInMu; ld++) {
      for(int sb = 0; sb < nbInMu; sb++) {
         if (ld == sb) continue;
         std::vector<reco::TransientTrack> ttV_D;
         ttV_D.clear();
         ttV_D.push_back(InitMuons[ld].Strack);
         ttV_D.push_back(InitMuons[sb].Strack);
         KalmanVertexFitter kvf(true);
         CachingVertex<5> vtx = kvf.vertex(ttV_D);
         float vtxNormChi = vtx.totalChiSquared()/vtx.degreesOfFreedom();
         InvariantMassFromVertex imfv;
         static const double muon_mass = 0.1056583;
         init_vtxNormChi.push_back(vtxNormChi);
         Measurement1D massK = imfv.invariantMass(vtx, muon_mass);
         init_vrtMass.push_back(massK.value());
       }
   }
 /*|-------------------------------------------------------|
   |                Tau Part                               |
   |-------------------------------------------------------|*/
  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauToken_, taus);

  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  const pat::MET &met = mets->front();
    pfMET_       = met.et();
    pfMETPhi_    = met.phi();
    pfMETSumEt_    = met.sumEt();

  vector<PATRONUS> InitTau;
  vector<float> MassDYtaumu;
  vector<float> MassDYtau;

  for (const pat::Tau &tau : *taus) {
    PATRONUS tauON;
    tauON.Spt = tau.pt();
    tauON.Seta = tau.eta();
    tauON.Jcharge = tau.charge();
    tauON.Sphi = tau.phi();
    tauON.Smass = tau.mass();
    tauON.dxy = tau.dxy();
    tauON.Sdz = tau.dxy();
    bool VTauIdIso = tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
    bool VtauByTightMuonRejection3 = tau.tauID("againstMuonTight3");
    bool VtauByMVA5LooseElectronRejection =  tau.tauID("againstElectronLooseMVA5");
    tauON.STauIdIso  = VTauIdIso;
    tauON.STauIdMuRj = VtauByTightMuonRejection3;
    tauON.STauIdMVA  = VtauByMVA5LooseElectronRejection;
    InitTau.push_back(tauON);    
  }
  sort(InitTau.begin(), InitTau.end(), lepSortingRule_D);
  for (unsigned mj = 0; mj < GoodRecoMu.size(); mj++) {
      for (unsigned tj = 0; tj < InitTau.size(); tj++) {
      //muon
      if(GoodRecoMu[mj].Spt < 18.0) continue;
      if(fabs(GoodRecoMu[mj].Seta) > 2.1) continue;
      if(GoodRecoMu[mj].IsoPf >= 0.25) continue;
      if(!GoodRecoMu[mj].isMedium) continue;
      if(GoodRecoMu[mj].dxy >= 0.045) continue;  
      if(GoodRecoMu[mj].Sdz  >= 0.2) continue;  

      //tau
      if(InitTau[tj].Spt < 20.0) continue;
      if(fabs(InitTau[tj].Seta) > 2.3) continue;
      if(InitTau[tj].STauIdIso == 0) continue; 
      if(InitTau[tj].STauIdMuRj == 0) continue;
      if(InitTau[tj].STauIdMVA == 0) continue;
      if(fabs(InitTau[tj].dxy) >= 0.05) continue;
    
      //mu-tau
      TLorentzVector Mu4Momentum,Tau4Momentum,Z4Momentum;
      bool  OS = GoodRecoMu[mj].Scharge * InitTau[tj].Jcharge < 0;   
      float tmass = TMass_F(GoodRecoMu[mj].Spt, GoodRecoMu[mj].Spt*cos(GoodRecoMu[mj].Sphi),GoodRecoMu[mj].Spt*sin(GoodRecoMu[mj].Sphi),pfMET_,pfMETPhi_);
      MuMetTranverseMass = tmass;
      //cout<<"tmass"<<tmass<<endl;
      if(tmass > 40) continue;
      float deltaRtauMu = delR(GoodRecoMu[mj].Seta,GoodRecoMu[mj].Sphi,InitTau[tj].Seta,InitTau[tj].Sphi);
      //if(!(Mu4Momentum.DeltaR(Tau4Momentum) > 0.5)) continue;    
      if(!(deltaRtauMu > 0.5)) continue;    
      //cout<<"Taumass"<<tmass<<endl;
      if (OS == 1){
        Mu4Momentum.SetPtEtaPhiM(GoodRecoMu[mj].Spt,GoodRecoMu[mj].Seta,GoodRecoMu[mj].Sphi,0.1056583);
        Tau4Momentum.SetPtEtaPhiM(InitTau[tj].Spt,InitTau[tj].Seta,InitTau[tj].Sphi,InitTau[tj].Smass);
        Z4Momentum=Mu4Momentum+Tau4Momentum;
        MassDYtaumuOS.push_back(Z4Momentum.M());
      }
      if (OS == 0){
        Mu4Momentum.SetPtEtaPhiM(GoodRecoMu[mj].Spt,GoodRecoMu[mj].Seta,GoodRecoMu[mj].Sphi,0.1056583);
        Tau4Momentum.SetPtEtaPhiM(InitTau[tj].Spt,InitTau[tj].Seta,InitTau[tj].Sphi,InitTau[tj].Smass);
        Z4Momentum=Mu4Momentum+Tau4Momentum;
        MassDYtaumuSS.push_back(Z4Momentum.M());
      }
    }
  }

  for (unsigned l = 0; l < InitTau.size(); l++) {
     tauPt.push_back(InitTau[l].Spt);
     tauEta.push_back(InitTau[l].Seta);
     tauPhi.push_back(InitTau[l].Sphi);
     tauCharge.push_back(InitTau[l].Jcharge);
     tauMass.push_back(InitTau[l].Smass);
     taudxy.push_back(InitTau[l].dxy);
     TauIdIso.push_back(InitTau[l].STauIdIso);
     TauIdMuRj.push_back(InitTau[l].STauIdMuRj);
     TauIdMVA.push_back(InitTau[l].STauIdMVA);
  } 
 /*|-------------------------------------------------------|
   |                Photons Part                           |
   |-------------------------------------------------------|*/

  edm::Handle<pat::PhotonCollection> photons;
  iEvent.getByToken(photonToken_, photons);

    for (const pat::Photon &pho : *photons) {

	ptPhoton_  .push_back(pho.pt());
	etaPhoton_ .push_back(pho.superCluster()->eta());
	phiPhoton_ .push_back(pho.superCluster()->phi());
     }


  mytree->Fill();
}
// ------------ method called once each job just before starting event loop  ------------
void 
TauEl_Both_b::beginJob()
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
    //hPUTrueTrue_    = fs->make<TH1F>("hPUTrueTrue","number of true pilepu", 50, 0, 50);
//------------------ Vertex------------------------------------------------------
    h_NbVtx       = fs->make<TH1I>("NVtx","NVtx",50,0,50);
    //h_NbGoodvtx   = fs->make<TH1I>("NGoodVtx","NGoodVtx",50,0,50);
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
    //DYGenMass_nom   = fs->make<TH1F>("DYGenMass_nom","DYGenMass_nom", 45, Array);
    //DYGenMass_den   = fs->make<TH1F>("DYGenMass_den","DYGenMass_den", 45, Array);
//---------------------------------------------------------------------------------
//--------------------Histogramms for electron variables-------------------------------
//----------------------------------------------------------------------------------
    DiEle_Vrt    = fs->make<TH1F>("DiEleVrt","",30,0,30);

    DiEle_Mass_fromVrt = fs->make<TH1F>("DiEle_frVrt","",45,Array);
    El_DYMass_hist   = fs->make<TH1F>("El_DYMass","DYMass", 45, Array);
    El_DYZpkM_hist   = fs->make<TH1F>("El_DYMass_Zpk","DYMass_Zpk", 60, 60, 120);

    El_DYPt_hist     = fs->make<TH1F>("El_DYPt","El_DYPt", 100, 0., 250.);
    El_DYPt_lead     = fs->make<TH1F>("El_DYLeadPt","El_DYLeadPt", 100, 0., 250.);
    El_DYPt_sub      = fs->make<TH1F>("El_DYSubPt","El_DYSubPt", 100, 0., 250.);

    El_DYEta_hist    = fs->make<TH1F>("El_DYEta","El_DYEta", 100, -3, 3.);
    El_DYEta_lead    = fs->make<TH1F>("El_DYLeadEta","El_DYLeadEta", 100, -3, 3.);
    El_DYEta_sub     = fs->make<TH1F>("El_DYSubEta","El_DYSubEta", 100, -3, 3.);
    El_DYPhi_hist    = fs->make<TH1F>("El_DYPhi","El_DYPhi", 100, -4., 4.);
    El_DYPhi_lead    = fs->make<TH1F>("El_DYLeadPhi","El_DYLeadPhi", 100, -4., 4.);
    El_DYPhi_sub     = fs->make<TH1F>("El_DYSubPhi","El_DYSubPhi", 100, -4., 4.);

//--------------------Tree Branches--------------------------------------------------------
//--------------------Gen Part--------------------------------------------------------
   mytree->Branch("GenMass"          ,&GenMass            ,"GenMass/F");
   mytree->Branch("GenMass_den"      ,&GenMass_den        ,"GenMass_den/F");
   mytree->Branch("GenMass_nom"      ,&GenMass_nom        ,"GenMass_nom/F");
   mytree->Branch("genPdgId"       ,&genPdgId);
   mytree->Branch("genStatus"     ,&genStatus);
   mytree->Branch("genMotherId"   ,&genMotherId);
   mytree->Branch("genPt"         ,&genPt);
   mytree->Branch("genEta"        ,&genEta);
   mytree->Branch("genPhi"        ,&genPhi);
   mytree->Branch("genCharge"     ,&genCharge);
   mytree->Branch("genPx"         ,&genPx);
   mytree->Branch("genPy"         ,&genPy);
   mytree->Branch("genPz"         ,&genPz);
   mytree->Branch("genP"          ,&genP);
   mytree->Branch("fromHProcessFinalState",&fromHProcessFinalState_);
   mytree->Branch("fromHProcessDecayed"   ,&fromHProcessDecayed_);
   mytree->Branch("isHProcess"            ,&isHProcess_);
   mytree->Branch("genPreRap"          ,&genPreRap);
   mytree->Branch("genPreEn"          ,&genPreEn);
   mytree->Branch("genPtP4"          ,&genPtP4);
    mytree->Branch("genEtaP4"          ,&genEtaP4);
   mytree->Branch("genPhiP4"          ,&genPhiP4);
   mytree->Branch("genPxP4"          ,&genPxP4);
   mytree->Branch("genPyP4"          ,&genPyP4);
   mytree->Branch("genPzP4"          ,&genPzP4);
   mytree->Branch("genPP4"          ,&genPP4);
//---------------------------------------------------------------

//------------------ Weight------------------------------------------------------
   mytree->Branch("mcWeight",&mcWeight_,"mcWeight/F");
//-------------------Pile up-----------------------------------------------
   mytree->Branch("puTrue"      ,&puTrue_);
//------------------ Mass------------------------------------------------------
   mytree->Branch("mmM_Sel"          ,&mmM_Sel            ,"mmM_Sel/F");
   mytree->Branch("mmM_Vrt_Good"     ,&mmM_Vrt_Good       ,"mmM_Vrt_Good/F");
   mytree->Branch("mmM_Vrt_Sel"      ,&mmM_Vrt_Sel        ,"mmM_Vrt_Sel/F");
//--------------------Eta & Phi--------------------------------------------------------------
   mytree->Branch("muPt_Good"        ,&muPt_Good);
   mytree->Branch("muPt_Sel"         ,&muPt_Sel);
   mytree->Branch("muEta_Good"       ,&muEta_Good);
   mytree->Branch("muEta_Sel"        ,&muEta_Sel);
   mytree->Branch("muPhi_Good"       ,&muPhi_Good);
   mytree->Branch("muPhi_Sel"        ,&muPhi_Sel);
//--------------------Charge-----------------------------------------------
   mytree->Branch("muChrg_Good"      ,&muChrg_Good);
   mytree->Branch("muChrg_Sel"       ,&muChrg_Sel);
//--------------------Indexes of good muons---------------------------------------
   mytree->Branch("gdInx_mu1"      ,&gdInx_mu1        ,"gdInx_mu1/I");
   mytree->Branch("gdInx_mu2"      ,&gdInx_mu2        ,"gdInx_mu2/I");

//--------------------Addition--------------------------------------------------------------
   mytree->Branch("Angle_Good"   ,&Angle_Good,"Angle_Good/D");
   mytree->Branch("Angle_Sel"   ,&Angle_Sel,"Angle_Sel/D");
   mytree->Branch("vtxNormChi_Good" ,&vtxNormChi_Good,"vtxNormChi_Good/F");
   mytree->Branch("minVrt_Sel" ,&minVrt_Sel,"minVrt_Sel/F");
//--------------------Trigger---------------------------------------------------
   mytree->Branch("HLT_match"  ,&HLT_match);
   mytree->Branch("HLT_name", &HLT_name);
   mytree->Branch("HLT_isaccept", &HLT_isaccept);
   mytree->Branch("HLTObj_nbObj",&HLTObj_nbObj);
   mytree->Branch("HLTObj_pt",&HLTObj_pt);
   mytree->Branch("HLTObj_eta",&HLTObj_eta);
   mytree->Branch("HLTObj_phi",&HLTObj_phi);
   mytree->Branch("HLTObj_collection", &HLTObj_collection);
//------------------- Init Tree---------------------------------------------------
    mytree->Branch("event_runNo",  &Run,   "event_runNo/I");
    mytree->Branch("event_evtNo",  &Event, "event_evtNo/I");
    mytree->Branch("event_lumi",   &lumi,  "event_lumi/I");
    mytree->Branch("event_bunch",  &bunch, "event_bunch/I");

    mytree->Branch("nbPv_dl", &nbPv_dl ,"nbPv_dl/I");

    mytree->Branch("Mu_isTightMuon",&Mu_isTightMuon);
    mytree->Branch("Mu_isMediumMuon",&Mu_isMediumMuon);
    mytree->Branch("Mu_isHighPtMuon",&Mu_isHighPtMuon);
    mytree->Branch("Mu_isGlobalMuon",&Mu_isGlobalMuon);
    mytree->Branch("Mu_isTrackerMuon",&Mu_isTrackerMuon);

    mytree->Branch("muPtP4"       ,&muPtP4);
    mytree->Branch("muEtaP4"      ,&muEtaP4);
    mytree->Branch("muPhiP4"      ,&muPhiP4);
    mytree->Branch("muPxP4"       ,&muPxP4);
    mytree->Branch("muPyP4 "      ,&muPyP4);
    mytree->Branch("muPzP4"       ,&muPzP4);
    mytree->Branch("muPP4"        ,&muPP4);
    mytree->Branch("muChargP4"    ,&muChargP4);

    mytree->Branch("muPtTuneP"    ,&muPtTuneP);
    mytree->Branch("muEtaTuneP"   ,&muEtaTuneP);
    mytree->Branch("muPhiTuneP"   ,&muPhiTuneP);
    mytree->Branch("muPxTuneP"    ,&muPxTuneP);
    mytree->Branch("muPyTuneP"    ,&muPyTuneP);
    mytree->Branch("muPzTuneP"    ,&muPzTuneP);
    mytree->Branch("muPTuneP"     ,&muPTuneP);
    mytree->Branch("muChargTuneP" ,&muChargTuneP);

    mytree->Branch("muglobalMuonHits"  ,&muglobalMuonHits);
    mytree->Branch("munbMatchSt"       ,&munbMatchSt);
    mytree->Branch("mupixelHits"       ,&mupixelHits);
    mytree->Branch("mutrackerLayers"   ,&mutrackerLayers);
    mytree->Branch("mudxy"             ,&mudxy);
    mytree->Branch("mupterror"         ,&mupterror);
    mytree->Branch("muisoTrk"          ,&muisoTrk);
    mytree->Branch("muisoPf"          ,&muisoPf);
    mytree->Branch("mudz"              ,&mudz);
//dimuon value
    mytree->Branch("init_vtxNormChi"   ,&init_vtxNormChi);
    mytree->Branch("init_vrtMass"      ,&init_vrtMass);
//------------------- Electron Tree---------------------------------------------------
    mytree->Branch("elPt_Good"         ,&elPt_Good);
    mytree->Branch("elEta_Good"        ,&elEta_Good);
    mytree->Branch("elPhi_Good"        ,&elPhi_Good);
    mytree->Branch("elChrg_Good"       ,&elChrg_Good);
    mytree->Branch("elP_Good"          ,&elP_Good);
    mytree->Branch("elEnSC_Good"       ,&elEnSC_Good);
    mytree->Branch("elEtaSC_Good"      ,&elEtaSC_Good);
    mytree->Branch("elPhiSC_Good"      ,&elPhiSC_Good);

    mytree->Branch("elPt_Sel"         ,&elPt_Sel);
    mytree->Branch("elEta_Sel"        ,&elEta_Sel);
    mytree->Branch("elPhi_Sel"        ,&elPhi_Sel);
    mytree->Branch("elChrg_Sel"       ,&elChrg_Sel);

//------------------ Mass------------------------------------------------------
   mytree->Branch("eeM_Sel"          ,&eeM_Sel            ,"eeM_Sel/F");
   mytree->Branch("El_mmM_Vrt_Good"  ,&El_mmM_Vrt_Good    ,"El_mmM_Vrt_Good/F");
   mytree->Branch("El_mmM_Vrt_Sel"   ,&El_mmM_Vrt_Sel     ,"El_mmM_Vrt_Sel/F");
//------------------ Additional------------------------------------------------------
   mytree->Branch("El_Angle_Good"    ,&El_Angle_Good   ,"El_Angle_Good/D");
   mytree->Branch("El_Angle_Sel"     ,&El_Angle_Sel    ,"El_Angle_Sel/D");
   mytree->Branch("El_vtxNormChi_Good" ,&El_vtxNormChi_Good,"El_vtxNormChi_Good/F");
   mytree->Branch("El_minVrt_Sel"    ,&El_minVrt_Sel   ,"El_minVrt_Sel/F");
//--------------------Indexes of good electrons---------------------------------------
   mytree->Branch("gdInx_el1"      ,&gdInx_el1        ,"gdInx_el1/I");
   mytree->Branch("gdInx_el2"      ,&gdInx_el2        ,"gdInx_el2/I");
//--------Init Tree--------------------------------------
   mytree->Branch("elisMedium"   ,&elisMedium);

   mytree->Branch("elPtP4"       ,&elPtP4);
   mytree->Branch("elEtaP4"      ,&elEtaP4);
   mytree->Branch("elPhiP4"      ,&elPhiP4);
   mytree->Branch("elPxP4"       ,&elPxP4);
   mytree->Branch("elPyP4 "      ,&elPyP4);
   mytree->Branch("elPzP4"       ,&elPzP4);
   mytree->Branch("elPP4"        ,&elPP4);
   mytree->Branch("elChargP4"    ,&elChargP4);
        
   mytree->Branch("elEnSC"       ,&elEnSC);
   mytree->Branch("elEtaSC"      ,&elEtaSC);
   mytree->Branch("elPhiSC"      ,&elPhiSC);

   mytree->Branch("El_init_vtxNormChi"   ,&El_init_vtxNormChi);
   mytree->Branch("El_init_vrtMass"      ,&El_init_vrtMass);

 /*|-------------------------------------------------------|
   |                Tau Part                               |
   |-------------------------------------------------------|*/
   mytree->Branch("tauPt"     ,&tauPt);
   mytree->Branch("tauEta"    ,&tauEta);
   mytree->Branch("tauPhi"    ,&tauPhi);
   mytree->Branch("tauCharge" ,&tauCharge);
   mytree->Branch("tauMass"   ,&tauMass);
   mytree->Branch("taudxy"    ,&taudxy);
   mytree->Branch("TauIdIso"  ,&TauIdIso);
   mytree->Branch("TauIdMuRj" ,&TauIdMuRj);
   mytree->Branch("TauIdMVA"  ,&TauIdMVA);
   mytree->Branch("MassDYtaumuSS"   ,&MassDYtaumuSS);
   mytree->Branch("MassDYtaumuOS"   ,&MassDYtaumuOS);
   mytree->Branch("pfMET_"    ,&pfMET_    ,"pfMET_/F");
   mytree->Branch("pfMETPhi_" ,&pfMETPhi_ ,"pfMETPhi_/F");
   mytree->Branch("MuMetTranverseMass"   ,&MuMetTranverseMass   ,"MuMetTranverseMass/F");
 /*|-------------------------------------------------------|
   |                Photon Part                            |
   |-------------------------------------------------------|*/
   mytree->Branch("ptPhoton"     ,&ptPhoton_);
   mytree->Branch("etaPhoton"     ,&etaPhoton_);
   mytree->Branch("phiPhoton"     ,&phiPhoton_);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TauEl_Both_b::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

void 
TauEl_Both_b::beginRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when ending the processing of a run  ------------

void 
TauEl_Both_b::endRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when starting to processes a luminosity block  ------------

void 
TauEl_Both_b::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method called when ending the processing of a luminosity block  ------------

void 
TauEl_Both_b::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauEl_Both_b::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void TauEl_Both_b::clearTree() {
  /*|-------------------------------------------------------|
    |                Run Info Part                           |
    |-------------------------------------------------------|*/     
	Run = -999;
	Event = -999;
	lumi = -999;
	bunch = -999;
	isRealData = -999;
	nbPv_del.clear();
	nbPv_dl = -999;
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
        genPdgId.clear();
        genStatus.clear();
        genMotherId.clear();
        genPt.clear();
        genEta.clear();
        genPhi.clear();
        genCharge.clear();
        genPx.clear();
        genPy.clear();
        genPz.clear();
        genP.clear();
        fromHProcessFinalState_.clear();
        fromHProcessDecayed_.clear();
        isHProcess_.clear();
        genPreRap.clear();
        genPreEn.clear();
        genPtP4.clear();
        genEtaP4.clear();
        genPhiP4.clear();
        genPxP4.clear();
        genPyP4.clear();
        genPzP4.clear();
        genPP4.clear();
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
	gdInx_mu1 = -999;
	gdInx_mu2 = -999;
	muPt_Good.clear();
	muPt_Sel.clear();
	muEta_Good.clear();
	muEta_Sel.clear();
	muPhi_Good.clear();
	muPhi_Sel.clear();
	muChrg_Good.clear();
	muChrg_Sel.clear();
   //--------Init Tree--------------------------------------
     Mu_isTightMuon.clear();
     Mu_isMediumMuon.clear();
     Mu_isGlobalMuon.clear();
     Mu_isHighPtMuon.clear();
     Mu_isTrackerMuon.clear();
     muPtP4.clear();
     muEtaP4.clear();
     muPhiP4.clear();
     muPxP4.clear();
     muPyP4.clear();
     muPzP4.clear();
     muPP4.clear();
     muChargP4.clear();
     muPtTuneP.clear();
     muEtaTuneP.clear();
     muPhiTuneP.clear();
     muPxTuneP.clear();
     muPyTuneP.clear();
     muPzTuneP.clear();
     muPTuneP.clear();
     muChargTuneP.clear();
     muglobalMuonHits.clear();
     munbMatchSt.clear();
     mupixelHits.clear();
     mutrackerLayers.clear();
     mudxy.clear();
     mupterror.clear();
     muisoTrk.clear();
     muisoPf.clear();
     mudz.clear();
//dimuon value
     init_vtxNormChi.clear();
     init_vrtMass.clear();

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
  /*|-------------------------------------------------------|
    |                Trigger Part                           |
    |-------------------------------------------------------|*/
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
    |                Electron Part                          |
    |-------------------------------------------------------|*/
        elPt_Good.clear();
        elEta_Good.clear();
        elPhi_Good.clear();
        elChrg_Good.clear();
        elP_Good.clear();
        elEnSC_Good.clear();
        elEtaSC_Good.clear();
        elPhiSC_Good.clear();

        eeM_Sel         = -999;
	El_mmM_Vrt_Good    = -999; 
	El_mmM_Vrt_Sel     = -999; 
	El_Angle_Good      = -999;
	El_Angle_Sel       = -999;
	El_vtxNormChi_Good = 999;
	El_minVrt_Sel      = 21;
	gdInx_el1 = -999;
	gdInx_el2 = -999;

	elPt_Sel.clear();
	elEta_Sel.clear();
	elPhi_Sel.clear();
	elChrg_Sel.clear();
   //--------Init Tree--------------------------------------

        elisMedium.clear();

        elPtP4.clear();
        elEtaP4.clear();
        elPhiP4.clear();
        elPxP4.clear();
        elPyP4.clear();
        elPzP4.clear();
        elPP4.clear();
        elChargP4.clear();
        
        elEnSC.clear();
        elEtaSC.clear();
        elPhiSC.clear();

        El_init_vtxNormChi.clear();
        El_init_vrtMass.clear();
 /*|-------------------------------------------------------|
   |                Tau Part                               |
   |-------------------------------------------------------|*/
        tauPt.clear();
        tauEta.clear();
        tauPhi.clear();
        tauCharge.clear();
        tauMass.clear();
        taudxy.clear();
        TauIdIso.clear();
        TauIdMuRj.clear();
        TauIdMVA.clear();
        MassDYtaumuSS.clear();
        MassDYtaumuOS.clear();
 /*|-------------------------------------------------------|
   |                MET Part                               |
   |-------------------------------------------------------|*/
        pfMET_ = -999;
        pfMETPhi_ = -999;
        MuMetTranverseMass = -999;
        pfMETSumEt_ = -999;
 /*|-------------------------------------------------------|
   |                Photon Part                            |
   |-------------------------------------------------------|*/
        ptPhoton_.clear();
        etaPhoton_.clear();
        phiPhoton_.clear();
}
//define this as a plug-in
DEFINE_FWK_MODULE(TauEl_Both_b);
bool TauEl_Both_b::isPassHLT(){ 
  int nbMatch = 0;
  std::string pyName1,pyName2,pyName3;
  if (isMC_){ 
    pyName1 = "HLT_IsoMu20_v1";
    pyName2 = "HLT_IsoTkMu20_v1";
    pyName3 = "HLT_Ele22_eta2p1_WP75_Gsf_v1";
  } else {
    pyName1 = "HLT_IsoMu20_v3";
    pyName2 = "HLT_IsoTkMu20_v4";
    pyName3 = "HLT_Ele22_eta2p1_WPLoose_Gsf_v2";
  }
  //std::cout<<" HLT_nb.size() = "<< HLT_nb.size() << endl;
  for(unsigned i=0; i < HLT_nb.size(); i++){
    if( (HLT_name[i] == pyName1 || 
         HLT_name[i] == pyName2 ||
         HLT_name[i] == pyName3 ) && HLT_isaccept[i] == 1 ) {

//"HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v(.*)"
//"HLT_Ele22_eta2p1_WP75_Gsf_v(.*)"
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

bool TauEl_Both_b::RecoHLTMuonMatching(float RecoEta,float RecoPhi){
  int nbMatch  = 0;
  float deltaR = -10000.0;
  std::string pyName1,pyName2,pyName3;
  if (isMC_){ 
    pyName1 = "HLT_IsoMu20_v1";
    pyName2 = "HLT_IsoTkMu20_v1";
    pyName3 = "HLT_Ele22_eta2p1_WP75_Gsf_v1";
  } else {
    pyName1 = "HLT_IsoMu20_v3";
    pyName2 = "HLT_IsoTkMu20_v4";
    pyName3 = "HLT_Ele22_eta2p1_WPLoose_Gsf_v2";
  }
   //std::cout<<" HLTObj_nbObj = "<< HLTObj_nbObj.size() << endl;
	
  for(unsigned i=0; i<HLTObj_nbObj.size(); i++){
    if( HLTObj_collection[i] == pyName1 || 
		HLTObj_collection[i] == pyName2 ||
                  HLTObj_collection[i] == pyName3){ 
//|| HLT_name[i] == "HLT_Ele22_eta2p1_WP75_Gsf_v1"

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
float TauEl_Both_b::delR(float eta1,float phi1,float eta2,float phi2){
	float mpi=3.14;
	float dp=std::abs(phi1-phi2);
	if (dp>mpi) dp-=float(2*mpi);
	return sqrt((eta1-eta2)*(eta1-eta2) + dp*dp);
	}
float TauEl_Both_b::TMass_F(float pt3lep, float px3lep, float py3lep, float met, float metPhi) {
    return sqrt(pow(pt3lep + met, 2) - pow(px3lep + met * cos(metPhi), 2) - pow(py3lep + met * sin(metPhi), 2));
}
