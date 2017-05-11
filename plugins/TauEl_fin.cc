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

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "TH1.h"
//#include "MSETools.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;
//
// class declaration
//

class TauEl_fin : public edm::EDAnalyzer {
   public:
      explicit TauEl_fin(const edm::ParameterSet&);
      ~TauEl_fin();

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
      virtual void endRun(edm::Run const& iRun, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      void clearTree();
      inline double getEffectiveAreaForElectrons(const double& ) const;//eta
      const reco::Candidate* mother(const reco::CandidateBaseRef& cand);
      int motherId(const reco::CandidateBaseRef& cand);
      int grandmotherId(const reco::CandidateBaseRef& cand);
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
      edm::EDGetTokenT<LHEEventProduct> lheEvtInfo_;
      edm::EDGetTokenT<LHERunInfoProduct> lheInfo_;

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
       float SEt;
       float SetaSC;
       float SphiSC;
       bool isMedium;
       bool STauIdIso;
       bool STauIdMuRj;
       bool STauIdMVA;
       float Smass;
       bool isMedium2;
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
         int mother0Id;
         int grmotherZId;
         int grmotherN;

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
       float p;
      };
      static bool lepSortingRuleGen(GENPARTICLE x, GENPARTICLE y) {return x.p4.Pt() > y.p4.Pt();}
  /*|-------------------------------------------------------|
    |                Histogramms Part                       |
    |-------------------------------------------------------|*/
      TH1D  *hWEvents_;
      TH1I  *pdgID_;
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
      std::vector<double> lheWts;
      double lheWtsMu;
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
    |                Gen Part                               |
    |-------------------------------------------------------|*/
      float  mcWeight_;
      float GenMass;
      float GenMass_den; 
      float GenMass_nom; 
      std::vector<int>    genPdgId;
      std::vector<int>    genStatus;
      std::vector<int>    genMotherId;
      std::vector<int>    mother0Id;
      std::vector<int>    genGrandMotherId;
      std::vector<int>    genGrandMotherN;
      std::vector<float>  genPt;
      std::vector<float>  genEta;
      std::vector<float>  genPhi;
      std::vector<int>    genCharge;
      std::vector<float>  genPx;
      std::vector<float>  genPy;
      std::vector<float>  genPz;
      std::vector<float>  genP;

      std::vector<int>    genPdgId_my;
      std::vector<int>    genStatus_my;
      std::vector<int>    genMotherId_my;
      std::vector<float>  genPt_my;
      std::vector<float>  genEta_my;
      std::vector<float>  genPhi_my;
      std::vector<int>    genCharge_my;
      std::vector<float>  genPx_my;
      std::vector<float>  genPy_my;
      std::vector<float>  genPz_my;
      std::vector<float>  genP_my;
      std::vector<float>  genPreRap_my;
      std::vector<float>  genPreEn_my;
      std::vector<float>  genPtP4_my;
      std::vector<float>  genEtaP4_my;
      std::vector<float>  genPhiP4_my;
      std::vector<float>  genPxP4_my;
      std::vector<float>  genPyP4_my;
      std::vector<float>  genPzP4_my;
      std::vector<float>  genPP4_my;

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
      float  emuM_Vrt_Sel; 
      double Angle_Good;
      double Angle_Sel;
      float  vtxNormChi_Good;
      float  minVrt_Sel;
      float  minVrt_Sel_emu;
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
     std::vector<bool> Mu_isLooseMuon;
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
        Int_t                 nPU_ch;
        Int_t                 puTrue_ch;
  /*|-------------------------------------------------------|
    |                PrimaryVertex Part                     |
    |-------------------------------------------------------|*/
	std::vector<int>    nbPv;
  	std::vector<int>    Nbdof;
        std::vector<float>  nbPv_del;
        Int_t nbPv_dl; // number of reconsrtucted primary vertices 
        //nbPv_dl;
        int                   nPV;
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

	std::vector<float>    elEt;
	std::vector<float>    elPt;
	std::vector<float>    elEta;
	std::vector<float>    elPhi;
	std::vector<float>    elPx;
	std::vector<float>    elPy;
	std::vector<float>    elPz;
	std::vector<float>    elP;
	std::vector<int>      elCharg;

        
	std::vector<float>    elEnSC;
	std::vector<float>    elEtaSC;
	std::vector<float>    elPhiSC;

	std::vector<float>    El_init_vtxNormChi;
	std::vector<float>    El_init_vrtMass;
	std::vector<float>    emu_commonVrt;
        std::vector<float>    isoDeltaBeta_;

        std::vector<float> r9_el;
        std::vector<float> isoPFUnc_el;
        std::vector<float> isoPFRho_el;
        std::vector<float> sigmaIEtaIEta_el;
        std::vector<float> hadronicOverEm_el;
        std::vector<bool> hasMatchedConversion_el;
        std::vector<int> MissedInnerTrackHits_el;
        std::vector<int> EleInnerTrackHits_el;
        std::vector<int> EleLostInnerTrackHits_el;
        std::vector<float> gsfTrackPt_el;
        std::vector<float> SCRawE_el;
        std::vector<float> eSCOverP_el; 
        std::vector<float> dPhiSCTrackAtVtx_el;
        std::vector<float> dEtaSCTrackAtVtx_el; 
        std::vector<float> epDiff_el; 
        std::vector<float> d0_el;
        std::vector<float> dz_el;
        std::vector<bool> isMediumOwn_El;
        float rho_;

        std::vector<bool> isEE_el;
        std::vector<bool> isEB_el;

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
TauEl_fin::TauEl_fin(const edm::ParameterSet& iConfig):
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
   lheEvtInfo_(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEvtInfo"))),
   lheInfo_(mayConsume<LHERunInfoProduct,InRun>(iConfig.getParameter<edm::InputTag>("externalLHEProducer"))),
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


TauEl_fin::~TauEl_fin()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TauEl_fin::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  for( auto & puInfoElement : *genPileupHandle){
     if( puInfoElement.getBunchCrossing() == 0 ){
	nPU_ch     = puInfoElement.getPU_NumInteractions();
	puTrue_ch  = puInfoElement.getTrueNumInteractions();
      }
    }
 } 
/*|-------------------------------------------------------|
  |                PrimaryVertex Part                     |
  |-------------------------------------------------------|*/
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  if (vertices->empty()) return;
  const reco::Vertex &PV = vertices->front();
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
  
  edm::Handle<LHEEventProduct> lhEvtInfo;
  iEvent.getByToken( lheEvtInfo_, lhEvtInfo);
  //const LHEEventProduct* Product = lhEvtInfo.product();

    //double lheweightMu = 1;
    for(unsigned int i = 0; i < lhEvtInfo->weights().size(); i++) {
        double lheweight = lhEvtInfo->weights()[i].wgt / lhEvtInfo->originalXWGTUP();
        lheWtsMu *= lhEvtInfo->weights()[i].wgt / lhEvtInfo->originalXWGTUP();
        lheWts.push_back(lheweight); 
        //lheWtsMu = lheweightMu; 
        //lheWts_wgt.push_back(lhEvtInfo->weights()[i].wgt); 
        //lheWts_origUp.push_back(lhEvtInfo->originalXWGTUP()); 
    }

   for (const reco::GenParticle &i_gen : *genPrunedColl) {           
      // ----------- Muons and Electrons -------------
      // set lepton 4-vector and ID
      fromHProcessFinalState_.push_back(i_gen.fromHardProcessFinalState());
      fromHProcessDecayed_.push_back(i_gen.fromHardProcessDecayed());
      isHProcess_.push_back(i_gen.isHardProcess());
      //mother0Id.push_back(i_gen.mother(0)->pdgId());

      GENPARTICLE aGenLepton;
      TLorentzVector lepP4GEN(i_gen.p4().Px(),i_gen.p4().Py(),i_gen.p4().Pz(),i_gen.p4().E());
      aGenLepton.p4 = lepP4GEN;
      aGenLepton.pdgId = i_gen.pdgId();
      aGenLepton.charge = i_gen.charge();
      aGenLepton.status = i_gen.status();
      aGenLepton.rap = i_gen.rapidity();
      aGenLepton.en  = i_gen.energy();
      aGenLepton.pt  = i_gen.pt();
      aGenLepton.eta = i_gen.eta();
      aGenLepton.phi = i_gen.phi();
      aGenLepton.px  = i_gen.px();
      aGenLepton.py  = i_gen.py();
      aGenLepton.pz  = i_gen.pz();
      aGenLepton.p  = i_gen.p();

      //filling initial tree variables
      genPdgId.push_back(aGenLepton.pdgId);
      genStatus.push_back(aGenLepton.status);
      genPtP4.push_back(aGenLepton.p4.Pt());
      //genEtaP4.push_back(aGenLepton.p4.Eta());
      genPhiP4.push_back(aGenLepton.p4.Phi());
      genCharge.push_back(aGenLepton.charge);
      genPxP4.push_back(aGenLepton.p4.Px());
      genPyP4.push_back(aGenLepton.p4.Py());
      genPzP4.push_back(aGenLepton.p4.Pz());
      genPP4.push_back(aGenLepton.p4.P());
      genPreRap.push_back(aGenLepton.rap);
      genPreEn.push_back(aGenLepton.en);
      genPx.push_back(aGenLepton.px);
      genPy.push_back(aGenLepton.py);
      genPz.push_back(aGenLepton.pz);
      genP.push_back(aGenLepton.p);
      genPt.push_back(aGenLepton.pt);
      genEta.push_back(aGenLepton.eta);
      genPhi.push_back(aGenLepton.phi);

      int pId = i_gen.pdgId();
      const reco::Candidate* mom = i_gen.mother();
      while (mom != 0 && mom->pdgId() == pId)
        mom = mom->mother();
      if (mom != 0) 
      aGenLepton.motherId = mom->pdgId();
      else aGenLepton.motherId = -999;

      if (mom == 0){ 
         aGenLepton.grmotherZId = -999;
         aGenLepton.grmotherN = -999;
      }
      else {
      const reco::Candidate *gm = mom->mother();
      if (gm != 0) aGenLepton.grmotherN = gm->pdgId();
      else aGenLepton.grmotherN = -999;


      while (gm != 0 && gm->pdgId() == mom->pdgId())
         gm = gm->mother();
      if (gm != 0) aGenLepton.grmotherZId = gm->pdgId();
      else aGenLepton.grmotherZId = -999;
      }

        genMotherId.push_back(aGenLepton.motherId);
        genGrandMotherId.push_back(aGenLepton.grmotherZId);
        genGrandMotherN.push_back(aGenLepton.grmotherN);

      //aGenLepton.mother0Id = -999;
      //aGenLepton.motherId = -999;
      //aGenLepton.grmotherZId = -999;




/*
      if (i_gen.numberOfMothers() != 0){
      const reco::Candidate * motherCand = i_gen.mother(0);
      while(motherCand!=0){
        //if(abs(i_gen.pdgId()) != 13) continue;
        aGenLepton.mother0Id = i_gen.mother(0)->pdgId();
        aGenLepton.motherId = i_gen.mother()->pdgId();
      }  
        //if (abs(aGenLepton.motherId) == 23 || abs(aGenLepton.motherId) ==22 || abs(aGenLepton.motherId) ==21 || abs(aGenLepton.motherId) ==9)
          aGenLepton.grmotherZId = i_gen.mother()->mother()->pdgId();
      }
        mother0Id.push_back(aGenLepton.mother0Id);
        genMotherId.push_back(aGenLepton.motherId);
        genMotherIdZph.push_back(aGenLepton.grmotherZId);
*/
      //if ( abs(i_gen.pdgId()) != 13 && abs(i_gen.pdgId()) != 11) continue;
      if ( i_gen.status() != 1 ) continue;

      //if (abs(i_gen.pdgId()) == 15) cout << "i_gen.pdgId() = " << i_gen.pdgId()<<endl;
      // Mother particles are of type reco::GenParticle


      const reco::Candidate * motherCand = i_gen.mother(0);
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
     genPdgId_my.push_back(myGenLeptons[l].pdgId);
     genStatus_my.push_back(myGenLeptons[l].status);
     genMotherId_my.push_back(myGenLeptons[l].motherId);
     genPtP4_my.push_back(myGenLeptons[l].p4.Pt());
     genEtaP4_my.push_back(myGenLeptons[l].p4.Eta());
     genPhiP4_my.push_back(myGenLeptons[l].p4.Phi());
     genCharge_my.push_back(myGenLeptons[l].charge);
     genPxP4_my.push_back(myGenLeptons[l].p4.Px());
     genPyP4_my.push_back(myGenLeptons[l].p4.Py());
     genPzP4_my.push_back(myGenLeptons[l].p4.Pz());
     genPP4_my.push_back(myGenLeptons[l].p4.P());
     genPreRap_my.push_back(myGenLeptons[l].rap);
     genPreEn_my.push_back(myGenLeptons[l].en);
     genPx_my.push_back(myGenLeptons[l].px);
     genPy_my.push_back(myGenLeptons[l].py);
     genPz_my.push_back(myGenLeptons[l].pz);
     genPt_my.push_back(myGenLeptons[l].pt);
     genEta_my.push_back(myGenLeptons[l].eta);
     genPhi_my.push_back(myGenLeptons[l].phi);
  }

}

    Handle<double> rho;
    iEvent.getByLabel(mSrcRho,rho);
    rho_ = *rho;

    edm::Handle<reco::BeamSpot> beamspot;
    iEvent.getByToken(beamSpotToken_,beamspot);
    //const reco::BeamSpot beamSpot = (*beamspot);
    edm::Handle<reco::ConversionCollection> conversions;
    iEvent.getByToken(conversionsToken_, conversions);

  /*|-------------------------------------------------------|
    |                Electron Part                          |
    |-------------------------------------------------------|*/

    edm::Handle<edm::View<reco::GsfElectron> > electrons;
    //edm::Handle<reco::GsfElectronCoreCollection> electrons;
    iEvent.getByToken(electronsToken_, electrons);

    edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
    iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);

    vector<PATRONUS> InitEl;
    //vector<int> passMediumId_;

    for (size_t i = 0; i < electrons->size(); ++i){
    const auto el = electrons->ptrAt(i);
    PATRONUS elON;
    TLorentzVector elONP4(el->p4().Px(), el->p4().Py(), el->p4().Pz(), el->p4().E() );
    elON.p4      = elONP4;
    elON.Jcharge = el->charge();
    elON.Sp      = el->p();
    elON.Spx     = el->px();
    elON.Spy     = el->py();
    elON.Spz     = el->pz();
    elON.Spt     = el->pt();
    elON.Seta    = el->eta();
    elON.Sphi    = el->phi();
    elON.SEt     = el->superCluster()->energy()/cosh(el->superCluster()->eta());
    elON.SEn     = el->superCluster()->energy();
    elON.SetaSC  = el->superCluster()->eta();
    elON.SphiSC  = el->superCluster()->phi();

    edm::ESHandle<TransientTrackBuilder> ttkb;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb);
    reco::TransientTrack ttv;
    reco::GsfTrackRef theTrack = el->gsfTrack();
    ttv = ttkb->build(theTrack);
    elON.Strack  = ttv;

    bool isPassMedium = (*medium_id_decisions)[el];
    elON.isMedium = isPassMedium;

// add a selection for muon
    GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
    isoDeltaBeta_.push_back((pfIso.sumChargedHadronPt + max<float>( 0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt))/(el->pt()));

// Compute isolation with delta beta correction for PU
     float isoChargedHadrons = pfIso.sumChargedHadronPt;
     float isoNeutralHadrons = pfIso.sumNeutralHadronEt;
     float isoPhotons        = pfIso.sumPhotonEt;
     //float isoChargedFromPU  = pfIso.sumPUPt;

     float electronIsoPFUnc   = (isoChargedHadrons + isoNeutralHadrons + isoPhotons)/el->pt();
     float electronIsoPFRho   = (isoChargedHadrons + std::max(isoNeutralHadrons + isoPhotons -  (*rho) * getEffectiveAreaForElectrons(el->eta()), 0.))/el->pt();

     float sigmaIetaIeta                  = el->full5x5_sigmaIetaIeta();
     float hadronicOverEm                 = el->hcalOverEcal(); //hOverE
     //float hadronicOverEm                 = el->hadronicOverEm();
     float dEtaSCTrackAtVtx               = el->deltaEtaSuperClusterTrackAtVtx(); //dEtaIn
     float dPhiSCTrackAtVtx               = el->deltaPhiSuperClusterTrackAtVtx(); //dPhiIn
     float epDifference                   = fabs( 1./el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()  ); //ooEmooP
     int   expectedMissingInnerHits       = el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS); //expectedMissingInnerHits
     float d0                             = (-1) * theTrack->dxy(nGoodVtxs->position() );
     float dz                             = theTrack->dz( nGoodVtxs->position() );
     bool hasMConv(false);
     //if(ConversionTools::hasMatchedConversion(*i_el, conversions, beamspot.position())) {
     const reco::GsfElectron* recoelectron=electrons->refAt(i).get();
     const pat::Electron* patelec=dynamic_cast<const pat::Electron*>(recoelectron);
     if (patelec->passConversionVeto()){
         hasMConv = true;
     }

     r9_el.push_back(el->r9());
     isoPFUnc_el.push_back(electronIsoPFUnc);
     isoPFRho_el.push_back(electronIsoPFRho);
     sigmaIEtaIEta_el.push_back(sigmaIetaIeta);
     hadronicOverEm_el.push_back(hadronicOverEm);
    
     bool matchConv = 0;
     if (hasMConv) {
       matchConv=1;
       }
     hasMatchedConversion_el.push_back(matchConv);
     MissedInnerTrackHits_el.push_back(expectedMissingInnerHits);
     EleInnerTrackHits_el.push_back(el->gsfTrack()->numberOfValidHits());
     EleLostInnerTrackHits_el.push_back(el->gsfTrack()->numberOfLostHits());
     gsfTrackPt_el.push_back(el->gsfTrack()->pt());
     SCRawE_el.push_back(el->superCluster()->rawEnergy());
     eSCOverP_el.push_back(el-> eSuperClusterOverP()); 
     dPhiSCTrackAtVtx_el.push_back(el->deltaPhiSuperClusterTrackAtVtx());
     dEtaSCTrackAtVtx_el.push_back(el->deltaEtaSuperClusterTrackAtVtx()); 
     epDiff_el.push_back(1./el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()  ); 

     d0_el.push_back(d0);
     dz_el.push_back(dz);
     isEE_el.push_back(el->isEE());
     isEB_el.push_back(el->isEB());

     InitEl.push_back(elON);

      bool isMediumOwn(false);
      elON.isMedium2 = isMediumOwn;
      if (el->isEB()) {

                if (sigmaIetaIeta            > 0.0101) continue;
                if (dEtaSCTrackAtVtx         > 0.0103) continue;
                if (dPhiSCTrackAtVtx         > 0.0336) continue; 
                if (hadronicOverEm           > 0.0876) continue;
                if (electronIsoPFRho         > 0.0766) continue;
                if (epDifference             > 0.0174) continue;
                if (fabs(d0)                 > 0.0118) continue;
                if (fabs(dz)                 > 0.373)  continue;
                if (expectedMissingInnerHits >= 2 )    continue;
                if (!hasMConv)                         continue;
                isMediumOwn = true;
  
         }
      // Endcap Electron ID
        if (el->isEE()) {
                if (sigmaIetaIeta            > 0.0283)  continue;
                if (dEtaSCTrackAtVtx         > 0.00733) continue;
                if (dPhiSCTrackAtVtx         > 0.114)   continue; 
                if (hadronicOverEm           > 0.0678)  continue;
                if (electronIsoPFRho         > 0.0678)  continue;
                if (epDifference             > 0.0898)  continue;
                if (fabs(d0)                 > 0.0739)  continue;
                if (fabs(dz)                 > 0.602)   continue;
                if (expectedMissingInnerHits >= 1 )     continue;
                if (!hasMConv)                          continue;
                isMediumOwn = true;
            
          }
    elON.isMedium2 = isMediumOwn;
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

     elEnSC.push_back(InitEl[l].SEn);
     elEtaSC.push_back(InitEl[l].SetaSC);
     elPhiSC.push_back(InitEl[l].SphiSC);
     elisMedium.push_back(InitEl[l].isMedium);
     isMediumOwn_El.push_back(InitEl[l].isMedium2);

     elEt.push_back(InitEl[l].SEt);
     elPt.push_back(InitEl[l].Spt);
     elEta.push_back(InitEl[l].Seta);
     elPhi.push_back(InitEl[l].Sphi);
     elPx.push_back(InitEl[l].Spx);
     elPy.push_back(InitEl[l].Spy);
     elPz.push_back(InitEl[l].Spz);
     elP.push_back(InitEl[l].Sp);
     elCharg.push_back(InitEl[l].Jcharge);

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

  vector<PATRONUS> InitMuons;
  const reco::Vertex &PV_na = vertices->front();

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  for (const pat::Muon &mu : *muons) {

    Mu_isMediumMuon.push_back(mu.isMediumMuon());
    Mu_isLooseMuon.push_back(mu.isLooseMuon());
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
    }
    if (trackerTrack.isNonnull()){
         const reco::HitPattern & inhit = trackerTrack->hitPattern();
         muON.pixelHits     = inhit.numberOfValidPixelHits();
         muON.trackerLayers = inhit.trackerLayersWithMeasurement();
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
  //cout<<"nbInMu =" << nbInMu << endl;
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
         //if (vtxNormChi < 20) cout<<"vtxNormChi =" << vtxNormChi <<endl;
         //cout<<"vtxNormChi =" << vtxNormChi <<endl;
         Measurement1D massK = imfv.invariantMass(vtx, muon_mass);
         init_vrtMass.push_back(massK.value());
       }
   }
// Emu part, filling common vertex tree
  for(unsigned l = 0; l < InitMuons.size(); l++) {
    for(unsigned s = 0; s < InitEl.size(); s++) {
         std::vector<reco::TransientTrack> ttV_em;
         ttV_em.clear();
         ttV_em.push_back(InitMuons[l].Strack);
         ttV_em.push_back(InitEl[s].Strack);
         KalmanVertexFitter kvf(true);
         CachingVertex<5> vtx_em = kvf.vertex(ttV_em);
         float vtxNormChi = vtx_em.totalChiSquared()/vtx_em.degreesOfFreedom();
         emu_commonVrt.push_back(vtxNormChi);
       }
   }

 /*|-------------------------------------------------------|
   |                MET Part                               |
   |-------------------------------------------------------|*/

  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  const pat::MET &met = mets->front();
    pfMET_       = met.et();
    pfMETPhi_    = met.phi();
    pfMETSumEt_  = met.sumEt();

 /*|-------------------------------------------------------|
   |                Tau Part                               |
   |-------------------------------------------------------|*/

  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauToken_, taus);
  vector<PATRONUS> InitTau;

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
TauEl_fin::beginJob()
{
/*
    Float_t Array[46] = {15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,
    64.,68.,72.,76.,81.,86.,91.,96.,101.,106.,110.,115.,
    120.,126.,133.,141.,150.,160.,171.,185.,
    200.,220.,243.,273.,320.,380.,440.,510.,
    600.,700.,830.,1000.,1200.,1500.,2000.,3000.};
*/
    edm::Service<TFileService> fs;
    mytree = fs->make<TTree>("tree","tr");
//------------------ Weight------------------------------------------------------
    hWEvents_   = fs->make<TH1D>("WEvents", "Weighted Events",6,0,6);
    hWEvents_->Sumw2();
//------------------ Vertex------------------------------------------------------
    h_NbVtx       = fs->make<TH1I>("NVtx","NVtx",50,0,50);

//*****************************************************************************************
//--------------------Tree Branches--------------------------------------------------------
//*****************************************************************************************

   mytree->Branch("event_runNo",  &Run,     "event_runNo/I");
   mytree->Branch("event_evtNo",  &Event,   "event_evtNo/I");
   mytree->Branch("event_lumi",   &lumi,    "event_lumi/I");
   mytree->Branch("event_bunch",  &bunch,   "event_bunch/I");
   mytree->Branch("rho",          &rho_,    "rho/F");
   mytree->Branch("nbPv_dl",      &nbPv_dl, "nbPv_dl/I");


//--------------------Gen Part--------------------------------------------------------

   mytree->Branch("mcWeight",&mcWeight_,"mcWeight/F");
   mytree->Branch("lheWts",      &lheWts);
   mytree->Branch("lheWtsMu",      &lheWtsMu,"lheWtsMu/D");

   mytree->Branch("fromHProcessDecayed"   ,&fromHProcessDecayed_);
   mytree->Branch("fromHProcessFinalState",&fromHProcessFinalState_);
   mytree->Branch("isHProcess"            ,&isHProcess_);
   mytree->Branch("genPdgId"              ,&genPdgId);
   mytree->Branch("genStatus"             ,&genStatus);
   mytree->Branch("genMotherId"           ,&genMotherId);
   mytree->Branch("mother0Id"             ,&mother0Id);
   mytree->Branch("genGrandMotherId"        ,&genGrandMotherId);
   mytree->Branch("genGrandMotherN"        ,&genGrandMotherN);
   mytree->Branch("genPtP4"               ,&genPtP4);
   mytree->Branch("genEtaP4"              ,&genEtaP4);
   mytree->Branch("genPhiP4"              ,&genPhiP4);
   mytree->Branch("genPxP4"               ,&genPxP4);
   mytree->Branch("genPyP4"               ,&genPyP4);
   mytree->Branch("genPzP4"               ,&genPzP4);
   mytree->Branch("genPP4"                ,&genPP4);
   mytree->Branch("genPreRap"             ,&genPreRap);
   mytree->Branch("genPreEn"              ,&genPreEn);
   mytree->Branch("genPt"                 ,&genPt);
   mytree->Branch("genEta"                ,&genEta);
   mytree->Branch("genPhi"                ,&genPhi);
   mytree->Branch("genCharge"             ,&genCharge);
   mytree->Branch("genPx"                 ,&genPx);
   mytree->Branch("genPy"                 ,&genPy);
   mytree->Branch("genPz"                 ,&genPz);
   mytree->Branch("genP"                  ,&genP);

   mytree->Branch("genPdgId_my"              ,&genPdgId_my);
   mytree->Branch("genStatus_my"             ,&genStatus_my);
   mytree->Branch("genMotherId_my"           ,&genMotherId_my);
   mytree->Branch("genPtP4_my"               ,&genPtP4_my);
   mytree->Branch("genEtaP4_my"              ,&genEtaP4_my);
   mytree->Branch("genPhiP4_my"              ,&genPhiP4_my);
   mytree->Branch("genPxP4_my"               ,&genPxP4_my);
   mytree->Branch("genPyP4_my"               ,&genPyP4_my);
   mytree->Branch("genPzP4_my"               ,&genPzP4_my);
   mytree->Branch("genPP4_my"                ,&genPP4_my);
   mytree->Branch("genPreRap_my"             ,&genPreRap_my);
   mytree->Branch("genPreEn_my"              ,&genPreEn_my);
   mytree->Branch("genPt_my"                 ,&genPt_my);
   mytree->Branch("genEta_my"                ,&genEta_my);
   mytree->Branch("genPhi_my"                ,&genPhi_my);
   mytree->Branch("genCharge_my"             ,&genCharge_my);
   mytree->Branch("genPx_my"                 ,&genPx_my);
   mytree->Branch("genPy_my"                 ,&genPy_my);
   mytree->Branch("genPz_my"                 ,&genPz_my);
   mytree->Branch("genP_my"                  ,&genP_my);

//-------------------Pile up-----------------------------------------------

   mytree->Branch("nPU_ch",   &nPU_ch,"nPU_ch/I");
   mytree->Branch("puTrue_ch",&puTrue_ch,"puTrue_ch/I");

//--------------------Trigger---------------------------------------------------

   mytree->Branch("HLT_name",          &HLT_name);
   mytree->Branch("HLT_isaccept",      &HLT_isaccept);
   mytree->Branch("HLTObj_nbObj",      &HLTObj_nbObj);
   mytree->Branch("HLTObj_pt",         &HLTObj_pt);
   mytree->Branch("HLTObj_eta",        &HLTObj_eta);
   mytree->Branch("HLTObj_phi",        &HLTObj_phi);
   mytree->Branch("HLTObj_collection", &HLTObj_collection);

//--------------------Electrons---------------------------------------------------

   mytree->Branch("elPtP4",         &elPtP4);
   mytree->Branch("elEtaP4",        &elEtaP4);
   mytree->Branch("elPhiP4",        &elPhiP4);
   mytree->Branch("elPxP4",         &elPxP4);
   mytree->Branch("elPyP4",         &elPyP4);
   mytree->Branch("elPzP4",         &elPzP4);
   mytree->Branch("elPP4",          &elPP4);
   mytree->Branch("elEnSC",         &elEnSC);
   mytree->Branch("elEtaSC",        &elEtaSC);
   mytree->Branch("elPhiSC",        &elPhiSC);
   mytree->Branch("elisMedium",     &elisMedium);
   mytree->Branch("isMediumOwn_El", &isMediumOwn_El);
   mytree->Branch("elEt",           &elEt);
   mytree->Branch("elPt",           &elPt);
   mytree->Branch("elEta",          &elEta);
   mytree->Branch("elPhi",          &elPhi);
   mytree->Branch("elPx",           &elPx);
   mytree->Branch("elPy",           &elPy);
   mytree->Branch("elPz",           &elPz);
   mytree->Branch("elP",            &elP);
   mytree->Branch("elCharg",        &elCharg);

   mytree->Branch("El_init_vtxNormChi",  &El_init_vtxNormChi);
   mytree->Branch("El_init_vrtMass",     &El_init_vrtMass);
        
   mytree->Branch("isoDeltaBeta", &isoDeltaBeta_);

   mytree->Branch("r9_el",                     &r9_el);
   mytree->Branch("isoPFUnc_el",               &isoPFUnc_el);
   mytree->Branch("isoPFRho_el",               &isoPFRho_el);
   mytree->Branch("sigmaIEtaIEta_el",          &sigmaIEtaIEta_el);
   mytree->Branch("hadronicOverEm_el",         &hadronicOverEm_el);
   mytree->Branch("hasMatchedConversion_el",   &hasMatchedConversion_el);
   mytree->Branch("MissedInnerTrackHits_el",   &MissedInnerTrackHits_el);
   mytree->Branch("EleInnerTrackHits_el",      &EleInnerTrackHits_el);
   mytree->Branch("EleLostInnerTrackHits_el",  &EleLostInnerTrackHits_el);
   mytree->Branch("gsfTrackPt_el",             &gsfTrackPt_el);
   mytree->Branch("SCRawE_el",                 &SCRawE_el);
   mytree->Branch("eSCOverP_el",               &eSCOverP_el); 
   mytree->Branch("dPhiSCTrackAtVtx_el",       &dPhiSCTrackAtVtx_el);
   mytree->Branch("dEtaSCTrackAtVtx_el",       &dEtaSCTrackAtVtx_el); 
   mytree->Branch("epDiff_el",                 &epDiff_el); 
   mytree->Branch("d0_el",                     &d0_el);
   mytree->Branch("dz_el",                     &dz_el);
   mytree->Branch("isEB_el",                   &isEB_el);
   mytree->Branch("isEE_el",                   &isEE_el);

//--------------------Muons---------------------------------------------------

    mytree->Branch("Mu_isTightMuon",  &Mu_isTightMuon);
    mytree->Branch("Mu_isMediumMuon", &Mu_isMediumMuon);
    mytree->Branch("Mu_isHighPtMuon", &Mu_isHighPtMuon);
    mytree->Branch("Mu_isGlobalMuon", &Mu_isGlobalMuon);
    mytree->Branch("Mu_isTrackerMuon", &Mu_isTrackerMuon);
    mytree->Branch("Mu_isLooseMuon",  &Mu_isLooseMuon);

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
    mytree->Branch("muisoPf"           ,&muisoPf);
    mytree->Branch("mudz"              ,&mudz);
//dimuon value
    mytree->Branch("init_vtxNormChi"   ,&init_vtxNormChi);
    mytree->Branch("init_vrtMass"      ,&init_vrtMass);
    mytree->Branch("emu_commonVrt"     ,&emu_commonVrt);

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

 /*|-------------------------------------------------------|
   |                MET Part                               |
   |-------------------------------------------------------|*/
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
TauEl_fin::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

void 
TauEl_fin::beginRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when ending the processing of a run  ------------

void 
TauEl_fin::endRun(edm::Run const& iRun, edm::EventSetup const&)
{
    edm::Handle<LHERunInfoProduct> run;
    typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
    iRun.getByToken(lheInfo_,run);
//    iRun.getByLabel( "externalLHEProducer", run );
    LHERunInfoProduct myLHERunInfoProduct = *(run.product());

    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
      //std::cout << iter->tag() << std::endl;
      std::vector<std::string> lines = iter->lines();
      for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
      //std::cout << lines.at(iLine);
      }
    }
}


// ------------ method called when starting to processes a luminosity block  ------------

void 
TauEl_fin::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method called when ending the processing of a luminosity block  ------------

void 
TauEl_fin::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauEl_fin::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void TauEl_fin::clearTree() {
  /*|-------------------------------------------------------|
    |                Run Info Part                           |
    |-------------------------------------------------------|*/     
	Run = -999;
	Event = -999;
	lumi = -999;
	bunch = -999;
	isRealData = -999;
	nbPv_dl = -999;
        rho_ = -999;
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
  /*|-------------------------------------------------------|
    |                Gen Part                           |
    |-------------------------------------------------------|*/
	mcWeight_   = -999;
        lheWts.clear();
        lheWtsMu = -999;

        fromHProcessFinalState_.clear();
        fromHProcessDecayed_.clear();
        isHProcess_.clear();
        genPdgId.clear();
        genStatus.clear();
        genMotherId.clear();
        mother0Id.clear();
        genGrandMotherId.clear();
        genGrandMotherN.clear();
        genPtP4.clear();
        genEtaP4.clear();
        genPhiP4.clear();
        genPxP4.clear();
        genPyP4.clear();
        genPzP4.clear();
        genPP4.clear();
        genPreRap.clear();
        genPreEn.clear();
        genPt.clear();
        genEta.clear();
        genPhi.clear();
        genCharge.clear();
        genPx.clear();
        genPy.clear();
        genPz.clear();
        genP.clear();

        genPdgId_my.clear();
        genStatus_my.clear();
        genMotherId_my.clear();
        genPtP4_my.clear();
        genEtaP4_my.clear();
        genPhiP4_my.clear();
        genPxP4_my.clear();
        genPyP4_my.clear();
        genPzP4_my.clear();
        genPP4_my.clear();
        genPreRap_my.clear();
        genPreEn_my.clear();
        genPt_my.clear();
        genEta_my.clear();
        genPhi_my.clear();
        genCharge_my.clear();
        genPx_my.clear();
        genPy_my.clear();
        genPz_my.clear();
        genP_my.clear();

  /*|-------------------------------------------------------|
    |                Reco Muon Part                           |
    |-------------------------------------------------------|*/
     Mu_isTightMuon.clear();
     Mu_isMediumMuon.clear();
     Mu_isGlobalMuon.clear();
     Mu_isHighPtMuon.clear();
     Mu_isTrackerMuon.clear();
     Mu_isLooseMuon.clear();
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
     emu_commonVrt.clear();

  /*|-------------------------------------------------------|
    |                PileUp part                            |
    |-------------------------------------------------------|*/
	nPU_ch   = -999;
	puTrue_ch = -999;
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

        elPtP4.clear();
        elEtaP4.clear();
        elPhiP4.clear();
        elPxP4.clear();
        elPyP4.clear();
        elPzP4.clear();
        elPP4.clear();
        elEnSC.clear();
        elEtaSC.clear();
        elPhiSC.clear();
        elisMedium.clear();
        isMediumOwn_El.clear();
        elEt.clear();
        elPt.clear();
        elEta.clear();
        elPhi.clear();
        elPx.clear();
        elPy.clear();
        elPz.clear();
        elP.clear();
        elCharg.clear();
        
        El_init_vtxNormChi.clear();
        El_init_vrtMass.clear();

        isoDeltaBeta_.clear();

        r9_el.clear();
        isoPFUnc_el.clear();
        isoPFRho_el.clear();
        sigmaIEtaIEta_el.clear();
        hadronicOverEm_el.clear();
        hasMatchedConversion_el.clear();
        MissedInnerTrackHits_el.clear();
        EleInnerTrackHits_el.clear();
        EleLostInnerTrackHits_el.clear();
        gsfTrackPt_el.clear();
        SCRawE_el.clear();
        eSCOverP_el.clear(); 
        dPhiSCTrackAtVtx_el.clear();
        dEtaSCTrackAtVtx_el.clear(); 
        epDiff_el.clear(); 
        d0_el.clear();
        dz_el.clear();
        isEB_el.clear();
        isEE_el.clear();

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
DEFINE_FWK_MODULE(TauEl_fin);
bool TauEl_fin::isPassHLT(){ 
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

bool TauEl_fin::RecoHLTMuonMatching(float RecoEta,float RecoPhi){
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
float TauEl_fin::delR(float eta1,float phi1,float eta2,float phi2){
	float mpi=3.14;
	float dp=std::abs(phi1-phi2);
	if (dp>mpi) dp-=float(2*mpi);
	return sqrt((eta1-eta2)*(eta1-eta2) + dp*dp);
	}
float TauEl_fin::TMass_F(float pt3lep, float px3lep, float py3lep, float met, float metPhi) {
    return sqrt(pow(pt3lep + met, 2) - pow(px3lep + met * cos(metPhi), 2) - pow(py3lep + met * sin(metPhi), 2));
}
double TauEl_fin::getEffectiveAreaForElectrons(const double& eta) const 
{
  double abseta = abs(eta);
    // values for 0.3 cone
    // values for 25ns runs
  if(abseta <= 1.)                            return 0.1752;
  else if(abseta > 1.    && abseta <= 1.479)  return 0.1862;
  else if(abseta > 1.479 && abseta <= 2.0)    return 0.1411;
  else if(abseta > 2.    && abseta <= 2.2)    return 0.1534;
  else if(abseta > 2.2   && abseta <=  2.3)   return 0.1903;
  else if(abseta > 2.3   && abseta <= 2.4)    return 0.2243;
  else if(abseta > 2.4   && abseta <= 2.5)    return 0.2687;
  else                                        return 9999;
}

const reco::Candidate* mother(const reco::CandidateBaseRef& cand) {
  int pId = cand->pdgId();
  const reco::Candidate* mom = cand->mother();
  while (mom != 0 && mom->pdgId() == pId)
    mom = mom->mother();
  return mom;
}
  
int motherId(const reco::CandidateBaseRef& cand) {
  const reco::Candidate* mom = mother(cand);
  if (mom != 0) return mom->pdgId();
  else return 0;
  /*    throw cms::Exception("motherId")
        << "+++ cannot find mother for particle! pdgId: " << pId
        << " status: " << cand->status() << " +++\n";*/
}

int grandmotherId(const reco::CandidateBaseRef& cand) {
  const reco::Candidate *mom = mother(cand);
  if (mom == 0) return 0;
  const reco::Candidate *gm = mom->mother();
  if (gm != 0) return gm->pdgId();
  else return 0;
  /* throw cms::Exception("grandmotherId")
     << "+++ cannot find grandmother for particle! pdgId: " << p->pdgId()
     << " status: " << cand->status() 
     << " mother pdgId: " << mId << " status: " << mom->status() << " +++\n";*/
}

