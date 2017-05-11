# DON'T FORGET TO CHANGE FOR MC AND DATA: 
# 1. GLOBAL TAG 
# 2. isMC flag
# 3. filename   

import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),
SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
	#'file:/home/uyevarou/Analiz/CMSSW_7_4_15/src/Kokosh/forMiniAod_an/VersionMC_another.root',
        #'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root',
        #'file:/afs/cern.ch/user/u/uyevarou/Nova/CMSSW_7_4_15/src/VersionMC_another.root',
	#'/store/data/Run2015D/SingleMuon/MINIAOD/PromptReco-v4/000/258/159/00000/700FDBE1-246C-E511-ABD5-02163E0143C4.root',
        #'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/00759690-D16E-E511-B29E-00261894382D.root',
        #'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/0EE9C450-7A6E-E511-854D-0025905A48F2.root',
        #'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/180D195E-766F-E511-9AF8-00266CF89498.root',
        '/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/02BB716A-A36D-E511-9280-002590743042.root',
        #'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/5C30C753-636D-E511-9FA2-008CFA197998.root',
        #'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v3/50000/323F5423-A26F-E511-A847-02163E00EA7B.root',
        #'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/7EC09472-D96C-E511-A113-38EAA7A30576.root',
        #'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/30F44BB2-D96C-E511-835F-008CFA0A5844.root',
        #'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_74X_mcRun2_asymptotic_v2-v1/00000/C6004B8D-F8B7-E511-BE25-0026B9277A25.root',
        #'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/1C4F406D-D96C-E511-A105-B083FED42A6E.root',
        #'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/76B75B02-F36C-E511-8384-02163E01698B.root',
         )
)

from Configuration.StandardSequences.GeometryRecoDB_cff import *
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.load("RecoTracker.Configuration.RecoTracker_cff")
process.load('Configuration/StandardSequences/EndOfProcess_cff')

process.GlobalTag.globaltag = cms.string('74X_mcRun2_asymptotic_v2')
#process.GlobalTag.globaltag = cms.string('74X_dataRun2_Prompt_v4')


from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.demo = cms.EDAnalyzer('TauEl_fin',
    isMC                = cms.untracked.bool(True), 
    vertices            = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot            = cms.InputTag("offlineBeamSpot"),
    muons               = cms.InputTag("slimmedMuons"),
    electrons           = cms.InputTag("slimmedElectrons"),
    gsfelectrons        = cms.InputTag("reducedEGamma"),
    jets                = cms.InputTag("slimmedJets"),
    fatjets             = cms.InputTag("slimmedJetsAK8"),
    taus                = cms.InputTag("slimmedTaus"),
    photons             = cms.InputTag("slimmedPhotons"),
    genJets             = cms.InputTag("slimmedGenJets"),
    mets                = cms.InputTag("slimmedMETs"),
    bits                = cms.InputTag("TriggerResults","","HLT"),
    conversions         = cms.InputTag("reducedEgamma:reducedConversions"),
    prescales          = cms.InputTag("patTrigger"),
    objects             = cms.InputTag("selectedPatTrigger"),  
    packedPFCandidates  = cms.InputTag("packedPFCandidates"),
    srcRho              = cms.InputTag('fixedGridRhoFastjetAll'),
    triggerMu           = cms.untracked.string("HLT_IsoMu20"),
    triggerTrkMu        = cms.untracked.string("HLT_IsoTkMu20"),
    packedGenParticles  = cms.InputTag("packedGenParticles"),
    prunedGenParticles  = cms.InputTag("prunedGenParticles"),
    pileupCollection    = cms.InputTag("slimmedAddPileupInfo"),
    lheEvtInfo          = cms.InputTag("externalLHEProducer"),
    externalLHEProducer = cms.InputTag("externalLHEProducer"),
    #
    # ID decisions (common to all formats)
    #
    # for 2015 recommended CMSSW_7_5_16+
    # https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns
    #eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    #eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    #eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    #eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
    # for 74X
    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
    eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-HEEPIDV51"),
)
process.primaryVertexFilter  = cms.EDFilter("VertexSelector",
      src = cms.InputTag('offlineSlimmedPrimaryVertices'),
      cut = cms.string('!isFake && ndof > 4.0 && position.Rho < 2.0 && abs(z) < 24'),
      filter = cms.bool(True)  ## otherwise it won't filter the events, just produce an empty vertex collection.
)
#process.Tracer = cms.Service("Tracer")
process.TFileService = cms.Service("TFileService",
          #fileName = cms.string('/afs/cern.ch/work/u/uyevarou/private/DY200/DY200in.root'),
          fileName = cms.string('MC_fin.root'),
          #fileName = cms.string('Data_mu_tauEl.root')
)

process.p = cms.Path(process.egmGsfElectronIDSequence * process.primaryVertexFilter * process.demo)
