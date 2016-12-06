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
	'/store/data/Run2015D/SingleMuon/MINIAOD/PromptReco-v4/000/258/159/00000/700FDBE1-246C-E511-ABD5-02163E0143C4.root',
        #'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/168/00000/60FF8405-EA26-E511-A892-02163E01387D.root'
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

process.GlobalTag.globaltag = cms.string('74X_dataRun2_Prompt_v4')

process.demo = cms.EDAnalyzer('Moge_highPtID_Data',
    vertices            = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot            = cms.InputTag("offlineBeamSpot"),
    muons               = cms.InputTag("slimmedMuons"),
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
    #packedGenParticles  = cms.InputTag("packedGenParticles"),
    #prunedGenParticles  = cms.InputTag("prunedGenParticles"),
    pileupCollection    = cms.InputTag("slimmedAddPileupInfo"),
)
#process.Tracer = cms.Service("Tracer")
process.TFileService = cms.Service("TFileService",
          fileName = cms.string('Run2015D_PromptReco-v4.root')
)

process.p = cms.Path(process.demo)