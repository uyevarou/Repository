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
	'file:/afs/cern.ch/user/u/uyevarou/Nova/CMSSW_7_4_15/src/VersionMC_another.root',
        #'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root',
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
#process.GlobalTag.globaltag = cms.string('MCRUN2_74_V9')

process.demo = cms.EDAnalyzer('Moge_tightID',
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
    packedGenParticles  = cms.InputTag("packedGenParticles"),
    prunedGenParticles  = cms.InputTag("prunedGenParticles"),
    pileupCollection    = cms.InputTag("slimmedAddPileupInfo"),
)
#process.Tracer = cms.Service("Tracer")
process.TFileService = cms.Service("TFileService",
          fileName = cms.string('RecoGenMu_tightID_c.root')
)

process.p = cms.Path(process.demo)