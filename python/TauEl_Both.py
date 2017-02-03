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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
	#'file:/home/uyevarou/Analiz/CMSSW_7_4_15/src/Kokosh/forMiniAod_an/VersionMC_another.root',
        #'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root',
        'file:/afs/cern.ch/user/u/uyevarou/Nova/CMSSW_7_4_15/src/VersionMC_another.root',
	#'/store/data/Run2015D/SingleMuon/MINIAOD/PromptReco-v4/000/258/159/00000/700FDBE1-246C-E511-ABD5-02163E0143C4.root',
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

process.demo = cms.EDAnalyzer('TauEl_Both_b',
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
          fileName = cms.string('MC_mu_tauEl.root'),
          #fileName = cms.string('Data_mu_tauEl.root')
)

process.p = cms.Path(process.egmGsfElectronIDSequence * process.primaryVertexFilter * process.demo)
