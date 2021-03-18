import FWCore.ParameterSet.Config as cms

process = cms.Process("GenAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/afs/cern.ch/user/s/soffi/public/10C23D4F-94BD-E811-9588-E0071B7B2320.root'
                )
                            )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple_tetraquarks.root")
)


process.GenAnalysis = cms.EDAnalyzer('TQGenAnalyzer',
                                     generatorInfo= cms.InputTag("generator"),
                                     genparticles    = cms.untracked.InputTag("prunedGenParticles", "", "PAT"),
                                     muons = cms.untracked.InputTag("slimmedMuons","","PAT"),
                                     electrons = cms.untracked.InputTag("slimmedElectrons","","PAT"),
                              )

process.p = cms.Path(process.GenAnalysis)
