import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing


process = cms.Process("Demo")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

opts = VarParsing.VarParsing('analysis')
opts.register('debug',                                      # option name
              False,                                        # default value
              VarParsing.VarParsing.multiplicity.singleton, # singleton or list
              VarParsing.VarParsing.varType.bool,           # type: string, int, float, bool
              "Print out debug info")                       # description
opts.register('sampleID',
              100,                                          # default value
              VarParsing.VarParsing.multiplicity.singleton, # singleton or list
              VarParsing.VarParsing.varType.int,            # type: string, int, float, bool
              "SampleID [1-99] Bkg, [100+] Sig")            # description
opts.parseArguments()

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
	SkipEvent = cms.untracked.vstring('ProductNotFound')
	 )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
	    'root://cmsxrootd.fnal.gov////store/mc/RunIISummer16MiniAODv2/ggXToYYTo2mu2e_m14_PseudoScalar_13TeV-pythia8-JHUGen/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/10C23D4F-94BD-E811-9588-E0071B7B2320.root'	
            #'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
                )
                            )



process.TFileService = cms.Service("TFileService",
      fileName = cms.string("histo.root"),
      closeFileFast = cms.untracked.bool(True)
  )
#process.demo = cms.EDAnalyzer('DemoAnalyzer'
#                              )

process.demo = cms.EDAnalyzer('DemoAnalyzer',
	                        verbose         = cms.untracked.bool(opts.debug),
				sampleID	= cms.untracked.int32(opts.sampleID),
				generatorInfo	= cms.InputTag("generator"),
				genjets		= cms.untracked.InputTag("ak4GenJets", "", "SIM"),
                                genparticles    = cms.untracked.InputTag("genParticles", "", "SIM"), #era SIM
				pileupInfo	= cms.untracked.InputTag("slimmedAddPileupInfo"),
				vertices	= cms.untracked.InputTag("g4SimHits", "", "SIM"),
          			tracks		= cms.untracked.InputTag("g4SimHits", "", "SIM"),
         )

process.p = cms.Path(process.demo)
