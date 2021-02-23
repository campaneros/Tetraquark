// -*- C++ -*-
//
// Package:    demo/DemoAnalyzer
// Class:      DemoAnalyzer
//
/**\class DemoAnalyzer DemoAnalyzer.cc demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mattia Campana
//         Created:  Tue, 09 Feb 2021 13:55:35 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"


#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TTree.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


///
//
using reco::TrackCollection;

struct tree_struc_{

  int				sample;
  int				run;
  int				lumi;
  long unsigned int		event;
  float				weight;
  int				ngenjets;
};

class DemoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


	bool	verbose_;
	int	sampleID_;
	      //unsigned int minTracks_;
      // ----------member data ---------------------------
	edm::EDGetTokenT<TrackCollection> 			tracksToken_;  //used to select what tracks to read from configuration file
	TH1D *demohisto;
	edm::Service<TFileService> fs;
	edm::EDGetTokenT<GenEventInfoProduct>			genInfoToken_;
       	edm::EDGetTokenT<std::vector<reco::GenJet> >		genjetToken_;
     	edm::EDGetTokenT<std::vector<reco::GenParticle>>	genparticleToken_;
    	edm::EDGetTokenT<TriggerResults>                        triggerResults_;
      	edm::EDGetTokenT<edm::View<PileupSummaryInfo> >		pileupToken_;
	edm::EDGetTokenT<std::vector<SimVertex> >		verticesToken_;








//////////////////////
//setup Tree/////////
////////////////////
      TTree* tree;
      tree_struc_ tree_;

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
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig)
 :
  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))
 // genInfoToken_(consumes<GenEventInfoProduct>(iConfig.getUntrackedParameter<edm::InputTag>("generatorInfo")))


//   minTracks_(iConfig.getUntrackedParameter<unsigned int>("minTracks",0))
{
   verbose_		= iConfig.getUntrackedParameter<bool>("verbose");
   sampleID_		= iConfig.getUntrackedParameter<int>("sampleID",0);
   genInfoToken_	= consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generatorInfo"));
   genjetToken_		= consumes<std::vector<reco::GenJet> >(iConfig.getUntrackedParameter<edm::InputTag>("genjets"));
   genparticleToken_    = consumes<std::vector<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("genparticles"));
   pileupToken_		= consumes<edm::View<PileupSummaryInfo> >(iConfig.getUntrackedParameter<edm::InputTag>("pileupInfo"));
   verticesToken_	= consumes<std::vector<SimVertex> >(iConfig.getUntrackedParameter<edm::InputTag>("vertices"));
 

 // verbose_		= iConfig.getUntrackedParameter<bool>("verbose");
  // sampleID_		= iConfig.getUntrackedParameter<int>("sampleID",0);
  ///edm::Service<TFileService> fs;
  // demohisto = fs->make<TH1D>("tracks" , "Tracks" , 100 , 0 , 5000 );
//   genInfoToken_ = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generatorInfo"));
   //now do what ever initialization is needed
 //  genInfoToken_	= consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generatorInfo"));
  //genparticleToken_    = consumes<std::vector<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("genparticles"));
 // triggerResults_ = consumes<TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("triggerreuslts")) 
  // pileupToken_		= consumes<edm::View<PileupSummaryInfo> >(iConfig.getUntrackedParameter<edm::InputTag>("pileupInfo"));
  // verticesToken_	= consumes<std::vector<SimVertex> >(iConfig.getUntrackedParameter<edm::InputTag>("vertices"));
   //tracksToken_		= consumes<std::vector<SimTrack> >(iConfig.getUntrackedParameter<edm::InputTag>("tracks"));


}


DemoAnalyzer::~DemoAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  /*edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(genInfoToken_, genInfo);  
 
  //edm::Handle<std::vector<reco::GenJet> > genjets;
  //iEvent.getByToken(genjetToken_, genjets);

  edm::Handle<std::vector<reco::GenParticle> > genparticles;
  iEvent.getByToken(genparticleToken_, genparticles);

  edm::Handle<edm::View<PileupSummaryInfo> > pileupInfo;
  iEvent.getByToken(pileupToken_, pileupInfo);

  edm::Handle<std::vector<SimVertex> > vertices;
  iEvent.getByToken(verticesToken_, vertices);

  edm::Handle<std::vector<SimTrack> > tracks;
  iEvent.getByToken(tracksToken_, tracks);*/
    Handle<TrackCollection> tracks;
    iEvent.getByToken(tracksToken_, tracks);
   // iEvent.getByToken("generalTracks", tracks);  
 //  if( minTracks_ <= tracks->size() ) {
   //     LogInfo ("Demo") << "number of tracks " << tracks->size();
     //}
// for(TrackCollection::const_iterator itTrack = tracks->begin();
   //     itTrack != tracks->end();
     //   ++itTrack) {
      // do something with track parameters, e.g, plot the charge.
      // int charge = itTrack->charge();
    //}

 //  Handle<std::vector<reco::GenParticle>> genParticles;
 //  iEvent.getByToken(genparticleToken_, genParticles);

  /* using namespace reco;

  for(size_t i = 0; i < genParticles->size(); ++ i) {
     const GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
     int st = p.status();  
     const Candidate * mom = p.mother();
     double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
     double vx = p.vx(), vy = p.vy(), vz = p.vz();
     int charge = p.charge();
     int n = p.numberOfDaughters();
     for(int j = 0; j < n; ++ j) {
       const Candidate * d = p.daughter( j );
       int dauId = d->pdgId();
       
	demohisto->Fill(dauId);
	//std::cout<< dauId  <<std::endl;
	// . . . 
           }
       //           // . . . 
        }
*/
 #ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
DemoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DemoAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
