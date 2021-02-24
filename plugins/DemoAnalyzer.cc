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
//using reco::TrackCollection;

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
//	edm::EDGetTokenT<TrackCollection> 			tracksToken_;  //used to select what tracks to read from configuration file
	TH1D *demohisto;
	edm::Service<TFileService> fs;
	edm::EDGetTokenT<GenEventInfoProduct>			genInfoToken_;
       	edm::EDGetTokenT<std::vector<reco::GenJet> >		genjetToken_;
     	edm::EDGetTokenT<std::vector<reco::GenParticle>>	genparticleToken_;
      	edm::EDGetTokenT<edm::View<PileupSummaryInfo> >		pileupToken_;
	edm::EDGetTokenT<std::vector<SimVertex> >		verticesToken_;
	edm::EDGetTokenT<std::vector<SimTrack> >		tracksToken_;







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
   //tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))
 // genInfoToken_(consumes<GenEventInfoProduct>(iConfig.getUntrackedParameter<edm::InputTag>("generatorInfo")))


//   minTracks_(iConfig.getUntrackedParameter<unsigned int>("minTracks",0))
{
   verbose_		= iConfig.getUntrackedParameter<bool>("verbose");
   sampleID_		= iConfig.getUntrackedParameter<int>("sampleID",0);
//   genInfoToken_	= consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"));
   genjetToken_		= consumes<std::vector<reco::GenJet> >(iConfig.getUntrackedParameter<edm::InputTag>("genjets"));
   genparticleToken_    = consumes<std::vector<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("genparticles"));
   pileupToken_		= consumes<edm::View<PileupSummaryInfo> >(iConfig.getUntrackedParameter<edm::InputTag>("pileupInfo"));
   verticesToken_	= consumes<std::vector<SimVertex> >(iConfig.getUntrackedParameter<edm::InputTag>("vertices"));
    tracksToken_		= consumes<std::vector<SimTrack> >(iConfig.getUntrackedParameter<edm::InputTag>("tracks"));

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
   if (verbose_) std::cout << "Done" << std::endl;

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
 
  
  edm::Handle<std::vector<reco::GenJet> > genjets;
  iEvent.getByToken(genjetToken_, genjets);

  edm::Handle<std::vector<reco::GenParticle> > genparticles;
  iEvent.getByToken(genparticleToken_, genparticles);

  edm::Handle<edm::View<PileupSummaryInfo> > pileupInfo;
  iEvent.getByToken(pileupToken_, pileupInfo);

  edm::Handle<std::vector<SimVertex> > vertices;
  iEvent.getByToken(verticesToken_, vertices);

  edm::Handle<std::vector<SimTrack> > tracks;
  iEvent.getByToken(tracksToken_, tracks);

  unsigned long int event = iEvent.id().event();   
  int run                 = iEvent.id().run();

  //int ngenjets= genjets->size();
  // ngenjets=0;
  std::vector<float>	genjet_pt;
  std::vector<float>	genjet_e;
  std::vector<float>	genjet_eta;
  std::vector<float>	genjet_phi;
  std::vector<int>  	genjet_ndaug;
  std::vector<int>  	genjet_nconst;
  std::vector<float>	genjet_vx;
  std::vector<float>	genjet_vy;
  std::vector<float>	genjet_vz;

  // --- genjet constituent info
  int genjet0_nconst = 0;
  int genjet1_nconst = 0;
  int genjet2_nconst = 0;
  int genjet3_nconst = 0;
  std::vector<int>  	genjet_i;
  std::vector<int>  	genjet_match;
  std::vector<float>	genjet0_const_st;
  std::vector<float>	genjet0_const_id;
  std::vector<float>	genjet0_const_pt;
  std::vector<float>	genjet0_const_pv;
  std::vector<float>	genjet0_const_theta;
  std::vector<float>	genjet0_const_mtheta;
  std::vector<float>	genjet0_const_vx;
  std::vector<float>	genjet0_const_vy;
  std::vector<float>	genjet0_const_vz;
  std::vector<float>	genjet1_const_st;
  std::vector<float>	genjet1_const_id;
  std::vector<float>	genjet1_const_pt;
  std::vector<float>	genjet1_const_pv;
  std::vector<float>	genjet1_const_theta;
  std::vector<float>	genjet1_const_mtheta;
  std::vector<float>	genjet1_const_vx;
  std::vector<float>	genjet1_const_vy;
  std::vector<float>	genjet1_const_vz;
  std::vector<float>	genjet2_const_st;
  std::vector<float>	genjet2_const_id;
  std::vector<float>	genjet2_const_pt;
  std::vector<float>	genjet2_const_pv;
  std::vector<float>	genjet2_const_theta;
  std::vector<float>	genjet2_const_mtheta;
  std::vector<float>	genjet2_const_vx;
  std::vector<float>	genjet2_const_vy;
  std::vector<float>	genjet2_const_vz;
  std::vector<float>	genjet3_const_st;
  std::vector<float>	genjet3_const_id;
  std::vector<float>	genjet3_const_pt;
  std::vector<float>	genjet3_const_pv;
  std::vector<float>	genjet3_const_theta;
  std::vector<float>	genjet3_const_mtheta;
  std::vector<float>	genjet3_const_vx;
  std::vector<float>	genjet3_const_vy;
  std::vector<float>	genjet3_const_vz;

  //if (genjets.isValid()){ // make sure have genjet collection
    //ngenjets = genjets->size();
  std::cout<<"test"<<std::endl;  
    int jetiter = 0;
    for (const auto & genjet_iter : *genjets){ // loop over genjets
 
      // standard information
 	genjet_pt.push_back(genjet_iter.pt()); 
	genjet_e.push_back(genjet_iter.energy()); 
	genjet_eta.push_back(genjet_iter.eta()); 
	genjet_phi.push_back(genjet_iter.phi());
	genjet_ndaug.push_back(genjet_iter.numberOfDaughters());
	genjet_vx.push_back(genjet_iter.vertex().x());
	genjet_vy.push_back(genjet_iter.vertex().y());
	genjet_vz.push_back(genjet_iter.vertex().z());
	genjet_i.push_back(jetiter); // note: probably not necessary
                                                              
	// constituent information
	std::vector<const reco::GenParticle*> genjet_const = genjet_iter.getGenConstituents();
	genjet_nconst.push_back(genjet_const.size()); 
 // }
}
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
