
// -*- C++ -*-
//
// Package:    GenAnalysis/TQGenAnalyzer
// Class:      TQGenAnalyzer
//
/**\class TQGenAnalyzer TQGenAnalyzer.cc GenAnalysis/TQGenAnalyzer/plugins/TQGenAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Livia Soffi
//         Created:  Mon, 01 Mar 2021 11:15:18 GMT
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
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SimDataFormats/Associations/interface/MuonToTrackingParticleAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include <SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h>

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEgamma/ElectronIdentification/interface/ElectronIDAlgo.h"


#include "DataFormats/Math/interface/deltaR.h"
// root include files
#include "TTree.h"
#include "TLorentzVector.h"

struct tree_struc_{
  int run;
  int lumi;
  long unsigned int event;
  int nEle;
  int nMu;
  int nEleReco;
  int nMuReco;
  float TQ_genMass;
  float TQ_recoMass;
  std::vector<float>            genLep_pt;
  std::vector<float>            genLep_eta;
  std::vector<float>            genLep_mass;
  std::vector<float>            genLep_phi;
  std::vector<int>            genLep_pdgId;
  std::vector<float>            genMom_pt;
  std::vector<float>            genMom_eta;
  std::vector<float>            genMom_mass;
  std::vector<float>            genMom_phi;
  std::vector<int>            genMom_pdgId;

  std::vector<float>            recoMu_pt;
  std::vector<float>            recoMu_eta;
  std::vector<float>            recoMu_mass;
  std::vector<float>            recoMu_phi;
  std::vector<float>            recoMu_dR;
  std::vector<int>              recoMu_charge;


  std::vector<float>            recoEle_pt;
  std::vector<float>            recoEle_eta;
  std::vector<float>            recoEle_mass;
  std::vector<float>            recoEle_phi;
  std::vector<float>            recoEle_dR;
  std::vector<int>              recoEle_charge;


  std::vector<float>            recoEle_pt_un;
  std::vector<float>            recoMu_pt_un;

};


class TQGenAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TQGenAnalyzer(const edm::ParameterSet&);
      ~TQGenAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  void initTreeStructure();
  void clearVectors(); 

      // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  edm::EDGetTokenT<GenEventInfoProduct>genInfoToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> >genparticleToken_;
  edm::EDGetTokenT<std::vector<pat::Muon> > muonToken_;
  edm::EDGetTokenT<std::vector<pat::Electron> > eleToken_;


  // setup tree;
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
TQGenAnalyzer::TQGenAnalyzer(const edm::ParameterSet& iConfig){
   //now do what ever initialization is needed
  genInfoToken_= consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generatorInfo"));
  genparticleToken_    = consumes<std::vector<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("genparticles"));
  muonToken_=consumes<std::vector<pat::Muon>>(iConfig.getUntrackedParameter<edm::InputTag>("muons"));
  eleToken_=consumes<std::vector<pat::Electron>>(iConfig.getUntrackedParameter<edm::InputTag>("electrons"));
}


TQGenAnalyzer::~TQGenAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TQGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // --- pick up handles 
  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(genInfoToken_, genInfo);  

  edm::Handle<std::vector<reco::GenParticle> > genparticles;
  iEvent.getByToken(genparticleToken_, genparticles);

  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByToken(eleToken_, electrons);


  std::vector<float>genLep_pt;
  std::vector<float>genLep_mass;
  std::vector<float>genLep_eta;
  std::vector<float>genLep_phi;
  std::vector<int>genLep_pdgId;

  std::vector<float>genMom_pt;
  std::vector<float>genMom_mass;
  std::vector<float>genMom_eta;
  std::vector<float>genMom_phi;
  std::vector<int>genMom_pdgId;

  std::vector<float>recoMu_pt;
  std::vector<float>recoMu_mass;
  std::vector<float>recoMu_eta;
  std::vector<float>recoMu_phi;
  std::vector<float>recoMu_dR;
  std::vector<int>recoMu_charge;


  std::vector<float>recoEle_pt;
  std::vector<float>recoEle_mass;
  std::vector<float>recoEle_eta;
  std::vector<float>recoEle_phi;
  std::vector<float>recoEle_dR;
  std::vector<int>recoEle_charge;


  std::vector<float>recoEle_pt_un;
  std::vector<float>recoMu_pt_un;

 // --- general event info 
  unsigned long int event = iEvent.id().event();   
  int run                 = iEvent.id().run();
  int lumi                = iEvent.luminosityBlock();
  int nEle=0;
  int nMu=0;
  int nMuReco=0;
  int nEleReco=0;
  int first_lepton=0;

  TLorentzVector* lep1=new TLorentzVector();
  TLorentzVector* lep2=new TLorentzVector();
  TLorentzVector* lep3=new TLorentzVector();
  TLorentzVector* lep0=new TLorentzVector();
  

  TLorentzVector* recolep1=new TLorentzVector();
  TLorentzVector* recolep2=new TLorentzVector();
  TLorentzVector* recolep3=new TLorentzVector();
  TLorentzVector* recolep0=new TLorentzVector();
  

  std::cout<< "---------------------------------------------------------"<<std::endl;

  for (const auto & genpar_iter : *genparticles){ // loop over genparticles
      if(abs(genpar_iter.pdgId())!=13 && abs(genpar_iter.pdgId())!=11)continue;

      if(genpar_iter.mother(0)->pdgId()<20 || genpar_iter.mother(0)->pdgId()>200) continue;
      
      //      std::cout<<"status: "<<genpar_iter.status()<<" pdgid: "<<genpar_iter.pdgId()<<" pt: "<<genpar_iter.pt()<<" mass: "<<genpar_iter.mass()<<" eta: "<<genpar_iter.eta()<<" phi: "<<genpar_iter.phi()<<std::endl;

      
      genLep_pt.push_back(genpar_iter.pt());
      genLep_mass.push_back(genpar_iter.mass());
      genLep_eta.push_back(genpar_iter.eta());
      genLep_phi.push_back(genpar_iter.phi());
      genLep_pdgId.push_back(genpar_iter.pdgId());
      

      genMom_pt.push_back(genpar_iter.mother(0)->pt());
      genMom_mass.push_back(genpar_iter.mother(0)->mass());
      genMom_eta.push_back(genpar_iter.mother(0)->eta());
      genMom_phi.push_back(genpar_iter.mother(0)->phi());
      genMom_pdgId.push_back(genpar_iter.mother(0)->pdgId());
      
      if(abs(genpar_iter.pdgId())==11)nEle++;
      if(abs(genpar_iter.pdgId())==13)nMu++;
      
    }

   // now look at muons and find the matching reco ones
  for(unsigned int i=0;i<genLep_pdgId.size();i++){
    if(abs(genLep_pdgId[i])!=13)continue;
       	float dR=999.;
	float dRMin=999.;
	float recopt=999.;
	float recoeta=999.;
	float recophi=999.;
	float recomass=999.;
        int recocharge=999;
        float dRApp=999.;
        bool last = true;


    for (const auto & mu_iter : *muons){
      float eta=mu_iter.eta();
      float phi=mu_iter.phi();
    
      
       dR= deltaR(eta,phi,genLep_eta[i],genLep_phi[i]); 
	if(dR<dRMin){
	 recopt=mu_iter.pt();
	 recoeta=mu_iter.eta();
	 recophi=mu_iter.phi();
	 recomass=mu_iter.mass();
	 dRMin=dR;
         recocharge=mu_iter.charge();
       }
   
     }
  
      if (first_lepton!=0 && recopt!=999 && ((recopt==recoMu_pt.back() && recophi==recoMu_phi.back() && recoeta==recoMu_eta.back()))){
       if (dR>recoMu_dR.back()) dRApp=dR;
       else {
         // std::cout<<"second last element  "<<recoMu_pt.back()<<std::endl;
        dRApp=recoMu_dR.back();
        recoMu_pt.push_back(recopt);
        recoMu_mass.push_back(recomass);
        recoMu_eta.push_back(recoeta);
        recoMu_phi.push_back(recophi);
        recoMu_dR.push_back(dRMin);
        recoMu_charge.push_back(recocharge);
        last=false;
       }
       std::cout<<"found"<<std::endl;

       for (const auto & mu_iter : *muons){
        float eta=mu_iter.eta();
         float phi=mu_iter.phi();


        if (last) dR= deltaR(eta,phi,genLep_eta[i],genLep_phi[i]);
        else dR=deltaR(eta,phi,genLep_eta[i-1],genLep_phi[i-1]);
          if(dR<dRMin && dR!=dRApp){
             recopt=mu_iter.pt();
             recoeta=mu_iter.eta();
             recophi=mu_iter.phi();
             recomass=mu_iter.mass();
             dRMin=dR;
             recocharge=mu_iter.charge();
          }
        }
     if (last){
      recoMu_pt.push_back(recopt);
      recoMu_mass.push_back(recomass);
      recoMu_eta.push_back(recoeta);
      recoMu_phi.push_back(recophi);
      recoMu_dR.push_back(dRMin);
      recoMu_charge.push_back(recocharge);
     }else{
    //std::cout<<"second last element  "<<recoMu_pt.end()[-2]<<std::endl;
    recoMu_pt.end()[-2]=recopt;
    recoMu_mass.end()[-2]=recomass;
    recoMu_phi.end()[-2]=recophi;
    recoMu_eta.end()[-2]=recoeta;
    recoMu_dR.end()[-2]=dRMin;
    recoMu_charge.end()[-2]=recocharge;
   }
   } else {
    recoMu_pt.push_back(recopt);
    recoMu_mass.push_back(recomass);
    recoMu_eta.push_back(recoeta);
    recoMu_phi.push_back(recophi);
    recoMu_dR.push_back(dRMin);
    recoMu_charge.push_back(recocharge);
       }
     
    first_lepton++;
     
   }
 
  first_lepton=0;

	
   // now look at electrons and find the matching reco ones
  for(unsigned int i=0;i<genLep_pdgId.size();i++){
   if(abs(genLep_pdgId[i])!=11)continue;
    float dR=999.;
    float dRMin=999.;
    float recopt=999.;
    float recoeta=999.;
    float recophi=999.;
    float recomass=999.;
    int recocharge=999;    
    float dRApp=999.;
    bool last = true; 
     
 
   for (const auto & el_iter : *electrons){ 
      float eta=el_iter.eta();
      float phi=el_iter.phi();

  
     
    dR= deltaR(eta,phi,genLep_eta[i],genLep_phi[i]);
      if(dR<dRMin){
	recopt=el_iter.pt();
	recoeta=el_iter.eta();
	recophi=el_iter.phi();
	recomass=el_iter.mass();
	dRMin=dR;
	recocharge=el_iter.charge();
     }

    }
   
   // std::cout<<"inx "<<first_lepton<<std::endl;
   // std::cout<<"index"<<i<<std::endl;
      if (first_lepton!=0 && recopt!=999 && ((recopt==recoEle_pt.back() && recophi==recoEle_phi.back() && recoeta==recoEle_eta.back()))){
       if (dR>recoEle_dR.back()) dRApp=dR;
       else {

    //      std::cout<<"second last element  "<<recoEle_pt.back()<<std::endl;
        dRApp=recoEle_dR.back();
        recoEle_pt.push_back(recopt);
        recoEle_mass.push_back(recomass);
        recoEle_eta.push_back(recoeta);
        recoEle_phi.push_back(recophi);
        recoEle_dR.push_back(dRMin);
        recoEle_charge.push_back(recocharge);
        last=false;	
       }
    
       std::cout<<"found"<<std::endl;
    

       for (const auto & el_iter : *electrons){
        float eta=el_iter.eta();
        float phi=el_iter.phi();



        if (last) dR= deltaR(eta,phi,genLep_eta[i],genLep_phi[i]);
        else dR=deltaR(eta,phi,genLep_eta[i-1],genLep_phi[i-1]);
          if(dR<dRMin && dR!=dRApp){
             recopt=el_iter.pt();
             recoeta=el_iter.eta();
             recophi=el_iter.phi();
             recomass=el_iter.mass();
             dRMin=dR;
             recocharge=el_iter.charge();
          }
        } 
     if (last){
      recoEle_pt.push_back(recopt);
      recoEle_mass.push_back(recomass);
      recoEle_eta.push_back(recoeta);
      recoEle_phi.push_back(recophi);
      recoEle_dR.push_back(dRMin);
      recoEle_charge.push_back(recocharge);
     }else{
    //std::cout<<"second last element  "<<recoEle_pt.end()[-2]<<std::endl;
    recoEle_pt.end()[-2]=recopt;
    recoEle_mass.end()[-2]=recomass;
    recoEle_phi.end()[-2]=recophi;
    recoEle_eta.end()[-2]=recoeta;    
    recoEle_dR.end()[-2]=dRMin;
    recoEle_charge.end()[-2]=recocharge;


   }

     
   } else {
   
      
     // dRApp=999.;
      recoEle_pt.push_back(recopt);
      recoEle_mass.push_back(recomass);
      recoEle_eta.push_back(recoeta);
      recoEle_phi.push_back(recophi);
      recoEle_dR.push_back(dRMin);
      recoEle_charge.push_back(recocharge);
  }

  first_lepton++;
  } 


  for (const auto & el_iter : *electrons){
  nEleReco++;

  float reco_ele=999.;
  reco_ele=el_iter.pt();
  recoEle_pt_un.push_back(reco_ele);
  }
  
  for (const auto & mu_iter : *muons){
  nMuReco++;

  float reco_mu=999.;
  reco_mu=mu_iter.pt();
  recoMu_pt_un.push_back(reco_mu);
  }

  // --- setup tree values
  initTreeStructure();

  // --- clear all vectors
  clearVectors();
  
    // --- fill trees
   tree_.run=run;
   tree_.event=event;
   tree_.lumi=lumi;
   tree_.nEle=nEle;
   tree_.nMu=nMu;
   tree_.nMuReco=nMuReco;
   tree_.nEleReco=nEleReco;  
   first_lepton=0;
 
   for(unsigned int i=0;i<genLep_pt.size();i++){
      tree_.genLep_pt.push_back(genLep_pt[i]);
      tree_.genLep_mass.push_back(genLep_mass[i]);
      tree_.genLep_eta.push_back(genLep_eta[i]);
      tree_.genLep_phi.push_back(genLep_phi[i]);
      tree_.genLep_pdgId.push_back(genLep_pdgId[i]);

      tree_.genMom_pt.push_back(genMom_pt[i]);
      tree_.genMom_mass.push_back(genMom_mass[i]);
      tree_.genMom_eta.push_back(genMom_eta[i]);
      tree_.genMom_phi.push_back(genMom_phi[i]);
      tree_.genMom_pdgId.push_back(genMom_pdgId[i]);
    

      if(i==0)lep0->SetPtEtaPhiM(genLep_pt[i],genLep_eta[i],genLep_phi[i],genLep_mass[i]);
      if(i==1)lep1->SetPtEtaPhiM(genLep_pt[i],genLep_eta[i],genLep_phi[i],genLep_mass[i]);
      if(i==2)lep2->SetPtEtaPhiM(genLep_pt[i],genLep_eta[i],genLep_phi[i],genLep_mass[i]);
      if(i==3)lep3->SetPtEtaPhiM(genLep_pt[i],genLep_eta[i],genLep_phi[i],genLep_mass[i]);
    }

    TLorentzVector TQ(*lep0+*lep1+*lep2+*lep3);
    tree_.TQ_genMass=TQ.M();

    for(unsigned int i=0;i<recoMu_pt.size();i++){
      tree_.recoMu_pt.push_back(recoMu_pt[i]);
      tree_.recoMu_mass.push_back(recoMu_mass[i]);
      tree_.recoMu_eta.push_back(recoMu_eta[i]);
      tree_.recoMu_phi.push_back(recoMu_phi[i]);   
      tree_.recoMu_dR.push_back(recoMu_dR[i]);   
      tree_.recoMu_charge.push_back(recoMu_charge[i]);   
 
      if(i==0)recolep0->SetPtEtaPhiM(recoMu_pt[i],recoMu_eta[i],recoMu_phi[i],recoMu_mass[i]);
      if(i==1)recolep1->SetPtEtaPhiM(recoMu_pt[i],recoMu_eta[i],recoMu_phi[i],recoMu_mass[i]);

}
    for(unsigned int i=0;i<recoEle_pt.size();i++){
      tree_.recoEle_pt.push_back(recoEle_pt[i]);
      tree_.recoEle_mass.push_back(recoEle_mass[i]);
      tree_.recoEle_eta.push_back(recoEle_eta[i]);
      tree_.recoEle_phi.push_back(recoEle_phi[i]);
      tree_.recoEle_dR.push_back(recoEle_dR[i]); 
      tree_.recoEle_charge.push_back(recoEle_charge[i]);
    

      if(i==0)recolep2->SetPtEtaPhiM(recoEle_pt[i],recoEle_eta[i],recoEle_phi[i],recoEle_mass[i]);
      if(i==1)recolep3->SetPtEtaPhiM(recoEle_pt[i],recoEle_eta[i],recoEle_phi[i],recoEle_mass[i]);

    }


  for(unsigned int i=0; i<recoEle_pt_un.size();i++){
  tree_.recoEle_pt_un.push_back(recoEle_pt_un[i]);

  }

  for(unsigned int i=0; i<recoMu_pt_un.size();i++){
  tree_.recoMu_pt_un.push_back(recoMu_pt_un[i]);

  }
    TLorentzVector recoTQ(*recolep0+*recolep1+*recolep2+*recolep3);
    tree_.TQ_recoMass=recoTQ.M();
    
    tree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
TQGenAnalyzer::beginJob()
{
  std::cout << "Starting job" << std::endl;

  // --- set up output tree
  tree = fs->make<TTree>("tree","tree");
  tree->Branch("run", &tree_.run, "run/I");
  tree->Branch("event", &tree_.event, "event/I");
  tree->Branch("lumi", &tree_.lumi, "lumi/I");
  tree->Branch("nEle", &tree_.nEle, "nEle/I");
  tree->Branch("nMu", &tree_.nMu, "nMu/I");
  tree->Branch("nEleReco", &tree_.nEleReco, "nEleReco/I");
  tree->Branch("nMuReco", &tree_.nMuReco, "nMuReco/I"); 
  tree->Branch("TQ_genMass", &tree_.TQ_genMass, "TQ_genMass/F");
  tree->Branch("TQ_recoMass", &tree_.TQ_recoMass, "TQ_recoMass/F");
  tree->Branch("genLep_pt",     &tree_.genLep_pt);
  tree->Branch("genLep_mass",      &tree_.genLep_mass);
  tree->Branch("genLep_eta",&tree_.genLep_eta);
  tree->Branch("genLep_phi",&tree_.genLep_phi);
  tree->Branch("genLep_pdgId",&tree_.genLep_pdgId);
  tree->Branch("genMom_pt",     &tree_.genMom_pt);
  tree->Branch("genMom_mass",      &tree_.genMom_mass);
  tree->Branch("genMom_eta",&tree_.genMom_eta);
  tree->Branch("genMom_phi",&tree_.genMom_phi);
  tree->Branch("genMom_pdgId",&tree_.genMom_pdgId);
  tree->Branch("recoMu_pt",     &tree_.recoMu_pt);
  tree->Branch("recoMu_pt_un",     &tree_.recoMu_pt_un);
  tree->Branch("recoMu_mass",      &tree_.recoMu_mass);
  tree->Branch("recoMu_eta",&tree_.recoMu_eta);
  tree->Branch("recoMu_phi",&tree_.recoMu_phi);
  tree->Branch("recoMu_dR",&tree_.recoMu_dR);
  tree->Branch("recoEle_pt",     &tree_.recoEle_pt);
  tree->Branch("recoEle_pt_un",     &tree_.recoEle_pt_un);
  tree->Branch("recoEle_mass",      &tree_.recoEle_mass);
  tree->Branch("recoEle_eta",&tree_.recoEle_eta);
  tree->Branch("recoEle_phi",&tree_.recoEle_phi);
  tree->Branch("recoEle_charge", &tree_.recoEle_charge);
  tree->Branch("recoEle_dR",&tree_.recoEle_dR);
  tree->Branch("recoMu_charge", &tree_.recoMu_charge);

}

// ------------ method called once each job just after ending the event loop  ------------
void
TQGenAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TQGenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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


void TQGenAnalyzer::initTreeStructure()
{
  
  tree_.run= -500.;
  tree_.lumi= -500.;
  tree_.event= 0.;


}

void TQGenAnalyzer::clearVectors()
{

  tree_.genLep_pt.clear();
  tree_.genLep_mass.clear();
  tree_.genLep_eta.clear();
  tree_.genLep_phi.clear();
  tree_.genLep_pdgId.clear();
  tree_.genMom_pt.clear();
  tree_.genMom_mass.clear();
  tree_.genMom_eta.clear();
  tree_.genMom_phi.clear();
  tree_.genMom_pdgId.clear();

  tree_.recoMu_pt.clear();
  tree_.recoMu_mass.clear();
  tree_.recoMu_eta.clear();
  tree_.recoMu_phi.clear();
  tree_.recoMu_dR.clear();
  tree_.recoMu_pt_un.clear();

  tree_.recoEle_pt.clear();
  tree_.recoEle_mass.clear();
  tree_.recoEle_eta.clear();
  tree_.recoEle_phi.clear();
  tree_.recoEle_dR.clear();
  tree_.recoEle_pt_un.clear();
  tree_.recoEle_charge.clear();
  tree_.recoMu_charge.clear();
}
//define this as a plug-in
DEFINE_FWK_MODULE(TQGenAnalyzer);
