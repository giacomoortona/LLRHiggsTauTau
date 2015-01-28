/** \class EleFiller
 *
 *  No description available.
 *
 *  $Date: 2013/05/24 15:42:42 $
 *  $Revision: 1.28 $
 *  \author N. Amapane (Torino)
 *  \author S. Bolognesi (JHU)
 *  \author C. Botta (CERN)
 *  \author S. Casasso (Torino)
 *  \author G. Ortona (LLR)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>

//#include <DataFormats/PatCandidates/interface/Electron.h>
//#include "DataFormats/VertexReco/interface/Vertex.h"
#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CutSet.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LeptonIsoHelper.h>
//#include "BDTId.h"

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;


//bool recomputeBDT = false;

class EleFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit EleFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~EleFiller(){
    //delete bdt;
  };  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  const edm::InputTag theCandidateTag;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > theGenTag ;
  int sampleType;
  int setup;
  const StringCutObjectSelector<pat::Electron, true> cut;
  const CutSet<pat::Electron> flags;
  //BDTId* bdt;
};


EleFiller::EleFiller(const edm::ParameterSet& iConfig) :
  theCandidateTag(iConfig.getParameter<InputTag>("src")),
  theGenTag(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genCollection"))),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  cut(iConfig.getParameter<std::string>("cut")),
  flags(iConfig.getParameter<ParameterSet>("flags"))//, 
  //bdt(0)
{
  //if (recomputeBDT) bdt = new BDTId;
  produces<pat::ElectronCollection>();
}


void
EleFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  

  // Get leptons and rho
  edm::Handle<pat::ElectronRefVector> electronHandle;
  iEvent.getByLabel(theCandidateTag, electronHandle);

  InputTag theRhoTag = LeptonIsoHelper::getEleRhoTag(sampleType,setup);
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(theRhoTag, rhoHandle); 
  double rho = *rhoHandle;

  edm::Handle<vector<Vertex> >  vertexs;
  iEvent.getByLabel("goodPrimaryVertices",vertexs);

  // Output collection
  auto_ptr<pat::ElectronCollection> result( new pat::ElectronCollection() );

  for (unsigned int i = 0; i< electronHandle->size(); ++i){
    //---Clone the pat::Electron
    pat::Electron l(*((*electronHandle)[i].get()));
    
    //--- PF ISO
    float PFChargedHadIso   = l.chargedHadronIso();
    float PFNeutralHadIso   = l.neutralHadronIso();
    float PFPhotonIso       = l.photonIso();

    //float fSCeta = fabs(l.eta()); 
    float fSCeta = fabs(l.superCluster()->eta());

    float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, l);

    //--- SIP, dxy, dz
    float IP      = fabs(l.dB(pat::Electron::PV3D));
    float IPError = l.edB(pat::Electron::PV3D);
    float SIP     = IP/IPError;

    float dxy = 999.;
    float dz  = 999.;
    const Vertex* vertex = 0;
    if (vertexs->size()>0) {
      vertex = &(vertexs->front());
      dxy = (l.gsfTrack()->dxy(vertex->position()));
      dz  = (l.gsfTrack()->dz(vertex->position()));
    } 

    
    // BDT value from PAT prodiuction
    float BDT = 0;
    //if (recomputeBDT) {
    //  BDT = bdt->compute(l);
    //} else {
    //BDT = l.electronID("mvaNonTrigV0");//RH
    //}
    

    float pt = l.pt();
    bool isBDT = (pt <= 10 && (( fSCeta < 0.8 && BDT > 0.47)  ||
			       (fSCeta >= 0.8 && fSCeta < 1.479 && BDT > 0.004) ||
			       (fSCeta >= 1.479               && BDT > 0.295))) || 
                 //Moriond13 eID cuts updated for the paper
		 //(pt >  10 && ((fSCeta < 0.8 && BDT > 0.5)  ||
		 //	       (fSCeta >= 0.8 && fSCeta < 1.479 && BDT > 0.12) || 
                 (pt >  10 && ((fSCeta < 0.8 && BDT > -0.34)  ||
			       (fSCeta >= 0.8 && fSCeta < 1.479 && BDT > -0.65) || 
			       (fSCeta >= 1.479               && BDT > 0.6)));
    
	//-- Missing hit  
	int missingHit = l.gsfTrack()->hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS);
    //--- Trigger matching
    //int HLTMatch=0; 
    
    //--- Embed user variables
    l.addUserFloat("PFChargedHadIso",PFChargedHadIso);
    l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso);
    l.addUserFloat("PFPhotonIso",PFPhotonIso);
    l.addUserFloat("combRelIsoPF",combRelIsoPF);
    l.addUserFloat("rho",rho);
    l.addUserFloat("SIP",SIP);
    l.addUserFloat("dxy",dxy);
    l.addUserFloat("dz",dz);
    l.addUserFloat("BDT",BDT);    
    l.addUserFloat("isBDT",isBDT);
    //l.addUserFloat("HLTMatch", HLTMatch);
    l.addUserFloat("missingHit", missingHit);

    //--- MC info
    const reco::GenParticle* genL= l.genParticleRef().get();
    float px=0,py=0,pz=0,E=0,fromH=0;
    if(genL){
      px =genL->p4().Px();
      py =genL->p4().Py();
      pz =genL->p4().Pz();
      E =genL->p4().E();

      //search if it comes from H
      Handle<edm::View<reco::GenParticle> > genHandle;
      iEvent.getByToken(theGenTag, genHandle);
      int nmot = genL->numberOfMothers();
      for (int im = 0; im<nmot&&fromH==0; ++im){
	for(unsigned int ipruned = 0; ipruned< genHandle->size(); ++ipruned){
	  if(ipruned==i)continue;
	  if(userdatahelpers::isAncestor(&(*genHandle)[ipruned],genL)){
	    fromH=1;
	    break;
	  }
	}
      }
    }
    l.addUserFloat("genPx",px);
    l.addUserFloat("genPy",py);
    l.addUserFloat("genPz",pz);
    l.addUserFloat("genE",E);
    l.addUserFloat("fromH",fromH);

//     MCHistoryTools mch(iEvent);
//     if (mch.isMC()) {
//       int MCParentCode = 0;
//       //      int MCParentCode = mch.getParentCode(&l); //FIXME: does not work on cmg
//       l.addUserFloat("MCParentCode",MCParentCode);
//     }

    //--- Check selection cut. Being done here, flags are not available; but this way we 
    //    avoid wasting time on rejected leptons.
    if (!cut(l)) continue;

    //--- Embed flags (ie flags specified in the "flags" pset)
    for(CutSet<pat::Electron>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      l.addUserFloat(flag->first,int((*(flag->second))(l)));
    }

    result->push_back(l);
  }
  iEvent.put(result);
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(EleFiller);

