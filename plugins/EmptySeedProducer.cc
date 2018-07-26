#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"

using namespace edm;
using namespace std;
using namespace reco;

class EmptySeedProducer: public edm::stream::EDProducer<> {
 public:
  explicit EmptySeedProducer(const edm::ParameterSet& cfg) {
		produces<ElectronSeedCollection>();
	}
  
private:
	void beginRun(const edm::Run & run,const edm::EventSetup&) override {}
	void produce(edm::Event& iEvent, const edm::EventSetup&) override {
		auto output_preid = std::make_unique<ElectronSeedCollection>();
		iEvent.put(std::move(output_preid));
	}	
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EmptySeedProducer);

