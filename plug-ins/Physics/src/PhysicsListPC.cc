#include "PhysicsListPC.hh"
#include "PhaseContrast.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmLivermorePolarizedPhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmStandardPhysicsSS.hh"

#include "G4ParticleDefinition.hh"

#include "GamosCore/GamosPhysics/PhysicsList/include/GmPhysicsGammaStandard.hh"
#include "GamosCore/GamosPhysics/PhysicsList/include/GmPhysicsGammaLowEner.hh"
#include "GamosCore/GamosPhysics/PhysicsList/include/GmPhysicsGammaPenelope.hh"
#include "GamosCore/GamosPhysics/PhysicsList/include/GmPhysicsElectronStandard.hh"
#include "GamosCore/GamosPhysics/PhysicsList/include/GmPhysicsElectronLowEner.hh"
#include "GamosCore/GamosPhysics/PhysicsList/include/GmPhysicsElectronPenelope.hh"
#include "GamosCore/GamosPhysics/PhysicsList/include/GmPhysicsPositronStandard.hh"
#include "GamosCore/GamosPhysics/PhysicsList/include/GmPhysicsPositronPenelope.hh"
#include "GamosCore/GamosPhysics/PhysicsList/include/GmPhysicsDecay.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4ProcessManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4BosonConstructor.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

//#include "G4EmProcessOptions.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "GamosCore/GamosPhysics/PhysicsList/include/GmPhysicsGammaStandard_XSChange.hh"
#include "GamosCore/GamosBase/Base/include/GmParameterMgr.hh"


#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4ProcessManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListPC::PhysicsListPC() : G4VModularPhysicsList()
{
  
  fEmName = G4String("local");
  G4cout << "<<< User Defined Physics List simulation engine: Phase Contrast included!"<<G4endl;

  defaultCutValue = 0.7*CLHEP::mm;
  G4int ver = 1;
  //SetVerboseLevel(ver);

   // EM Standar Physics
  fEMPhysicsList = new G4EmStandardPhysics_option4();
  RegisterPhysics(fEMPhysicsList);
  G4EmParameters::Instance()->SetBuildCSDARange(true);
  G4cout << "fEMPhysicsList & fPCPhysicsList defined" <<G4endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListPC::~PhysicsListPC()
{
  delete fEMPhysicsList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PhysicsListPC::SetCuts()
{
  // Use default cut values gamma and e processes
  SetCutsWithDefault();   
  G4cout << "<<< SetCutsWithDefault:" << defaultCutValue <<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListPC::ConstructParticle()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  
  // gamma
  G4Gamma::GammaDefinition();
  
  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();

  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();  

// mesons
  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();

// barions
  G4BaryonConstructor bConstructor;
  bConstructor.ConstructParticle();

// ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();
  //fEMPhysicsList->ConstructParticle();
  
  G4cout << "<<< Particles cosntructed"<<G4endl;
}


void PhysicsListPC::ConstructProcess()
{
  AddTransportation();
  fEMPhysicsList->ConstructProcess();
  //fPCPhysicsList->ConstructProcess();
  AddPhaseContrast();
  
}

void PhysicsListPC::AddPhaseContrast(){
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleTable::G4PTblDicIterator* particleIterator = particleTable->GetIterator();

  particleIterator->reset();

  while ((*particleIterator)())
  {
    G4ParticleDefinition* particle = particleIterator->value();

    if (particle == G4Gamma::Gamma()) 
    {
      G4ProcessManager* pManager = particle->GetProcessManager();
      pManager -> AddDiscreteProcess(new PhaseContrast);
      G4cout << "PhaseContrast Discrete Process added" << G4endl;
    }   
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListPC::ReplacePhysicsList(const G4String& name)
{
  if (verboseLevel>-1) {
    G4cout << "PhysicsListPC::AddGmEMStandardPhysics: <" << name << ">" << G4endl;
  }

  if (name == fEmName) return;

  if (name == "emstandard_opt0"
      || name == "emstandard_option0") {
    
    fEmName = name;
    delete fEMPhysicsList;
    fEMPhysicsList = new G4EmStandardPhysics();

  } else if (name == "emstandard_opt1" 
	     || name == "emstandard_option1") {

    fEmName = name;
    delete fEMPhysicsList;
    fEMPhysicsList = new G4EmStandardPhysics_option1();
    
  } else if (name == "emstandard_opt2"
    || name == "emstandard_option2") {

    fEmName = name;
    delete fEMPhysicsList;
    fEMPhysicsList = new G4EmStandardPhysics_option2();
    
  } else if (name == "emstandard_opt3"
	     || name == "emstandard_option3") {
    
    fEmName = name;
    delete fEMPhysicsList;
    fEMPhysicsList = new G4EmStandardPhysics_option3();

  } else if (name == "emstandard_opt4"
	     || name == "emstandard_option4") {

    fEmName = name;
    delete fEMPhysicsList;
    fEMPhysicsList = new G4EmStandardPhysics_option4();
    
  } else if (name == "emlivermore") {

    fEmName = name;
    delete fEMPhysicsList;
    fEMPhysicsList = new G4EmLivermorePhysics();
    
  } else if (name == "emlivermorepolarized") {

    fEmName = name;
    delete fEMPhysicsList;
    fEMPhysicsList = new G4EmLivermorePolarizedPhysics();
    
  } else if (name == "emlowEP") {

    fEmName = name;
    delete fEMPhysicsList;
    fEMPhysicsList = new G4EmLowEPPhysics();
    
  } else if (name == "empenelope") {

    fEmName = name;
    delete fEMPhysicsList;
    fEMPhysicsList = new G4EmPenelopePhysics();
    
  } else if (name == "emstandard_GS") {

    fEmName = name;
    delete fEMPhysicsList;
    fEMPhysicsList = new G4EmStandardPhysicsGS();
    
  } else if (name == "emstandard_WVI") {

    fEmName = name;
    delete fEMPhysicsList;
    fEMPhysicsList = new G4EmStandardPhysicsWVI();
    
  } else if (name == "emstandard_SS") {

    fEmName = name;
    delete fEMPhysicsList;
    fEMPhysicsList = new G4EmStandardPhysicsSS();
    
  } else {

    G4cout << "PhysicsListPC::AddGmEMStandardPhysics: <" << name << ">"
           << " is not defined"
           << G4endl;
	  exit(1);
  }

}
