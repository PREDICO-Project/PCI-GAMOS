#ifndef PhysicsListPC_HH
#define PhysicsListPC_HH 

#include "G4VModularPhysicsList.hh"
#include "globals.hh"


class PhysicsListPC: public G4VModularPhysicsList
{
  public:

   PhysicsListPC();
 	virtual ~PhysicsListPC();
    
   void ConstructProcess();
   void ConstructParticle();
 	void SetCuts();
   void AddPhaseContrast();
 	void ReplacePhysicsList(const G4String& name);


 private:

	G4String fEmName;
   G4VPhysicsConstructor* fEMPhysicsList; 
 
};

#endif