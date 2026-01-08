#include "PhysicsListPC.hh"

#ifdef ROOT5
#include "Reflex/PluginService.h"

PLUGINSVC_FACTORY(PhysicsListPC,G4VUserPhysicsList*())


#else

#include "SEAL_Foundation/PluginManager/PluginManager/ModuleDef.h"
#include "GamosCore/GamosPhysics/PhysicsList/include/GmPhysicsFactory.hh"


DEFINE_SEAL_MODULE ();
DEFINE_GAMOS_PHYSICS (PhysicsListPC);

#endif