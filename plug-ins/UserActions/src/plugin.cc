
#include "GmGetImage.hh"
#include "GmGetWavefront.hh"
//#include "GmHuygens.hh"


#ifdef ROOT5
#include "Reflex/PluginService.h"

PLUGINSVC_FACTORY(GmGetImage,GmUserAction*())
PLUGINSVC_FACTORY(GmGetWavefront,GmUserAction*())
//PLUGINSVC_FACTORY(GmHuygens,GmUserAction*())

#else 

#include "PluginManager/ModuleDef.h"
#include "GamosCore/GamosUserActionMgr/include/GmUserActionFactory.hh"
#include "GamosCore/GamosBase/Base/include/GmFilterFactory.hh"

DEFINE_SEAL_MODULE ();

DEFINE_GAMOS_USER_ACTION(GmGetImage);
DEFINE_GAMOS_USER_ACTION(GmGetWavefront);
//DEFINE_GAMOS_USER_ACTION(GmHuygens);
#endif
