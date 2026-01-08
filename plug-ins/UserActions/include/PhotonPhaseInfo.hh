#ifndef PhotonPhaseInfo__hh
#define PhotonPhaseInfo__hh

#include "G4VUserTrackInformation.hh"

class PhotonPhaseInfo : public G4VUserTrackInformation {
public:
    PhotonPhaseInfo() : phaseAccumulated(0.0) {}
    void SetPhase(G4double phase) {phaseAccumulated = phase;}
    G4double GetPhase() const {return phaseAccumulated;}

private:
    G4double phaseAccumulated;
};

#endif
