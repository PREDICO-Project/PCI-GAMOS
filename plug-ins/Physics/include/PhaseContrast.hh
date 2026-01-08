#ifndef PhaseContrast_HH
#define PhaseContrast_HH

// Include Geant4 and other necessary headers
#include "globals.hh"
#include "templates.hh"
#include "geomdefs.hh"
#include "Randomize.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4TransportationManager.hh"
#include "G4Gamma.hh"

// The PhaseContrast class inherits from G4VDiscreteProcess, which is a type of "process" in Geant4.
// A "process" is a set of physical or geometric actions that can occur to a particle
// during its trajectory (for example, interaction with matter, decay, refraction, etc.)

class PhaseContrast : public G4VDiscreteProcess
{

public:

    PhaseContrast(const G4String& processName = "PhaseContrast", G4ProcessType type = fOptical);

    ~PhaseContrast();

public: 
    G4bool IsApplicable(const G4ParticleDefinition& aParticleType);

    G4double GetMeanFreePath(const G4Track& ,
                    G4double ,
                    G4ForceCondition* condition);

    G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step&  aStep);

    G4double GetRindex(G4Material *Material, const G4double Energy);

private:
    // Snell_Law: Applies Snell's law to the photon's direction.
    void Snell_Law();

    void G4Swap(G4double* a, G4double* b) const;
    void G4Swap(G4Material* a, G4Material* b) const;
    void G4VectorSwap(G4ThreeVector* vec1, G4ThreeVector* vec2) const;
    void ReadFile(const std::string& materialName, std::vector<G4double>& Energy_vec, std::vector<G4double>& delta_vec);

    std::vector<G4Material*> GetAllMaterials();

    void InitializeRINDEX(const std::vector<G4Material*>& materials);

    // For materials: calculates and assigns the combined RINDEX.
    void FillRINDEX(G4Material* material);

    G4double InterpolateDelta(G4double energy, const std::vector<G4double>& Energy_vec, const std::vector<G4double>& delta_vec);
    
    void GetPerturbedNormal(G4ThreeVector& normal, double sigma_rad);

private:
    // Materials before and after the boundary
    G4Material* Material1;
    G4Material* Material2;

    // Material names (not essential if not used directly)
    G4String Material1Name;
    G4String Material2Name;

    // Refractive indices of the two materials (Rindex1 and Rindex2)
    G4double Rindex1;
    G4double Rindex2;

    // Variables used in the calculation of incidence and refraction angles.
    G4double cos1, cos2, sin1, sin2;
    G4double sigma;

    // Photon momentum (thePhotonMomentum)
    G4double thePhotonMomentum;
    // Photon direction before (OldMomentum) and after (NewMomentum) the change
    G4ThreeVector OldMomentum;
    G4ThreeVector NewMomentum;

    // Surface normals
    G4ThreeVector theGlobalNormal;
    G4ThreeVector theFacetNormal;

    // Global list of all environment materials
    std::vector<G4Material*> allMaterials;

    // Geometric tolerance
    G4double kCarTolerance;
};

// Inline method: indicates that it only applies to gammas (G4Gamma)
inline
G4bool PhaseContrast::IsApplicable(const G4ParticleDefinition& aParticleType)
{
   // Returns true if the particle is a gamma photon
   return  ( &aParticleType == G4Gamma::Gamma() );
}

// Swaps the value of two G4double
inline
void PhaseContrast::G4Swap(G4double* a, G4double* b) const
{
  G4double temp = *a;
  *a = *b;
  *b = temp;
}

// Swaps pointers to materials
inline
void PhaseContrast::G4Swap(G4Material* a, G4Material* b) const
{
   G4Material* temp = a;
   a = b;
   b = temp;
}

// Swaps two 3D vectors
inline
void PhaseContrast::G4VectorSwap(G4ThreeVector* vec1, G4ThreeVector* vec2) const
{
  G4ThreeVector temp = *vec1;
  *vec1 = *vec2;
  *vec2 = temp;
}

#endif
