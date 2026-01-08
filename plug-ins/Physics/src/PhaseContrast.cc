#include "G4ios.hh"
#include "G4Gamma.hh"
#include "PhaseContrast.hh"
#include "G4GeometryTolerance.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4TransportationManager.hh"
#include "G4RandomDirection.hh"
//#include "CLHEP/Random/G4RandGauss.h"

#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <algorithm> // for lower_bound

// Constructor: initializes variables and calls the function to get all materials and their indices.
PhaseContrast::PhaseContrast(const G4String& processName, G4ProcessType type) 
  : G4VDiscreteProcess(processName, type)
{
  // Initial values
  Material1 = NULL;
  Material2 = NULL;
  Rindex1 = Rindex2 =1.;
  cos1 = cos2 = sin1 = sin2 = 0.;
  sigma = 0.1;

#ifdef DEBUG
  G4cout << " PhaseContrast initialize" << G4endl;
#endif

  // Get global geometric tolerance
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  // Get all defined materials in the environment
  allMaterials = GetAllMaterials();

  // Initialize RINDEX for each material (reads files, calculates deltas, etc.)
  InitializeRINDEX(allMaterials);
}

// Destructor (does nothing special)
PhaseContrast::~PhaseContrast() {}

// PostStepDoIt is called after each photon step. If it crosses a boundary, Snell's law is applied.
G4VParticleChange *PhaseContrast::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep){
    
    // Initialize the particle change object (aParticleChange)
    aParticleChange.Initialize(aTrack);

    // Keep the same photon velocity
    aParticleChange.ProposeVelocity(aTrack.GetVelocity());

    // Get start and end points of the step
    G4StepPoint* pPreStepPoint = aStep.GetPreStepPoint();    
    G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

    // If the end of the step is not a geometric boundary, do nothing special.
    if (pPostStepPoint->GetStepStatus() != fGeomBoundary) {
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep); 
    }

    // If the step is extremely small, do nothing either
    if (aTrack.GetStepLength() <= kCarTolerance/2) {
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }

    // Get materials before and after the step
    Material1 = pPreStepPoint->GetMaterial(); 
    Material2 = pPostStepPoint->GetMaterial();

    // Get the dynamic particle (photon)
    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

    // Total photon momentum
    thePhotonMomentum = aParticle->GetTotalMomentum();

    //G4cout << "Pre Snell Photon Momentum: " << thePhotonMomentum << G4endl;

    // Previous photon direction
    OldMomentum = aParticle->GetMomentumDirection();

    // Get optical properties of the first material
    G4MaterialPropertiesTable* aMaterial1PropertiesTable = Material1->GetMaterialPropertiesTable();
    G4MaterialPropertyVector* RindexProp1 = aMaterial1PropertiesTable->GetProperty("RINDEX");

    G4bool bIsOutOfRange;
    // Refractive index of the first material
    Rindex1 = RindexProp1->GetValue(thePhotonMomentum, bIsOutOfRange);

    // Refractive index of the second material
    G4MaterialPropertiesTable* aMaterial2PropertiesTable = Material2->GetMaterialPropertiesTable();
    G4MaterialPropertyVector* RindexProp2 = aMaterial2PropertiesTable->GetProperty("RINDEX");
    Rindex2 = RindexProp2->GetValue(thePhotonMomentum, bIsOutOfRange); 

    // Global position at the end of the step
    G4ThreeVector theGlobalPoint = pPostStepPoint->GetPosition();
    // Navigator for global-local transformations
    G4Navigator * theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    // Convert to local coordinates
    G4ThreeVector theLocalPoint = theNavigator->GetGlobalToLocalTransform().TransformPoint(theGlobalPoint);

    // Get the surface normal
    G4ThreeVector theLocalNormal;
    G4bool valid;
    theLocalNormal = theNavigator->GetLocalExitNormal(&valid);
    GetPerturbedNormal(theLocalNormal, sigma);

    // If valid, invert it
    if (valid) {
        theLocalNormal = -theLocalNormal;
        GetPerturbedNormal(theLocalNormal, sigma);
    } else {
        G4cerr << " Err in PhaseContrast/PostStepDoIt(). "<< G4endl;
    }

    // Convert the normal to global coordinates
    theGlobalNormal = theNavigator->GetLocalToGlobalTransform().TransformAxis(theLocalNormal);

    // To add roughness ¿?
    // G4double deltaX = CLHEP::G4RandGauss::shoot(0,0.01);
    // G4double deltaY = CLHEP::G4RandGauss::shoot(0,0.01);
    // G4double deltaZ = CLHEP::G4RandGauss::shoot(0,0.01);
    // theGlobalNormal.setX(theGlobalNormal.x() + deltaX);
    // theGlobalNormal.setY(theGlobalNormal.y() + deltaY);
    // theGlobalNormal.setZ(theGlobalNormal.z() + deltaZ);

    // If the photon is going in the same direction as the normal, invert the normal
    if (OldMomentum * theGlobalNormal > 0.0) {
        theGlobalNormal = -theGlobalNormal;
    }

    // Apply Snell's law to calculate the new photon direction
    Snell_Law(); 

    // Ensure the vector is normalized
    NewMomentum = NewMomentum.unit();
    thePhotonMomentum = aParticle->GetTotalMomentum();
    //G4cout << "Post Snell Photon Momentum: " << thePhotonMomentum << G4endl;

    // Propose the new photon direction after crossing the boundary
    aParticleChange.ProposeMomentumDirection(NewMomentum);

    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// Apply Snell's law: calculates if there is refraction or total internal reflection.
void PhaseContrast::Snell_Law() {
    G4bool Swap = false;
    G4bool Through = false;

    // If Through were true, we would swap materials, etc., but here it is initially false.
    if (Through) {
        Swap = !Swap;
        Through = false;
        theGlobalNormal = -theGlobalNormal;
        G4Swap(Material1, Material2);
        G4Swap(&Rindex1, &Rindex2);
    }

    theFacetNormal = theGlobalNormal;
    G4double PdotN = OldMomentum * theFacetNormal;
    cos1 = -PdotN; // cosine of the incidence angle

    // Compute sin1 and sin2 based on the refractive indices
    if (std::abs(cos1) < 1.0 - kCarTolerance) {
        sin1 = std::sqrt(1.-cos1*cos1);
        sin2 = sin1 * Rindex1 / Rindex2; 
    } else {
        sin1 = 0.0;
        sin2 = 0.0;
    }

    // If sin2 >= 1.0, total internal reflection occurs
    if (sin2 >= 1.0) {
        if (Swap) Swap = !Swap;
        PdotN = OldMomentum * theFacetNormal;
        NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal; // reflection

    } else if (sin2 < 1.0) {
        // Refraction
        if (cos1 > 0.0) {
            cos2 = std::sqrt(1.-sin2*sin2);
        } else {
            cos2 = -std::sqrt(1.-sin2*sin2);
        }
        Through = true;
        if (sin1 > 0.0) {
            G4double alpha = cos1 - cos2*(Rindex2/Rindex1);
            NewMomentum = OldMomentum + alpha*theFacetNormal;
            NewMomentum = NewMomentum.unit();
        } else {
            // Normal incidence (sin1=0)
            NewMomentum = OldMomentum;
        }   
    }
    OldMomentum = NewMomentum.unit();

    /*/
    G4cout << "Rindex1: " << Rindex1 << G4endl;
    G4cout << "Rindex2: " << Rindex2 << G4endl;
    G4cout << "delta1: "<< 1-Rindex1 << G4endl;
    G4cout << "delta2: "<< 1-Rindex2 << G4endl;
    G4cout << "Material1 " << Material1 << G4endl;
    G4cout << "Material2 " << Material2 << G4endl;
    /*/
}

// GetMeanFreePath: force it to be infinite so that this process is evaluated at boundaries, not at fixed distances.
G4double PhaseContrast::GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* condition) {
    *condition = Forced;
    return DBL_MAX;
}

// Global cache to avoid re-reading files repeatedly.
static std::unordered_map<std::string, std::pair<std::vector<G4double>, std::vector<G4double>>> fileCache;

// ReadFile: reads a file with data (Energy, delta, beta). Only delta is used.
void PhaseContrast::ReadFile(const std::string& materialName, std::vector<G4double>& Energy_vec, std::vector<G4double>& delta_vec) 
{
    std::string name = materialName;
    // If the name starts with G4_, remove that part
    if (name.substr(0, 3) == "G4_") {
        name.erase(0,3);
    }

    // Check if this material was already read (cache)
    auto it = fileCache.find(name);
    if (it != fileCache.end()) {
        Energy_vec = it->second.first;
        delta_vec = it->second.second;
        return;
    }

    // Open the file with complex refractive index data
    std::ifstream infile("data/complex_refractive_index/"+name+".txt");
    if (!infile.is_open()) {
        G4cerr << "Error: Could not open file for: " << name << G4endl;
        return;
    }

    std::string line;
    G4double Energy, delta, beta;

    // Reserve approximate space
    Energy_vec.reserve(1000);
    delta_vec.reserve(1000);

    // Read line by line
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        if (!(iss >> Energy >> delta >> beta)) {
#ifdef DEBUG
            G4cout << "Line format error in file for material: " << name << G4endl;
#endif
            continue;
        }
        // Store Energy in keV, delta directly
        Energy_vec.push_back(Energy * keV);
        delta_vec.push_back(delta);
    }
    infile.close();

    // Store in the cache
    fileCache[name] = std::make_pair(Energy_vec, delta_vec);
}

// Get all available materials
std::vector<G4Material*> PhaseContrast::GetAllMaterials()
{
  std::vector<G4Material*> mates;
  const G4MaterialTable* matTab = G4Material::GetMaterialTable();
  mates.reserve(matTab->size());

  // Rename variable to 'mat' to avoid conflict with CLHEP's 'm'
  for(const auto& mat : *matTab) {
    mates.push_back(mat);
  }
  return mates;
}

// Initialize RINDEX for each material by calling FillRINDEX_el or FillRINDEX_compound as appropriate.
void PhaseContrast::InitializeRINDEX(const std::vector<G4Material*>& materials)
{
  size_t nMats = materials.size();
  for (size_t ii = 0; ii < nMats; ii++ ) {
      G4Material* mate = materials[ii];

#ifdef DEBUG
      G4cout << "Mat Vector: " << mate->GetElementVector()->size() << G4endl;
#endif

      FillRINDEX(mate);
      
  }
}


// If the material is a compound: calculate a mixture of deltas of its elements.
void PhaseContrast::FillRINDEX(G4Material* material)
{
    G4String materialName = material->GetName();
    G4double densityCompound = material->GetDensity()/(g/cm3);
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4double* fractions = material->GetFractionVector();

    G4int nElements = material->GetNumberOfElements();

    // Define a unified energy range for the calculation
    G4double minEnergy = 1.0 *keV;
    G4double maxEnergy = 150.0 * keV;
    int sizeUnifiedEnergies = 250;
    std::vector<G4double> unifiedEnergies;
    unifiedEnergies.reserve(sizeUnifiedEnergies);
    G4double step = (maxEnergy - minEnergy) / (sizeUnifiedEnergies - 1);

    for (int i = 0; i < sizeUnifiedEnergies; i++){
        unifiedEnergies.push_back(minEnergy + i * step);
    }

    // Vectors to accumulate the average delta of the compound
    //std::vector<G4double> fractionSum(sizeUnifiedEnergies,0.0);
    G4double totalMassFraction = 0.0;
    std::vector<G4double> massFractions(nElements,0.0);
    std::vector<G4double> deltaSum(sizeUnifiedEnergies,0.0);
    

    // For each element, read its data and add its contribution
    for( int ii = 0; ii < nElements; ii++ ) {
      const G4Element* element = (*elementVector)[ii];
      G4String elementName = element->GetName(); 
      G4double fraction = fractions[ii];

      // Retrieve Density of the element 

      G4Material* elementMaterial = G4Material::GetMaterial("G4_"+elementName);
      if (!elementMaterial){
        G4NistManager* nist = G4NistManager::Instance();
        elementMaterial = nist->FindOrBuildMaterial("G4_"+elementName);
      }

      G4double elementDensity = elementMaterial->GetDensity() / (g/cm3);

      std::vector<G4double> Energy_vec, delta_vec;
      ReadFile(elementName, Energy_vec, delta_vec);

      for (int j = 0; j < sizeUnifiedEnergies; ++j){
        G4double energyVal = unifiedEnergies[j];
        // Interpolate delta at this energy
        G4double interpolatedDelta = InterpolateDelta(energyVal, Energy_vec, delta_vec);
        deltaSum[j] += interpolatedDelta * (fraction * densityCompound / elementDensity);
      }
    }


    // Now calculate the effective delta for the compound
    std::vector<G4double> Rindex_vec;
    Rindex_vec.reserve(sizeUnifiedEnergies);
    for (int k = 0; k < sizeUnifiedEnergies; ++k){
        G4double averagedDelta = deltaSum[k];
        Rindex_vec.push_back(1.0 - averagedDelta);
    }
    
    // Material Properties Table
    G4MaterialPropertiesTable* matPropTbl = new G4MaterialPropertiesTable();
    matPropTbl->AddProperty("RINDEX",unifiedEnergies.data(),Rindex_vec.data(),(G4int)unifiedEnergies.size());
    material->SetMaterialPropertiesTable(matPropTbl);
      

#ifdef DEBUG
    G4cout << "Material: " << materialName << " RINDEX property added" << G4endl;
#endif
}

// InterpolateDelta: given a vector of energies and deltas, and a specific energy, find delta by linear interpolation.
G4double PhaseContrast::InterpolateDelta(G4double energy, const std::vector<G4double>& Energy_vec, const std::vector<G4double>& delta_vec)
{
    G4int size = (G4int)Energy_vec.size();

    if (size == 0){
        G4cerr << "Empty File" << G4endl;
        return 0.0;
    }

    // If energy is out of range, return the nearest boundary value
    if (energy <= Energy_vec.front()) return delta_vec.front();
    if (energy >= Energy_vec.back()) return delta_vec.back();

    // Find where 'energy' fits in
    auto lower = std::lower_bound(Energy_vec.begin(), Energy_vec.end(), energy);
    if (lower == Energy_vec.begin()) return delta_vec.front();
    if (lower == Energy_vec.end()) return delta_vec.back();

    // Linear interpolation
    G4int i = (G4int)(lower - Energy_vec.begin()) - 1;
    G4double E1 = Energy_vec[i];
    G4double E2 = Energy_vec[i+1];
    G4double d1 = delta_vec[i];
    G4double d2 = delta_vec[i+1];

    G4double slope = (d2 - d1)/(E2 - E1);
    return d1 + slope*(energy - E1);
}

// GetRindex: currently unused. Left in case it is needed.
G4double PhaseContrast::GetRindex(G4Material* Material, const G4double Energy)
{
   (void)Material;
   (void)Energy;
   return 1.0;
}


void PhaseContrast::GetPerturbedNormal(G4ThreeVector& normal, double sigma_rad) {
    // Random vector to rotate the normal
    G4ThreeVector rand_vec = G4RandomDirection();
    
    // Projection 
    G4ThreeVector tangent_component = rand_vec - (rand_vec.dot(normal)) * normal;
    tangent_component = tangent_component.unit();

    // Deviation angle
    double theta = G4RandGauss::shoot(0.0, sigma_rad);

    // Rotate the normal
    normal = std::cos(theta) * normal + std::sin(theta) * tangent_component;
    normal = normal.unit();
}
