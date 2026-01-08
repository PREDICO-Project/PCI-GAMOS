#ifndef GmGetImage__hh
#define GmGetImage__hh

#include "GamosCore/GamosUserActionMgr/include/GmUserRunAction.hh"
#include "GamosCore/GamosUserActionMgr/include/GmUserEventAction.hh"
#include "GamosCore/GamosUserActionMgr/include/GmUserSteppingAction.hh"
#include "GamosCore/GamosUserActionMgr/include/GmUserTrackingAction.hh"
#include "GamosCore/GamosSD/include/GmHitsEventMgr.hh"

class GmGetImage : public GmUserRunAction, public GmUserEventAction, public GmUserSteppingAction
{
public:
    GmGetImage();
    ~GmGetImage(){};

    virtual void BeginOfRunAction( const G4Run* run);
    virtual void UserSteppingAction(const G4Step* aStep);    
    virtual void EndOfRunAction(const G4Run*);
    virtual void EndOfEventAction( const G4Event* );
    

private:
        //metodos
        void GetParams();

        //void InitializeMap(std::map<double, double>& map, const std::string& file_name);
        //double GetMapValue(std::map<double, double> map, double stepEnergy);

        void InitializeVector(std::vector<std::pair<double, double>>& vec, const std::string& file_name);
        double GetInterpolatedValue(const std::vector<std::pair<double, double>>& vec, double stepEnergy );

        G4bool AcceptStep(const G4Step*);
        G4float CalculateGridTransmissionProbability(const G4ThreeVector& position, const G4ThreeVector& direction, const G4Step* aStep);        
        void FillMatrix_MCD(GmRecHit* hit);
        void FillMatrix_Virtual(const G4Step* aStep);
      
        int EnergyToCharge(double energyDeposited);
        void WriteText(const std::string& name);
        void WriteMHD(const std::string& name);
        void WriteDCM(const std::string& name);
        G4Material* LoadCustomMaterial(const std::string& materialName, const std::string& filePath);

        //Atributos
        std::vector<std::pair<double, double>> gridInterspaceMu, gridStripMu, efficiencyMap;
        std::vector<std::vector<double>> output;
        G4int numPixelsX, numPixelsY, totalEvents;
        G4double pixelSizeX, pixelSizeY, zstop, pairCreationEnergyEV,SensibilityOffset, SensibilityFactor, detectorThickness;
        G4String outputFilename, apply_grid, detectorModel, outputSel, outputFolder, outputFormat, detectorMaterial;
        G4float gridRatio, sdd, gridFreq, gridStripThickness, gap;
        G4Material* DetectorMaterial;
        
protected:
  G4ThreeVector postR;
  G4ThreeVector preR;
};

#endif
