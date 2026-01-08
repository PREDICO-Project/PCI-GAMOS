#ifndef GmGetWavefront__hh
#define GmGetWavefront__hh

#include "GamosCore/GamosUserActionMgr/include/GmUserRunAction.hh"
#include "GamosCore/GamosUserActionMgr/include/GmUserEventAction.hh"
#include "GamosCore/GamosUserActionMgr/include/GmUserSteppingAction.hh"
#include "GamosCore/GamosUserActionMgr/include/GmUserTrackingAction.hh"
#include "GamosCore/GamosGenerator/include/GmGenerator.hh"

class GmGetWavefront : public GmUserRunAction, public GmUserSteppingAction, public GmGenerator
{
public:
    GmGetWavefront();
    ~GmGetWavefront(){};

    virtual void BeginOfRunAction( const G4Run* run);
    //virtual void BeginOfEventAction(const G4Event* event);
    virtual void UserSteppingAction(const G4Step* aStep);    
    virtual void EndOfRunAction(const G4Run* run);
    
private:
        //methods
        void GetParams();
        void FillWavefront(const G4Step* aStep);
        void WriteWavefront();
        void PropagateWavefrontFresnel();
        void CheckNyquist(double px, double py, double wavelength, double dz);
        void ComputeIntensity();
        void WriteToMHD(const std::string& name, std::vector<std::vector<double>> data, double pSizeX, double pSizeY);
        void ApplyAbsorptionGrating(double period, const G4ThreeVector& displacement, const G4ThreeVector& axis, double theta);
        G4ThreeVector Rotate(const G4ThreeVector& point, const G4ThreeVector& axis, double theta);
        void DownsampleIntensity();
        std::vector<std::complex<double>> PadWavefrontReplication(const std::vector<std::vector<double>>& Real, const std::vector<std::vector<double>>& Imag, int padX, int padY);
        void ApplyGaussianBlur(std::vector<std::vector<double>>& waveReal,std::vector<std::vector<double>>& waveImag,int NX, int NY);

           
        G4int numPixelsX, numPixelsY, totalEvents, detnumPixelsX, detnumPixelsY, particleCount, numPixelsX_corr, numPixelsY_corr;
        G4double pixelSizeX, pixelSizeY, zstop, propagationDistance, periodG2, rotAngle, detpixelSizeX, detpixelSizeY, magnification, eventsSimulated, energyAcummulated, effectiveEnergy;
        //G4double* phaseAccumulated;
        G4String outputFilename, outputSel, outputFolder, outputFormat, distDirectionName;
        std::vector<std::vector<double>> waveReal, waveImag, WaveTot,intensity, phaseVec, grid, downsampledIntensity, sumPhase, sumAmp;
        G4bool debugMode, doTalbot;
        G4ThreeVector disp, rotAxis;
        std::vector<G4double> dispValues, dispVP;
        std::vector<G4double> axisValues;
        
protected:
  G4ThreeVector postR;
  G4ThreeVector preR;
};



#endif
