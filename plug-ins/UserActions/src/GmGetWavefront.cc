#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"
#include "G4GeometryTolerance.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "Randomize.hh"
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <complex>
#include <fftw3.h>

#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmimgle/dcmimage.h"
#include "dcmtk/dcmdata/dcfilefo.h"
#include "dcmtk/dcmdata/dcdeftag.h"

#include "GamosCore/GamosBase/Base/include/GmParameterMgr.hh"
#include "GamosCore/GamosGenerator/include/GmGenerator.hh"
#include "GamosCore/GamosGenerator/include/GmParticleSource.hh"
#include "GamosCore/GamosGenerator/include/GmVGenerDistDirection.hh"

#include "GmGetWavefront.hh"
#include "PhotonPhaseInfo.hh"

constexpr double MIN_U2 = 1.0e-8;
constexpr double PI = CLHEP::pi;
std::complex<double> I(0.,1.);

//std::vector<std::vector<double>> waveReal(numPixelsY, std::vector<double>(numPixelsX, 0.0));
//std::vector<std::vector<double>> waveImag(numPixelsY, std::vector<double>(numPixelsX, 0.0));



//---------------------------------------------------------------------
GmGetWavefront::GmGetWavefront() : GmUserRunAction(), GmUserSteppingAction(), GmGenerator()
{
}
//----------------------------------------------------------------------
void GmGetWavefront::BeginOfRunAction(const G4Run* run) {
    
    GmGenerator* gmGener = (GmGenerator*)(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    std::vector<GmParticleSource*> sources = gmGener->GetSources();
    GmParticleSource* primSource = sources[0];
    //GmVGenerDistTime* distTime = primSource->GetTimeDistribution();
    //GmVGenerDistPosition* distPos = primSource->GetPositionDistribution();
    GmVGenerDistDirection* distDir = primSource->GetDirectionDistribution();
    distDirectionName = distDir->GetName();
    //G4cout << "DistDirectionName: "  << distDirectionName << G4endl;
   
    GetParams();
    totalEvents = run->GetNumberOfEventToBeProcessed();
   
    G4cout << "=======================" << G4endl;
    G4cout << "= Wavefront Retrieval =" << G4endl;
    G4cout << "=======================" << G4endl;
    G4cout << "   Output selection: " << outputSel << G4endl;
    G4cout << "   Output Format: " << outputFormat << G4endl;
    G4cout << "   Output Folder: " << outputFolder << G4endl;
    G4cout << "   Output Filename: " << outputFilename << G4endl;
    G4cout << "   Total Events: " << totalEvents << G4endl;
    G4cout << "=======================" << G4endl;
    G4cout << "=       Params        = " << G4endl;
    G4cout << "=======================" << G4endl; 
    G4cout << "   Z Stop Position: " << zstop << " mm " << G4endl;
    G4cout << "   Virtual Plane Center: (" << dispVP[0]<<", " << dispVP[1] << ") mm" << G4endl;
    G4cout << "   Propagation Distance: " << propagationDistance << " mm " << G4endl;
    G4cout << "   Number of Pixels on X: " << numPixelsX << G4endl;
    G4cout << "   Number of Pixels on Y: " << numPixelsY << G4endl;
    G4cout << "   Pixel Size on X: " << pixelSizeX / CLHEP::um << " um" << G4endl;
    G4cout << "   Pixel Size on Y: " << pixelSizeY / CLHEP::um << " um" << G4endl;
    G4cout << "   Magnification: " << magnification << G4endl;
    G4cout << "==========================" << G4endl;
    if(doTalbot){
        G4cout << "=======================" << G4endl;
        G4cout << "=       Talbot        = " << G4endl;
        G4cout << "=======================" << G4endl; 
        G4cout << "   G2 Period: " << periodG2 / CLHEP::um << " um " << G4endl;
        G4cout << "   G2 displacement: " << disp << G4endl;
        G4cout << "   Rotation axis: " << rotAxis << G4endl;
        G4cout << "   G2 Rotation angle: " << rotAngle / CLHEP::deg << " deg" << G4endl;
        G4cout << "==========================" << G4endl;
    }
}

//----------------------------------------------------------------------
void GmGetWavefront::GetParams(){

    GmParameterMgr* paramMgr = GmParameterMgr::GetInstance();

    // Get detector parameters from GetImage.in file
    detnumPixelsX = paramMgr->GetNumericValue("GetWavefront:NumPixelsX", 1000);
    detnumPixelsY = paramMgr->GetNumericValue("GetWavefront:NumPixelsY", 1000);
    detpixelSizeX = paramMgr->GetNumericValue("GetWavefront:PixelSizeX", 0.01); //mm
    detpixelSizeY = paramMgr->GetNumericValue("GetWavefront:PixelSizeY", 0.01); //mm

    numPixelsX = paramMgr->GetNumericValue("GetWavefront:gridX", 10000);
    numPixelsY = paramMgr->GetNumericValue("GetWavefront:gridY", 10000);
    pixelSizeX = paramMgr->GetNumericValue("GetWavefront:gridSizeX", 0.0001); //mm
    pixelSizeY = paramMgr->GetNumericValue("GetWavefront:gridSizeY", 0.0001); //mm

    zstop = paramMgr->GetNumericValue("GetWavefront:zVirtualPlane",650);
    dispVP = paramMgr->GetVNumericValue("GetWavefront:VirtualPlaneDisplacement", std::vector<G4double>{0,0});
    propagationDistance = paramMgr->GetNumericValue("GetWavefront:PropDistance",650);
    

    // Output File Parameters
    outputFolder = paramMgr->GetStringValue("GetWavefront:ResultsFolder", "output/");	
	outputFilename = paramMgr->GetStringValue("GetWavefront:OutputFilename", "test");
    outputFormat = paramMgr->GetStringValue("GetWavefront:OutputFormat", "Text");
    debugMode = paramMgr->GetBooleanValue("GetWavefront:debugMode", 0);

    // Talbot 
    doTalbot = paramMgr->GetBooleanValue("DoTalbot:Talbot", 0);
    periodG2 = paramMgr->GetNumericValue("DoTalbot:PeriodG2", 0.005); //mm;
    dispValues = paramMgr->GetVNumericValue("DoTalbot:displacement", std::vector<G4double>{0,0,0});
    axisValues = paramMgr->GetVNumericValue("DoTalbot:rotationAxis", std::vector<G4double>{0,0,0});
    rotAngle = paramMgr->GetNumericValue("DoTalbot:rotationAngle", 0);

    if (distDirectionName == "GmGenerDistDirectionConst") {
        magnification = 1.0;
    }
    else{
        magnification = (zstop + propagationDistance) / zstop;
    }

    numPixelsX_corr = static_cast<G4int>(std::round(magnification * numPixelsX));
    numPixelsY_corr = static_cast<G4int>(std::round(magnification * numPixelsY));

    disp.setX(dispValues[0]);
    disp.setY(dispValues[1]);
    disp.setZ(dispValues[2]);
    rotAxis.setX(axisValues[0]);
    rotAxis.setY(axisValues[1]);
    rotAxis.setZ(axisValues[2]);

    // Adjust energy matrix size from deterctor size
    downsampledIntensity.resize(detnumPixelsY, std::vector<double>(detnumPixelsX, 0.0));
    intensity.resize(numPixelsY, std::vector<double>(numPixelsX, 0.0));
    waveReal.resize(numPixelsY, std::vector<double>(numPixelsX, 0.0));
    waveImag.resize(numPixelsY, std::vector<double>(numPixelsX, 0.0));
    //WaveTot.resize(numPixelsY, std::vector<double>(numPixelsX, 0.0));
    phaseVec.resize(numPixelsY, std::vector<double>(numPixelsX, 0.0));
    //sumPhase.resize(numPixelsY, std::vector<double>(numPixelsX, 0.0));
    //sumAmp.resize(numPixelsY, std::vector<double>(numPixelsX, 0.0));

    energyAcummulated = 0.0;
    particleCount = 0;

}

//----------------------------------------------------------------------
void GmGetWavefront::UserSteppingAction(const G4Step* aStep)
{
  G4Track* track = aStep->GetTrack();
  PhotonPhaseInfo* phaseInfo = static_cast<PhotonPhaseInfo*>(track->GetUserInformation());

 
  // Get the phase accumulated of the photon
  if (!phaseInfo){
    phaseInfo = new PhotonPhaseInfo();
    track->SetUserInformation(phaseInfo);
  }

  G4Material* material = aStep->GetPreStepPoint()->GetMaterial();
  

  G4MaterialPropertiesTable* materialProperties = material->GetMaterialPropertiesTable();
   if (!materialProperties){
    G4cerr << "Error: Material " << material->GetName() << "does not have a MaterialProperties" << G4endl;
    return;
  }

  G4MaterialPropertyVector* rindexProp = materialProperties->GetProperty("RINDEX");

  if (rindexProp == nullptr){
    G4cerr << "Error: Material " << material->GetName() << " does not have a RINDEX property" << G4endl;
    return;
  }

  G4bool bIsOutOfRange;
  G4double energy = aStep->GetPreStepPoint()->GetTotalEnergy();
  G4double rindex = rindexProp->GetValue(energy, bIsOutOfRange);

  if (bIsOutOfRange){
    G4cerr << "Error: Energy " << energy << " is out of range for material " << material->GetName() << G4endl;
    return;
  }

  G4double delta = 1.0 - rindex;

  // Wavelength in mm
  G4double wavelength = (CLHEP::h_Planck * CLHEP::c_light / energy);

  //G4cout << "Wavlength: " << wavelength << "mm" << G4endl;
  //G4cout << "Energy: " << energy/CLHEP::eV << "Delta: " << delta << G4endl;

  // Step Length in mm
  G4double stepLength = aStep->GetTrack()->GetStepLength();

  // Phase
  G4double phase = -2.0 * CLHEP::pi * delta * stepLength / wavelength;
  phaseInfo->SetPhase(phaseInfo->GetPhase() + phase);

  // Get Positions at the PreStep and PostStep Points
  preR= aStep->GetPreStepPoint()->GetPosition();
  postR= aStep->GetPostStepPoint()->GetPosition();
  if (debugMode){
  G4cout << "Material: " << material->GetName() << " Delta: " << delta << G4endl;
  G4cout << "Energy (keV): " << energy/CLHEP::keV << " Wavelength (mm): " << wavelength/ CLHEP::mm << " Step (mm): " << stepLength/CLHEP::mm<< G4endl;
  G4cout << "Phase (rad): " << phase << " Accumulated Phase (rad): " << phaseInfo->GetPhase() << G4endl;
  G4cout << "Pre: " << preR << " Post: " << postR << G4endl;
  } 
  // If PreStep position is larger than zStop, the particle has crossed the detector plane.
  if ( preR.z() > zstop) {track->SetTrackStatus(fStopAndKill);}

  // If PreStep position is lower than zStop and PostStep position is larger, the particle is crossing the detector plane.
  if ( (preR.z()-zstop) * (postR.z()-zstop) < 0. || fabs(postR.z()-zstop) < G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() ) {FillWavefront(aStep);}

}

//----------------------------------------------------------------------
void GmGetWavefront::EndOfRunAction(const G4Run* run) {

    eventsSimulated = run->GetNumberOfEvent();

    //WriteToMHD("phase", phaseVec, pixelSizeX , pixelSizeY);
    //WriteToMHD("imag_before", waveImag, pixelSizeX , pixelSizeY);
    //WriteToMHD("real_before", waveReal, pixelSizeX , pixelSizeY);
    //WriteToMHD("total_before", WaveTot, pixelSizeX , pixelSizeY);
    if (propagationDistance != 0.){PropagateWavefrontFresnel();}
    if(doTalbot){ApplyAbsorptionGrating(periodG2, disp, rotAxis, rotAngle);}
    //WriteToMHD("imag", waveImag, pixelSizeX , pixelSizeY);
    //WriteToMHD("real", waveReal, pixelSizeX , pixelSizeY);
    ComputeIntensity();
}


//----------------------------------------------------------------------
 
void GmGetWavefront::FillWavefront(const G4Step* aStep) {

    G4Track* track = aStep->GetTrack();
    PhotonPhaseInfo* phaseInfo = dynamic_cast<PhotonPhaseInfo*>(track->GetUserInformation());

    if (!phaseInfo) {
        G4cerr << "Error: Phase information not found in track" << G4endl;
        return;
    }
    // Detector Size
    double detectorSizeX = numPixelsX*pixelSizeX;  
    double detectorSizeY = numPixelsY*pixelSizeY; 

    // Interpolation
    double x = preR.x()+(postR.x()-preR.x())*(zstop-preR.z())/(postR.z()-preR.z());
    double y = preR.y()+(postR.y()-preR.y())*(zstop-preR.z())/(postR.z()-preR.z()); 

    // Obtain the coordinates of the pixel
    double X_moved = x + 0.5 * (detectorSizeX - pixelSizeX) - dispVP[0]; 
    int XCoord = std::floor(X_moved / pixelSizeX + 0.5); 

    double Y_moved = y + 0.5 * (detectorSizeY - pixelSizeY) - dispVP[1];
    int YCoord = std::floor(Y_moved / pixelSizeY  + 0.5); // Obtain the Y coordinate of the pixel

    if (XCoord < 0 || XCoord >= numPixelsX || YCoord < 0 || YCoord >= numPixelsY) return;

    
    double PreEnergy = aStep->GetPreStepPoint()->GetKineticEnergy(); // MeV
    //G4cout << PreEnergy/CLHEP::keV << G4endl; 
    //if (PreEnergy != 20.){
      //  G4cout << "Energy: " << PreEnergy << G4endl;
    //}
    energyAcummulated += PreEnergy;
    particleCount += 1;
    
    double phaseAccumulated = phaseInfo->GetPhase();
    //G4cout << phaseAccumulated << G4endl;

    // Wavefront construction
    waveReal[YCoord][XCoord] += std::sqrt(PreEnergy*1000)*std::cos(phaseAccumulated);
    waveImag[YCoord][XCoord] -= std::sqrt(PreEnergy*1000)*std::sin(phaseAccumulated);
    //double A = std::sqrt(PreEnergy*1000);
    //std::complex<double> inc = std::polar(A, -phaseAccumulated);

    //waveReal[YCoord][XCoord] -= inc.real();
    //waveImag[YCoord][XCoord] -= inc.imag();

    //sumPhase[YCoord][XCoord] += phaseAccumulated;
    //sumAmp[YCoord][XCoord] += A;
    //WaveTot[YCoord][XCoord] += PreEnergy*1000*(std::cos(phaseAccumulated)*std::cos(phaseAccumulated)+std::sin(phaseAccumulated)*std::sin(phaseAccumulated));
   
    // TEST
    //waveReal[YCoord][XCoord] += std::sqrt(PreEnergy);
    //waveImag[YCoord][XCoord] += 0;
    
    if (phaseVec[YCoord][XCoord] == 0.0 ) {
    phaseVec[YCoord][XCoord] = phaseAccumulated;}          
}

//----------------------------------------------------------------------
void GmGetWavefront::WriteWavefront(){
    std::ofstream file_real(outputFolder+outputFilename+"_real.out");
    std::ofstream file_imag(outputFolder+outputFilename+"_imag.out");

    for (int y = 0; y < numPixelsY; y++){
        for (int x = 0; x < numPixelsX; x++){
            file_real << waveReal[y][x] << " ";    
            file_imag << waveImag[y][x] << " ";   
        }
        file_real << "\n";  
        file_imag << "\n"; 
    }
    file_real.close();
    file_imag.close();
}

//----------------------------------------------------------------------
void GmGetWavefront::PropagateWavefrontFresnel(){

    
    //G4cout << "Total Energy: " << energyAcummulated << G4endl;
    using namespace std::chrono;
    auto start_total = high_resolution_clock::now();
    // TODO: Test which version is more efficient in time and memory
    int pad = 0;
    int paddedX = numPixelsX + 2 * pad;
    int paddedY = numPixelsY + 2 * pad;
    int N = paddedX * paddedY;

    fftw_plan plan_forward, plan_inverse;

    effectiveEnergy =  energyAcummulated / particleCount;
    G4cout << "Effective Energy: " << effectiveEnergy/CLHEP::keV << " keV" << G4endl;
    double lambda = 1.23984193 / (effectiveEnergy/CLHEP::eV) *CLHEP::um; // TODO: Energy must be obtained through the simulation or something; 
    double lambda_m = lambda / CLHEP::m;
    double dz = propagationDistance / magnification;
    double dz_m = dz/CLHEP::m;
    
    double fsX = 1.0 / (pixelSizeX/magnification * paddedX);
    double fsY = 1.0 / (pixelSizeY/magnification * paddedY);
    double fx_m = fsX * 1000;
    double fy_m = fsY * 1000.0;


    G4cout << "Wavelength: " << lambda << " mm"<< G4endl;
    G4cout << "Prop Distance: " << dz << " mm" << G4endl;
    G4cout << "fx: " << fsX << " 1/mm" << G4endl;
    G4cout << "fy: " << fsY << " 1/mm" << G4endl;

    auto start_fill = high_resolution_clock::now();
    
    //ApplyGaussianBlur(waveReal, waveImag, numPixelsX, numPixelsY);
    std::vector<std::complex<double>> wavefront = PadWavefrontReplication(waveReal, waveImag, pad, pad);


    auto end_fill = high_resolution_clock::now();
    G4cout << "Time filling in array : " << duration_cast<milliseconds>(end_fill-start_fill).count() << " ms" << G4endl; 
  
    // FFT plans
    auto start_plan = high_resolution_clock::now();
        //plan_forward = fftw_plan_dft_2d(numPixelsY, numPixelsX, reinterpret_cast<fftw_complex*>(in.data()), reinterpret_cast<fftw_complex*>(out.data()), FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward = fftw_plan_dft_2d(paddedY, paddedX, reinterpret_cast<fftw_complex*>(wavefront.data()), reinterpret_cast<fftw_complex*>(wavefront.data()), FFTW_FORWARD, FFTW_ESTIMATE);

        //plan_inverse = fftw_plan_dft_2d(numPixelsY, numPixelsX, reinterpret_cast<fftw_complex*>(out.data()), reinterpret_cast<fftw_complex*>(in.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
    plan_inverse = fftw_plan_dft_2d(paddedY, paddedX, reinterpret_cast<fftw_complex*>(wavefront.data()), reinterpret_cast<fftw_complex*>(wavefront.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
    
    auto end_plan = high_resolution_clock::now();
    G4cout << "Time creating FFT plans: " << duration_cast<milliseconds>(end_plan-start_plan).count() << " ms" << G4endl; 

    // Execute forward FFT
    auto start_forward = high_resolution_clock::now();
    fftw_execute(plan_forward);
    auto end_forward = high_resolution_clock::now();
    G4cout << "Time FFT forward: " << duration_cast<milliseconds>(end_forward-start_forward).count() << " ms" << G4endl; 
  
    auto start_prop = high_resolution_clock::now();

    for (int y = 0; y < paddedY; y++){
        // Spatial Frequencies
        double fy = (y < paddedY / 2) ? y*fsY : (y-paddedY)*fsY;
        for (int x = 0; x < paddedX; x++){
            int index = y * paddedX + x;

            // Spatial Frequencies
            double fx = (x < paddedX / 2) ? x * fsX : (x-paddedX) * fsX;
            // Fresnel Propagator
            
            //std::complex<double> propFactor = std::exp(-I * PI * lambda_m * dz_m * (fx_m*fx_m + fy_m*fy_m));
            std::complex<double> propFactor = std::exp(-I * PI * lambda/CLHEP::mm * dz/CLHEP::mm * (fx*fx + fy*fy));
            wavefront[index] *= propFactor;                
        }
    }
   
    
    auto end_prop = high_resolution_clock::now();
    G4cout << "Time propagation: " << duration_cast<milliseconds>(end_prop-start_prop).count() << " ms" << G4endl; 
  
    auto start_inverse = high_resolution_clock::now();
    fftw_execute(plan_inverse);
    auto end_inverse = high_resolution_clock::now();
    G4cout << "Time FFT inverse: " << duration_cast<milliseconds>(end_inverse-start_inverse).count() << " ms" << G4endl; 
  
    auto start_write = high_resolution_clock::now();
    for (int y = 0; y < numPixelsY; y++){
        for (int x = 0; x < numPixelsX; x++){
            int src_x = x + pad;
            int src_y = y + pad;
            int index = src_y * paddedX + src_x;
            //int index = y * numPixelsX + x;
            //double sign = ((x + y) % 2 == 0) ? 1.0 : -1.0;
            //intensity[y][x] =  std::norm(wavefront[index]) / N; // Problems if doTalbot = true
            //wavefront[index] *= std::exp(I * 2.0 * PI * lambda/CLHEP::mm * dz/CLHEP::mm);    
            waveReal[y][x] = wavefront[index].real()  /N;   
            waveImag[y][x] = wavefront[index].imag()  /N;     
        }
    }


    auto end_total = high_resolution_clock::now();
    G4cout << "Time Write wavefront: " << duration_cast<milliseconds>(end_total-start_write).count() << " ms" << G4endl;
    
    G4cout << "Time Total prop: " << duration_cast<milliseconds>(end_total-start_total).count() << " ms" << G4endl; 
  
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_inverse);
    fftw_cleanup();
}

//----------------------------------------------------------------------
void GmGetWavefront::CheckNyquist(double px, double py, double wavelength, double dz){

    double maxPixelSize =  std::sqrt(wavelength*dz);
    if (px > maxPixelSize || py > maxPixelSize){
        G4cerr << "Propagator is not sampled above the Nyquist frequency" << G4endl;
        G4cerr << "Grid Size: (" << px << ", " << py << ")" << G4endl; 
        G4cerr << "Max. Grid Size: (" << maxPixelSize<< ", " << maxPixelSize << ")" << G4endl;  
    }
}
//----------------------------------------------------------------------
void GmGetWavefront::ApplyAbsorptionGrating(double period, const G4ThreeVector& displacement, const G4ThreeVector& axis, double theta) {

    //G4ThreeVector gratingDir = axis.orthogonal().unit();
    G4ThreeVector gratingDir = G4ThreeVector(1,0,0);

    double period_px = period / pixelSizeX * magnification;
    G4cout << "period: " << period  << " um, " << period_px << " pixels" << G4endl;
    double displacement_px = displacement.x() / pixelSizeX;

    for (int y = 0; y < numPixelsY; y++){ 
        for (int x = 0; x < numPixelsX; x++){
            
            //G4ThreeVector pos(x*pixelSizeX/magnification, y*pixelSizeY/magnification, 0);
            //double xPos = (x - numPixelsX/2.0) * pixelSizeX / magnification;
            //double yPos = (y - numPixelsY/2.0) * pixelSizeY / magnification;
            //G4cout << xPos << G4endl;
            //G4ThreeVector pos(xPos, yPos, 0);
            //G4cout << "displacement: " << displacement << G4endl;
            //pos -=displacement;
            //G4cout << "pos before rotation: " << pos << G4endl;

            //G4ThreeVector rotatedPoint = Rotate(pos, axis, theta);
            
            
            //G4cout << "pos after rotation: " << rotatedPoint << G4endl;
            //double proj = rotatedPoint.dot(gratingDir);

            //double proj = pos.dot(gratingDir);

            //int newX = static_cast<int>(rotatedPoint.x() / pixelSizeX);
            //int newY = static_cast<int>(rotatedPoint.y() / pixelSizeY);

            //double modProj = std::fmod(proj, period);
            //if (modProj < 0) modProj += period;

            //double modPositionX = std::fmod(rotatedPoint.x(), period);

            //if (modPositionX < 0){ modPositionX += period; }
            //if (modPositionY < 0){ modPositionY += period; }
            //int stripeIndex = static_cast<int>(std::floor(proj / (period / 2.0)));
            //double transmission = (stripeIndex % 2 == 0) ? 1.0 : 0.0;
            //double transmission = (modProj < period/2 ) ? 1.0 : 0.0;

            //double modProj = std::fmod(proj, period);
            double modProj = fmod(x - displacement_px, period_px);
            if (modProj < 0) modProj += period_px;

            double transmission = (modProj < period_px/2 ) ? 1.0 : 0.0;

            waveReal[y][x] *= transmission;
            waveImag[y][x] *= transmission; 
        }   
    }
}


G4ThreeVector GmGetWavefront::Rotate(const G4ThreeVector& point, const G4ThreeVector& axis, double theta) {

    G4ThreeVector v = axis.unit();
    double cosT = std::cos(theta/CLHEP::rad);
    double sinT = std::sin(theta/CLHEP::rad);
    return point*cosT + (v.cross(point)) * sinT + v*(v.dot(point)) * (1-cosT);
}

// NOT WORKING PROPERLY
/*/void GmGetWavefront::PropagateWavefrontFresnelByBlocks(){
    int blockSizeX = 128; 
    int blockSizeY = 128; 
    int overlap = blockSizeX / 4;
    int stepSizeX = blockSizeX - overlap;
    int stepSizeY = blockSizeY - overlap;
    int paddedSizeX = 2 * blockSizeX;
    int paddedSizeY = 2 * blockSizeY;
    int N = paddedSizeX * paddedSizeY;

    fftw_complex *in, *out, *out2;
    fftw_plan plan_forward, plan_inverse;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    plan_forward = fftw_plan_dft_2d(paddedSizeY, paddedSizeX, in, out, FFTW_FORWARD, FFTW_MEASURE);
    plan_inverse = fftw_plan_dft_2d(paddedSizeY, paddedSizeX, out, out2, FFTW_BACKWARD, FFTW_MEASURE);


    double lambda = 1.23984193 / (20 * 1000) / 1000 * CLHEP::mm;
    double dz = propagationDistance;
    double k = 2 * PI / lambda;
    double fsX = 1.0 / (pixelSizeX * paddedSizeX);
    double fsY = 1.0 / (pixelSizeY * paddedSizeY);

    std::vector<std::vector<double>> hannWindow(blockSizeY, std::vector<double>(blockSizeX));
    for (int y = 0; y < blockSizeY; y++) {
        for (int x = 0; x < blockSizeX; x++) {
            hannWindow[y][x] = 0.5*(1 - std::cos(2 * PI * x / (blockSizeX - 1))) * 0.5*(1 - std::cos(2 * PI * y / (blockSizeY - 1)));
        }
    }

    for (int by = 0; by < numPixelsY; by += stepSizeY) {
        for (int bx = 0; bx < numPixelsX; bx += stepSizeX) {
            memset(in, 0, sizeof(fftw_complex) * paddedSizeX * paddedSizeY);
            
            for (int y = 0; y < blockSizeY; y++) {
                for (int x = 0; x < blockSizeX; x++) {
                   
                    if ((by + y) < numPixelsY && (bx + x) < numPixelsX) {
                        in[y * paddedSizeX + x][0] = waveReal[by + y][bx + x]*hannWindow[y][x];
                        in[y * paddedSizeX + x][1] = -waveImag[by + y][bx + x]*hannWindow[y][x];
                    } 
                }
            }
           
            fftw_execute(plan_forward);
            for (int y = 0; y < paddedSizeY; y++) {
                double fy = (y < paddedSizeY / 2) ? y * fsY : (y - paddedSizeY) * fsY;
                for (int x = 0; x < paddedSizeX; x++) {
                   
                    double fx = (x < paddedSizeX / 2) ? x * fsX : (x - paddedSizeX) * fsX;
                    std::complex<double> propFactor = std::exp(-I * PI * lambda * dz * (fx * fx + fy * fy));
                    std::complex<double> outComplex(out[y * paddedSizeX + x][0], out[y * paddedSizeX + x][1]);
                    outComplex *= propFactor;
                    out[y * paddedSizeX + x][0] = outComplex.real();
                    out[y * paddedSizeX + x][1] = outComplex.imag();
                }
            }

            fftw_execute(plan_inverse);

            for (int y =0; y < blockSizeY; y++) {
                for (int x =0; x < blockSizeX; x++) {
                    
                    if ((by + y) < numPixelsY && (bx + x) < numPixelsX) {
                        waveReal[by + y][bx + x] += out2[y * paddedSizeX + x][0] *hannWindow[y][x]/ (paddedSizeX * paddedSizeY);
                        waveImag[by + y][bx + x] += -out2[y * paddedSizeX + x][1]*hannWindow[y][x]/ (paddedSizeX * paddedSizeY);
                    }
                }
            }

            
            
        }
    }
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_inverse);
    fftw_free(in);
    fftw_free(out);
    fftw_free(out2);
}
/*/
//----------------------------------------------------------------------
void GmGetWavefront::ComputeIntensity(){
    using namespace std::chrono;
    auto start_total = high_resolution_clock::now();
    //double totalEnergy = 0.0;
    for (int y = 0; y < numPixelsY; y++){
        for (int x = 0; x < numPixelsX; x++){
            intensity[y][x] = (waveReal[y][x]*waveReal[y][x] + waveImag[y][x]*waveImag[y][x]) / (magnification*magnification);
            //totalEnergy += intensity[y][x];
        }
    }
    //G4cout << "Total Energy: " << totalEnergy << G4endl;
    DownsampleIntensity();
    
    //downsampledIntensity = downsample(intensity, detnumPixelsY, detnumPixelsX);
    double pSizeX = static_cast<double>(numPixelsX) / static_cast<double>(detnumPixelsX) * pixelSizeX/CLHEP::mm / magnification;
    double pSizeY = static_cast<double>(numPixelsY) / static_cast<double>(detnumPixelsY) * pixelSizeY/CLHEP::mm / magnification;
    //G4cout << pSizeX << G4endl;
    WriteToMHD("intensity", downsampledIntensity, pSizeX , pSizeY);
    //WriteToMHD("intensity", intensity, pixelSizeX/CLHEP::mm / magnification ,pixelSizeX/CLHEP::mm / magnification);    
    auto end_total = high_resolution_clock::now();
    G4cout << "Time write Intensity and MHD: " << duration_cast<milliseconds>(end_total-start_total).count() << " ms" << G4endl;  
}
//----------------------------------------------------------------------
void GmGetWavefront::DownsampleIntensity() {

    int blockSizeX = numPixelsX / detnumPixelsX;
    int blockSizeY = numPixelsY / detnumPixelsY;

    //double normalization = 1.0 / (blockSizeX * blockSizeY);

    if (blockSizeX == 1. && blockSizeY == 1.){
        for (int y = 0; y < numPixelsY; y++) {
            for (int x = 0; x < numPixelsX; x++) {
                downsampledIntensity[y][x] = intensity[y][x];
            }
        }    
    } 

    //G4cout << blockSizeX << G4endl;
    else {
        for (int iy = 0; iy < detnumPixelsY; iy++){ 
            for (int ix = 0; ix < detnumPixelsX; ix++){
                double sum = 0.0;
                int count = 0;

                for (int dy = 0; dy < blockSizeY; dy++){ 
                    for (int dx = 0; dx < blockSizeX; dx++){
                        int originalX = ix* blockSizeX + dx;
                        int originalY = iy* blockSizeY + dy;

                        if (originalX < numPixelsX && originalY < numPixelsY){
                            sum += intensity[originalY][originalX];
                            count++;
                        }
                    }
                }
                downsampledIntensity[iy][ix] =  sum/count;
            }
        }
    }
}
//----------------------------------------------------------------------
void GmGetWavefront::WriteToMHD(const std::string& name, std::vector<std::vector<double>> data, double pSizeX, double pSizeY) {
    int rows = data.size(); // Y, vertical
    int cols = data[0].size(); // X, horizontal

    std::string rawFilename = outputFilename +"_" + name +".raw";
    std::string mhdFilename = outputFilename +"_" + name +".mhd";

    std::ofstream rawFile(outputFolder+rawFilename, std::ios::binary);
    if (!rawFile) {
        G4cerr << "Error: File " << rawFilename << "can't be opened" << G4endl;
        return;}

    for (int y = 0; y < rows; ++y) {
        std::vector<float> row(cols);
        for (int x = 0; x < cols; ++x) {
            row[x] = data[y][x];
            
            }
        rawFile.write(reinterpret_cast<const char*>(row.data()),cols * sizeof(float));
        }
    rawFile.close();

    std::ofstream mhdFile(outputFolder+mhdFilename);
    if (!mhdFile){
        G4cerr << "Error: File " << mhdFilename << "can't be opened" << G4endl;
        return;
    }
    
    mhdFile << "ObjectType = Image\n";
    mhdFile << "BinaryData = True\n";
    mhdFile << "BinaryDataByteOrderMSB = False\n";
    mhdFile << "CompressedData = False\n";
    mhdFile << "NDims = 2\n";
    mhdFile << "DimSize = " << cols << " " << rows << "\n";
    mhdFile << "ElementSpacing = " << pSizeX << " " << pSizeY << "\n";
    mhdFile << "ElementType = MET_FLOAT\n";
    mhdFile << "DSO = " << zstop <<" mm" << "\n";
    mhdFile << "DSD = " << zstop+propagationDistance <<" mm" << "\n";
    mhdFile << "GeomMagnification = " << magnification << "\n";
    mhdFile << "SourceDirectionDistribution = " << distDirectionName << "\n";
    mhdFile << "EffectiveEnergy = " << effectiveEnergy/CLHEP::eV << " eV" << "\n";
    mhdFile << "EventsSimulated = " << eventsSimulated << "\n";
    mhdFile << "ElementDataFile = " << rawFilename << "\n";

    mhdFile.close();
    
}


std::vector<std::complex<double>> GmGetWavefront::PadWavefrontReplication(
    const std::vector<std::vector<double>>& Real,
    const std::vector<std::vector<double>>& Imag,
    int padX, int padY) {

    int origY = Real.size();
    int origX = Real[0].size();
    
    int newY = origY + 2 * padY;
    int newX = origX + 2 * padX;

    std::vector<std::complex<double>> padded(newY * newX);

    // DEBUG
    std::vector<std::vector<double>> wavefront_imag;
    wavefront_imag.resize(newY, std::vector<double>(newX, 0.0));
    std::vector<std::vector<double>> wavefront_real;
    wavefront_real.resize(newY, std::vector<double>(newX, 0.0));

    auto mirror = [](int i, int len){
        if (i < 0) return -i - 1;
        if (i >= len) return 2 * len - i -1;
        return i;
        
    };

    for (int y = 0; y < newY; y++) {
        int srcY = mirror(y-padY, origY);
        for (int x = 0; x < newX; x++) {
            int srcX = mirror(x - padX, origX);
            double realVal = Real[srcY][srcX];
            double imagVal = Imag[srcY][srcX];
            padded[y * newX + x] = std::complex<double>(realVal, -imagVal);
            // DEBUG!!
            //wavefront_real[y][x] = realVal;
            //wavefront_imag[y][x] = imagVal;
        }
    }

    // DEBUG!!
    //WriteToMHD("real_padded", wavefront_real, pixelSizeX , pixelSizeY);
    //WriteToMHD("imag_padded", wavefront_imag, pixelSizeX , pixelSizeY);

    return padded;
}


