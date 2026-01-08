#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"
#include "G4EmCalculator.hh"
#include "G4Gamma.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>


#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmimgle/dcmimage.h"
#include "dcmtk/dcmdata/dcfilefo.h"
#include "dcmtk/dcmdata/dcdeftag.h"
#include "GmGetImage.hh"
#include "GamosCore/GamosBase/Base/include/GmParameterMgr.hh"
#include "G4GeometryTolerance.hh"

constexpr double MIN_U2 = 1.0e-8;
constexpr double PI = CLHEP::pi;



//---------------------------------------------------------------------
GmGetImage::GmGetImage() : GmUserRunAction(), GmUserEventAction(), GmUserSteppingAction()
{
}


//----------------------------------------------------------------------
void GmGetImage::BeginOfRunAction(const G4Run* run) {

    GetParams();
   
    InitializeVector(gridInterspaceMu,"InterStripMu.txt");
    InitializeVector(gridStripMu,"StripMu.txt");

    totalEvents = run->GetNumberOfEventToBeProcessed();

    if (detectorModel == "VD") {
        DetectorMaterial = G4NistManager::Instance()->FindOrBuildMaterial(detectorMaterial);
        if (!DetectorMaterial) {
            DetectorMaterial = LoadCustomMaterial(detectorMaterial, "geom/MyMaterials.txt");
            //G4cout << detectorMaterial << " Material Added" << G4endl;
        }
        //InitializeVector(efficiencyMap,"efficiency.txt");
    }
    
    G4cout << "=======================" << G4endl;
    G4cout << "=  Simulation Param.  =" << G4endl;
    G4cout << "=======================" << G4endl;
    G4cout << "   Output selection: " << outputSel << G4endl;
    G4cout << "   Output Format: " << outputFormat << G4endl;
    G4cout << "   Output Folder: " << outputFolder << G4endl;
    G4cout << "   Output Filename: " << outputFilename << G4endl;
    G4cout << "   Total Events: " << totalEvents << G4endl;

    G4cout << "=======================" << G4endl;
    G4cout << "= " << detectorModel << " Detector ON = " << G4endl;
    G4cout << "=======================" << G4endl;
   
    G4cout << "   Z Stop Position: " << zstop << " mm " << G4endl;
    G4cout << "   Number of Pixels on X: " << numPixelsX << G4endl;
    G4cout << "   Number of Pixels on Y: " << numPixelsY << G4endl;
    G4cout << "   Pixel Size on X: " << pixelSizeX << " mm" << G4endl;
    G4cout << "   Pixel Size on Y: " << pixelSizeY << " mm" << G4endl;
    G4cout << "   Detector Thickness: " << detectorThickness << " mm" << G4endl;
    G4cout << "   Detector Material: " << detectorMaterial << G4endl;
   
    G4cout << "==========================" << G4endl;
    G4cout << " " << G4endl;
    G4cout << "==========================" << G4endl;
    G4cout << "=====  Grid Params  ======" << G4endl;
    G4cout << "==========================" << G4endl;
    G4cout << "   Use Scatter-Grid: "             << apply_grid << G4endl;
    G4cout << "   Grid Ratio: " << gridRatio << G4endl; // What Grid Ratio is?
    G4cout << "   Source-Detector Distance (SDD): " << sdd << "mm" << G4endl;
    G4cout << "   Grid Frequency: " << gridFreq << " per unit" << G4endl;
    G4cout << "   Grid Strip Thickness: " << gridStripThickness << "mm" << G4endl;
    G4cout << "   Gap grid displacement: " << gap << "mm" << G4endl;
    G4cout << "==========================" << G4endl;

    /*/
    const G4MaterialTable* matTable = G4Material::GetMaterialTable();
    for (const auto& material : *matTable) {
        G4cout << "Material: " << material->GetName() << G4endl;    
    }
    /*/
}

//----------------------------------------------------------------------
void GmGetImage::GetParams(){

    GmParameterMgr* paramMgr = GmParameterMgr::GetInstance();

    // Get detector parameters from GetImage.in file
    detectorModel = paramMgr->GetStringValue("GetImage:DetectorModel", "VD");
    outputSel = paramMgr->GetStringValue("GetImage:OutputType", "Energy");   
    numPixelsX = paramMgr->GetNumericValue("GetImage:NumPixelsX", 1000);
    numPixelsY = paramMgr->GetNumericValue("GetImage:NumPixelsY", 1000);
    pixelSizeX = paramMgr->GetNumericValue("GetImage:PixelSizeX", 0.1);
    pixelSizeY = paramMgr->GetNumericValue("GetImage:PixelSizeY", 0.1);
    detectorThickness = paramMgr->GetNumericValue("GetImage:DetectorThickness", 0.2);
    detectorMaterial = paramMgr->GetStringValue("GetImage:DetectorMaterial", "aSe");
    zstop = paramMgr->GetNumericValue("GetImage:zStop",650);

    // Adjust energy matrix size from deterctor size
    output.resize(numPixelsY, std::vector<double>(numPixelsX, 0.0));

    // Get grid parameters from WriteEnergy.in file
    apply_grid = paramMgr->GetStringValue("Grid:UseGrid", "Yes");
    gridRatio = paramMgr->GetNumericValue("Grid:Ratio", 1.0);
    sdd = paramMgr->GetNumericValue("Grid:SourceDetectorDistance", 700.0); // Is it really necessary¿?
    gridFreq = paramMgr->GetNumericValue("Grid:Frequency", 1.0);
    gridStripThickness = paramMgr->GetNumericValue("Grid:StripThickness", 0.1);
    gap = paramMgr-> GetNumericValue("Grid:Gap", 0.0);

    // Charge-Energy parameters
    pairCreationEnergyEV = paramMgr->GetNumericValue("GetImage:PairCreationEnergy", 50.0); // eV

    // Sensibility Curve Parameters Initialization
    SensibilityOffset = 0.0;
    SensibilityFactor = 1.0;

    // Output File Parameters
    outputFolder = paramMgr->GetStringValue("GetImage:ResultsFolder", "output/");	
	outputFilename = paramMgr->GetStringValue("GetImage:OutputFilename", "energy");
    outputFormat = paramMgr->GetStringValue("GetImage:OutputFormat", "MHD/RAW");

    // Exceptions

    if( (detectorModel != "VD" && detectorModel != "MCD")) {
    G4Exception("GetImage:DetectorModel",
		"Wrong argument",
		FatalErrorInArgument,
		"Detector Model is not well defined, please use VD or MCD");
    }

    if( (outputSel != "Energy" && outputSel != "Charge" && outputSel != "Gray")) {
    G4Exception("GetImage:OutputSelection",
		"Wrong argument",
		FatalErrorInArgument,
		"Output Selection is not well defined, please use Energy, Charge");
    }

    if( (outputFormat != "MHD/RAW" && outputFormat != "Text" && outputFormat != "DCM")) {
    G4Exception("GetImage:OutputFormat",
		"Wrong argument",
		FatalErrorInArgument,
		"Output Format is not well defined, please use MHD/RAW, DCM or Text");
    }

    
}

//----------------------------------------------------------------------
void GmGetImage::UserSteppingAction(const G4Step* aStep)
{

  //G4cout << "Step " << G4endl;
  // Get Positions at the PreStep and PostStep Points
  preR= aStep->GetPreStepPoint()->GetPosition();
  postR= aStep->GetPostStepPoint()->GetPosition();
  
  //G4cout << "pre Z: " << preR.z() << G4endl;  
  
  // If PreStep position is larger than zStop, the particle has crossed the detector plane, the Track is killed. Valid for both detector's model.
  if ( preR.z() > zstop) {
      //G4cout << "pre Z: " << preR.z() << G4endl;
  	  aStep->GetTrack()->SetTrackStatus(fStopAndKill);
  }

  
  // If PreStep position is lower than zStop and PostStep position is larger, the particle is crossing the detector plane.
  if ( (preR.z()-zstop) * (postR.z()-zstop) < 0. || fabs(postR.z()-zstop) < G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() ) {
    //G4cout << "Step Accepted: " << preR.z() << G4endl;
    
    //Apply antiscatter-grid
    if (detectorModel == "VD"){
      if (apply_grid == "Yes"){
          if (AcceptStep(aStep)){FillMatrix_Virtual(aStep);}
      }else{FillMatrix_Virtual(aStep);}  
    }
    if (detectorModel == "MCD"){
      if (apply_grid == "Yes"){
          if (!AcceptStep(aStep)){aStep->GetTrack()->SetTrackStatus(fStopAndKill);}
      }  
    }
  }
}

//----------------------------------------------------------------------
void GmGetImage::EndOfEventAction( const G4Event* ) 
{
  // Just for MCD detector Model
  if (detectorModel=="MCD"){  
      std::vector<GmRecHit*> rhits = GmHitsEventMgr::GetInstance()->GetAllRecHits();

      for (auto& hit : rhits) {
          FillMatrix_MCD(hit);
      }
  }
}

//----------------------------------------------------------------------
void GmGetImage::EndOfRunAction(const G4Run*) {
    // Valid for both detector's Models
    //WriteEnergyMatrixToFile(filename);
    if ( outputFormat == "MHD/RAW"){WriteMHD(outputFilename);}
    if ( outputFormat == "Text"){WriteText(outputFilename);}
    if ( outputFormat == "DCM"){WriteDCM(outputFilename);}
}



//----------------------------------------------------------------------
void GmGetImage::InitializeVector(std::vector<std::pair<double, double>>& vec, const std::string& file_name) {
    
    // Get the path to the directory with the files
    std::string directory_path_str = "plug-ins/resources/";
    std::filesystem::path full_path = std::filesystem::current_path() / directory_path_str / file_name;

    // Read the file and store the values in the map

    std::ifstream infile(full_path);
    if (!infile.is_open()) {
        G4cerr << "Couldn't read file: " << full_path << G4endl;
        return;
    }

    double Energy, value;
    vec.clear();

   
    while((infile >> Energy >> value)) {
        vec.emplace_back(Energy, value);
        continue;
    }
    
    std::sort(vec.begin(), vec.end());
    G4cout << "=======================" << G4endl;
    G4cout << file_name << " took from: " << full_path << G4endl;
    G4cout << "=======================" << G4endl; 
}

double GmGetImage::GetInterpolatedValue(const std::vector<std::pair<double, double>>& vec, double stepEnergy ){

     // Check if the map is empty to avoid execution errors.
    if (vec.empty()) {return stepEnergy;}

    // Find first element larger than stepEnergy.
    auto it = std::lower_bound(vec.begin(), vec.end(), stepEnergy, [](const std::pair<double, double>& a, double value) {return a.first < value;});
    
    // If it is the beginning, the first value is used.
    if (it == vec.begin()) {return stepEnergy * it->second;}
    
    // If it is the end, the last value is used.
    if (it == vec.end()) {return stepEnergy * std::prev(it)->second;}
    

    // Now, it points to the first element whose key is higher than stepEnergy.
    // Check if the previous element is closer to stepEnergy.
    auto prev_it = std::prev(it);

    double x1 = prev_it->first, y1 = prev_it->second;
    double x2 = it->first, y2 = it->second;
    double interpolated = y1 + (y2 - y1)*(stepEnergy - x1) / (x2 - x1);

    return interpolated;
}

//----------------------------------------------------------------------
G4bool GmGetImage::AcceptStep(const G4Step* aStep) {

    G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition(); // Get the position of the Step
    G4ThreeVector direction = aStep->GetPostStepPoint()->GetMomentumDirection(); // Get the momentum (direction) of the Step

    // Compute transmission probability
    G4double transmissionProbability = CalculateGridTransmissionProbability(position, direction, aStep);

    // Decide if accept or reject the Step based on the transmission probability
	if (G4UniformRand() < 1 - transmissionProbability) return false; // Reject Step
    return true;
}

//----------------------------------------------------------------------
G4float GmGetImage::CalculateGridTransmissionProbability(const G4ThreeVector& position, const G4ThreeVector& direction, const G4Step* aStep) {

    //const double PI = 3.14159265358979323846;
    // * gridRatio is the ratio between the height of the strips and the distance between the strips
    // * gridFreq is the inverse of the Width of a grid unit cell
    // * position is in mm
    // * prefix 2 means on the rectangular system of reference (d2, D2, u2)
    
    // This function computes the Transmission Probability based on the work of Day and Dance [Phys Med Biol 28, p. 1429-1433 (1983)]
    // The origin is at the center on the detector. 
    // The x-ray source is shifted, so it is mandatory to correct the x position on the calculus of grid_angle.
    

    G4double grid_angle, u, w, gridStripMu_E, gridInterspaceMu_E;

    //Look u coeficients
    double PreEnergy = aStep->GetPreStepPoint()->GetKineticEnergy() *1000; // The energy at PreStep point in keV

    gridStripMu_E = GetInterpolatedValue(gridStripMu,PreEnergy)/ 10; // Change to mm⁻1
    gridInterspaceMu_E = GetInterpolatedValue(gridInterspaceMu,PreEnergy)/ 10; // Change to mm⁻1

    //G4cout << "energy: " << aStep->GetPreStepPoint()->GetKineticEnergy()  << "Mu_Pb: " << gridStripMu_aux << G4endl;

    
    grid_angle = (0.5 * PI) - atan2(position.x()-gap,sdd); // Angle in the XZ oblique coordinate system. If grid_angle = pi/2 the grid is not oblique.
    //grid_angle = PI / 2;
    // director cosines
    u = direction.x(); // The grid is periodic along X-axis
    w = direction.z(); // The particles travel along Z-axis
 
    // based on Day&Dance1983
    
    G4float C = 1.0f / gridFreq; // Width of a grid "unit cell" 
    G4float d2 = gridStripThickness / std::sin(grid_angle); // thickness of the strip region (mm) 
    G4float D2 = C - d2; // thickness of the interspace region (mm)
    G4float h = std::fabs(gridRatio) * D2; // Height of the grid (mm)

	// Transformation between coordinates in oblique and rectangular coodinate systems
    G4double u2 = std::fabs(u - w / std::tan(grid_angle)); // Eq 1 of Day&Dance1983
    if (u2 < 1.0e-9) {
        u2 = 1.0e-8;
    } // To avoid NaN

    G4double P = (h / w) * u2; // Absolute value of the oblique projection of the photon path (eq. 4 Day&Dance1983)
    G4double n = std::floor(P * gridFreq); // Number of complete "unit cells" penetrated by the photon path (eq. 7 Day&Dance1983)
    G4float q = P - n * C; // The rest ¿?
    G4double alpha = u2 / (gridStripMu_E - gridInterspaceMu_E); // (eq. 9 Day&Dance1983)
    G4double inv_alpha = 1.0 / alpha;

    G4float A = std::exp(-gridInterspaceMu_E * h / w - d2 * n * inv_alpha); // Transmission is T = A*exp(-q/alpha)
    G4float H = 0.0f;

    if (q >= D2) {
        H = 1.0f; // Step function H(q-D2)
    }

    G4float B = (std::fabs(q - D2) + 2.0f * alpha) * std::exp((H * (D2 - q)) * inv_alpha) + (std::fabs(d2 - q) - 2.0f * alpha) * std::exp((-0.5f * (d2 + q - std::fabs(d2 - q))) * inv_alpha); // (eq. 12 Day&Dance1983) 
    
    return (A * B * gridFreq); // Day&Dance1983

    //return 0.2;
}

//----------------------------------------------------------------------
void GmGetImage::FillMatrix_Virtual(const G4Step* aStep) {

    //G4cout << "Filling VD Matrix" << G4endl;
    // Detector Size
    double detectorSizeX = numPixelsX*pixelSizeX;  
    double detectorSizeY = numPixelsY*pixelSizeY; 

    // Interpolation
    double x = preR.x()+(postR.x()-preR.x())*(zstop-preR.z())/(postR.z()-preR.z());
    double y = preR.y()+(postR.y()-preR.y())*(zstop-preR.z())/(postR.z()-preR.z()); 

    // Obtain the coordinates of the pixel
    double X_moved = x + 0.5 * (detectorSizeX - pixelSizeX); 
    int XCoord = static_cast<int>(X_moved / pixelSizeX); 

    double Y_moved = y + 0.5 * (detectorSizeY - pixelSizeY);
    int YCoord = static_cast<int>(Y_moved / pixelSizeY); // Obtain the Y coordinate of the pixel

    if (XCoord < 0 || XCoord >= numPixelsX || YCoord < 0 || YCoord >= numPixelsY) return;

    // The energy at PreStep point in keV
    G4double PreEnergy = aStep->GetPreStepPoint()->GetKineticEnergy(); 

    // Efficiency Factor
    static G4EmCalculator emCal;
    //G4cout << DetectorMaterial << G4endl;
    G4double lambda_phot = emCal.ComputeMeanFreePath(PreEnergy, G4Gamma::Gamma(), "phot", DetectorMaterial);
    G4double lambda_compt = emCal.ComputeMeanFreePath(PreEnergy, G4Gamma::Gamma(), "compt", DetectorMaterial);
    G4double lambda_conv = emCal.ComputeMeanFreePath(PreEnergy, G4Gamma::Gamma(), "conv", DetectorMaterial);
    G4double lambda_rayl = emCal.ComputeMeanFreePath(PreEnergy, G4Gamma::Gamma(), "Rayl", DetectorMaterial);
    G4double lambda_tot = 1/(1/lambda_phot + 1/lambda_compt + 1/lambda_conv+ 1/lambda_rayl);

    /*/
    G4cout << "MFP Phot (cm): " << lambda_phot/cm << G4endl;
    G4cout << "MFP Compt (cm): " << lambda_compt/cm << G4endl;
    G4cout << "MFP Conv (cm): " << lambda_conv/cm << G4endl;
    G4cout << "MFP Rayl (cm): " << lambda_rayl/cm << G4endl;
    G4cout << "MFP Total (cm): " << lambda_tot/cm << G4endl;
    /*/

    G4double efficiencyFactor = lambda_tot/lambda_phot*(1 - std::exp(-detectorThickness/lambda_tot));
    //G4cout << "Efficiency factor for " << PreEnergy << " MeV and " << detectorThickness/mm << " mm: " << efficiencyFactor << G4endl; 
    //double efficiencyFactor = GetInterpolatedValue(efficiencyMap,PreEnergy); 

    if (G4UniformRand() > efficiencyFactor) return;

    // Energy deposited in eV
    double EnergyDeposited = PreEnergy / eV;

    
    if(outputSel == "Energy"){output[YCoord][XCoord] += EnergyDeposited;}
    else if(outputSel == "Charge"){output[YCoord][XCoord] += EnergyToCharge(EnergyDeposited);}
    else if(outputSel == "Gray"){
        output[YCoord][XCoord] += EnergyToCharge(EnergyDeposited);
        SensibilityOffset = 49.07;
        SensibilityFactor = 4.2e-04;
        }
    }       

//----------------------------------------------------------------------
void GmGetImage::FillMatrix_MCD(GmRecHit* hit) {
    double detectorSizeX = numPixelsX*pixelSizeX; // Define el tamaño del detector en X
    double detectorSizeY = numPixelsY*pixelSizeY; // Define el tamaño del detector en Y

    double X_moved = hit->GetPosition().x() + 0.5 * (detectorSizeX - pixelSizeX);
    int XCoord = static_cast<int>(X_moved / pixelSizeX);

    double Y_moved = hit->GetPosition().y() + 0.5 * (detectorSizeY - pixelSizeY);
    int YCoord = static_cast<int>(Y_moved / pixelSizeY);

    if (XCoord < 0 || XCoord >= numPixelsX) return;
    if (YCoord < 0 || YCoord >= numPixelsY) return;

    double hitEnergy = hit->GetEnergy() * 1000000;
    if(outputSel == "Energy"){output[YCoord][XCoord] += hitEnergy;}
    else if(outputSel == "Charge"){output[YCoord][XCoord] += EnergyToCharge(hitEnergy);}
    else if(outputSel == "Gray"){
        output[YCoord][XCoord] += EnergyToCharge(hitEnergy);
        SensibilityOffset = 48.62;
        SensibilityFactor = 5.9e-4;
        }
}

//----------------------------------------------------------------------
void GmGetImage::WriteText(const std::string& name) {
    std::ofstream file(outputFolder+name+".out");
    if (!file.is_open()) {
        G4cerr << "The file: " << name << " couldn't be found/opened." << G4endl;
        return;
    }

    for (int y = 0; y < numPixelsY; ++y) {
        for (int x = 0; x < numPixelsX; ++x) {
            file << SensibilityOffset + SensibilityFactor*output[y][x];
            if (x < numPixelsX - 1) {
                file << " "; // Split the output values with a space
            }
        }
        file << "\n"; // New line at the end of each row of the matrix
    }

    file.close();
    G4cout << "Text File written on: " << name << G4endl;
}

//----------------------------------------------------------------------
void GmGetImage::WriteMHD(const std::string& name) {
    int rows = numPixelsY;
    int cols = numPixelsX;

    std::string rawFilename = name +".raw";
    std::string mhdFilename = name +".mhd";

    std::ofstream rawFile(outputFolder+rawFilename, std::ios::binary);
    if (!rawFile) {
        G4cerr << "Error: File " << rawFilename << "can't be opened" << G4endl;
        return;}

    for (int y = 0; y < numPixelsY; ++y) {
        std::vector<float> row(numPixelsX);
        for (int x = 0; x < numPixelsX; ++x) {
            row[x] = SensibilityOffset + SensibilityFactor * output[y][x];
            
            }
        rawFile.write(reinterpret_cast<const char*>(row.data()),numPixelsX * sizeof(float));
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
    mhdFile << "ElementSpacing = " << pixelSizeX << " " << pixelSizeY << "\n";
    mhdFile << "ElementType = MET_FLOAT\n";
    mhdFile << "ElementDataFile = " << rawFilename << "\n";
    mhdFile.close();
    
}

//----------------------------------------------------------------------
void GmGetImage::WriteDCM(const std::string& name) {
    int rows = numPixelsY;
    int cols = numPixelsX;

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << pixelSizeX << "\\" << pixelSizeY;
    std::string ImagerPixelSpacing = oss.str();

    std::ostringstream oss2;
    oss2 << std::fixed << std::setprecision(0) << totalEvents;
    std::string totalEventsStr = oss2.str();
    
    //double pixelSpacing[] = {pixelSizeX, pixelSizeY}

    DcmFileFormat dcmfileformat;
    DcmDataset *dataset = dcmfileformat.getDataset();

    dataset->putAndInsertString(DCM_SOPClassUID, UID_SecondaryCaptureImageStorage);
    dataset->putAndInsertString(DCM_PatientName, "Phantom");
    dataset->putAndInsertString(DCM_PatientID, "123");
    dataset->putAndInsertString(DCM_ImagerPixelSpacing, ImagerPixelSpacing.c_str());
    dataset->putAndInsertUint16(DCM_Rows, rows);
    dataset->putAndInsertUint16(DCM_Columns, cols);
    dataset->putAndInsertUint16(DCM_BitsAllocated, 16);
    dataset->putAndInsertUint16(DCM_BitsStored, 16);
    dataset->putAndInsertUint16(DCM_HighBit, 15);
    dataset->putAndInsertUint16(DCM_PixelRepresentation, 0);
    dataset->putAndInsertString(DCM_Exposure, totalEventsStr.c_str());
    
    std::vector<Uint16> pixelData(rows * cols);
    for (int y = 0; y < numPixelsY; ++y) {
        for (int x = 0; x < numPixelsX; ++x) {
            double value = SensibilityOffset + SensibilityFactor * output[y][x];
            pixelData[y *cols+x]  = static_cast<Uint16>(std::max(0.0, std::min(65535.0, value)));
            
            }
        }

    dataset->putAndInsertUint16Array(DCM_PixelData, pixelData.data(), rows * cols);

    OFCondition status = dcmfileformat.saveFile((outputFolder+name+".dcm").c_str(), EXS_LittleEndianExplicit);
    if (status.bad()) {
        G4cerr << "Error: DCM file can't be saved" << G4endl;
        return;
    }
    
} 

//----------------------------------------------------------------------
int GmGetImage::EnergyToCharge(double energyDeposited){
    //double eVtoCharge = CLHEP::RandGauss::shoot(pairCreationEnergyEV, SwankFactor);
    double eVtoCharge = pairCreationEnergyEV;
    return std::floor(energyDeposited /eVtoCharge);
}

G4Material* GmGetImage::LoadCustomMaterial(const std::string& materialName, const std::string& filePath){

    std::filesystem::path full_path = std::filesystem::current_path() / filePath;
    std::ifstream file(full_path.string());
   // G4cout << full_path << G4endl;
    if (!file.is_open()){
        G4Exception("LoadCustomMaterial", "FileNotFound", FatalException, ("Could not open material file: " + full_path.string()).c_str());
    }
    std::string line;
    bool found = false;
    std::string name;
    double density = 0.0;
    int nComponents = 0;
    std::vector<std::pair<G4Material*, double>> components;
    G4NistManager* nist = G4NistManager::Instance();


    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string tag;
        iss >> tag;

        if (tag == ":MIXT"){
            iss >> std::quoted(name) >> density >> nComponents;
            //G4cout << name << density << nComponents << G4endl;
            if (name != materialName) {
                found = false;
                continue;
            }
            found = true;
            components.clear();
            for (int i = 0; i < nComponents && std::getline(file, line); ++i) {
                std::istringstream elemStream(line);
                std::string elemName;
                double fraction;
                elemStream >> elemName >> fraction;
                G4Material* elem = nist->FindOrBuildMaterial(elemName);
                //G4cout << "Element: " << elem->GetName() << G4endl;
                static G4EmCalculator emCal;

                /*/
                G4double lambda_phot = emCal.ComputeMeanFreePath(12*keV, G4Gamma::Gamma(), "phot", elem);
                G4double lambda_compt = emCal.ComputeMeanFreePath(12*keV, G4Gamma::Gamma(), "compt", elem);
                G4double lambda_conv = emCal.ComputeMeanFreePath(12*keV, G4Gamma::Gamma(), "conv", elem);
                G4double lambda_tot = 1/(1/lambda_phot + 1/lambda_compt + 1/lambda_conv);
                G4cout << "Element MFP Phot: " << lambda_phot/mm << G4endl;
                G4cout << "Element MFP Total: " << lambda_tot/mm << G4endl;
                /*/
                if (!elem) {
                    G4cerr << "ERROR: Could not find element '" << elemName << "'" << G4endl;
                }
                components.emplace_back(elem, fraction);
            }
            break;
        }
        
    }
    if (!found || components.empty()) {
        G4cerr << "Material " << materialName << "not found in file!" << G4endl;
        return nullptr;
    }
    
    G4Material* mat = new G4Material(materialName, density*g/cm3, components.size());
    //G4cout << materialName << density << components.size() << G4endl;
    for (const auto& comp : components) {

        //G4cout << "First: "<< comp.first << "Second: " << comp.second << G4endl;

        mat->AddMaterial(comp.first, comp.second);
    }

    return mat;
}

