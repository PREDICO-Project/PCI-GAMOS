#########################
## AAPM Breast Phantom ##
#########################

# Breast Skin semi-circular cylinder of thickness 50 mm and radius 100mm. Centered in the y-direction, chest wall at the edge of scoring plane. 595 mm below x-ray source

# Breast tissue semi-circular cylinder, concentric with the breast skin, 46mm thickness, 98mm radius, ie 2mm of skin envelope. Centered in breast skin


#Rot Matrix
:ROTM RM1 0. 0. 0.
:ROTM RM2 -180. 0. 0.
:ROTM RM3 0. 0. 90.


################
##   BREAST   ##
################

:P ZBreast 620.
:P XBreastCenter -70.
:P YBreastCenter 0.

:P tissueRadius 98.
:P skinRadius 100.

:P tissueThickness 46.
:P skinThickness 50.

:P skinDepth 2.0

###########################
##  Compression paddle   ##
###########################

:P compDimX 140.
:P compDimY 260.
:P compDimZ 2.

#############
##  Body   ##
#############

:P bodyDimX 170.
:P bodyDimY 300.
:P bodyDimZ 300.

##################
##  Dosimeters  ##
##################

:P dosDimX 20.
:P dosDimY 20.
:P dosDimZ 10.

########################
##   VOL definition   ##
########################

# Reference planes
#:VOLU refPlaneXY BOX 100. 150. 0.01 G4_Galactic
#:PLACE refPlaneXY 1 world RM1 0. 0. 620.
#:VIS refPlaneXY OFF

# Breast skin 
:VOLU skinExternal TUBS $tissueRadius $skinRadius 0.5*$skinThickness 0 180 "Breast_Skin"
:PLACE skinExternal 1 world RM3 $XBreastCenter $YBreastCenter $ZBreast
:COLOUR skinExternal 1.0 0.7688 0.2
:VIS skinExternal ON

:VOLU skinTop TUBS 0. $skinRadius 0.5*$skinDepth 0 180 "Breast_Skin"
:PLACE skinTop 1 world RM3 $XBreastCenter $YBreastCenter 595.0+0.5*$skinDepth
:COLOUR skinTop 1.0 0.7688 0.2
:VIS skinTop ON

:VOLU skinBottom TUBS 0. $skinRadius 0.5*$skinDepth 0 180 "Breast_Skin"
:PLACE skinBottom 1 world RM3 $XBreastCenter $YBreastCenter $ZBreast+0.5*$tissueThickness+0.5*$skinDepth
:COLOUR skinBottom 1.0 0.7688 0.2
:VIS skinBottom ON

# Breast Tissue
:VOLU tissue TUBS 0. $tissueRadius 0.5*$tissueThickness 270 180 "Breast_Tissue"
:PLACE tissue 1 world RM1 $XBreastCenter $YBreastCenter $ZBreast
:COLOUR tissue 1.0 0.34 0.2
:VIS tissue ON

# Compression paddles
:VOLU paddleTop BOX 0.5*$compDimX 0.5*$compDimY 0.5*$compDimZ "AAPM_PMMA" 
:PLACE paddleTop 1 world RM1 $XBreastCenter+0.5*$compDimX $YBreastCenter $ZBreast-0.5*$tissueThickness-$skinDepth-0.5*$compDimZ
:COLOUR paddleTop 0.6353 0.2 1.0
:VIS paddleTop OFF

:VOLU paddleBottom BOX 0.5*$compDimX 0.5*$compDimY 0.5*$compDimZ "AAPM_PMMA" 
:PLACE paddleBottom 1 world RM1 $XBreastCenter+0.5*$compDimX $YBreastCenter $ZBreast+0.5*$tissueThickness+$skinDepth+0.5*$compDimZ
:COLOUR paddleBottom 0.6353 0.2 1.0
:VIS paddleBottom OFF

# Body
:VOLU body BOX 0.5*$bodyDimX 0.5*$bodyDimY 0.5*$bodyDimZ G4_WATER
:PLACE body 1 world RM2 $XBreastCenter-0.5*$bodyDimX 0. $ZBreast
:COLOUR body 0.2 0.698 1.0
:VIS body OFF

# Scoring Plane
:VOLU scoringPlane BOX 0.5*$compDimX 0.5*$compDimY 0.5*$compDimZ "AAPM_AIR"
:PLACE scoringPlane 1 world RM1 $XBreastCenter+0.5*$compDimX $YBreastCenter $ZBreast+0.5*$tissueThickness+$skinDepth+15.0+0.5*$compDimZ
:COLOUR scoringPlane 1.0 0. 0.
:VIS scoringPlane ON

#################################
##   SCORING VOIs definition   ##
#################################

# Breast Scoring VOIS definition
:VOLU voi BOX 0.5*$dosDimX 0.5*$dosDimY 0.5*$dosDimZ "Breast_Tissue"
:PLACE voi 1 tissue RM1 50. 50. 0.
:PLACE voi 2 tissue RM1 20. 0. 0.
:PLACE voi 3 tissue RM1 50. 0. 0.
:PLACE voi 4 tissue RM1 80. 0. 0.
:PLACE voi 5 tissue RM1 50. -50. 0.
:PLACE voi 6 tissue RM1 50. 0. 15.
:PLACE voi 7 tissue RM1 50. 0. -15.
:COLOUR voi 1. 0. 0.
#:VIS voi OFF 

# Breast Scoring Dosimeter into the VOIS
:VOLU dosimeterVOI BOX 0.5*$dosDimX 0.5*$dosDimY 0.5*$dosDimZ "Breast_Tissue"
:PLACE dosimeterVOI 1 voi RM1 0. 0. 0.

###################################
##   SCORING ROIsFF definition   ##
###################################

#:VOLU roiFF BOX 0.5*$dosDimX 0.5*$dosDimY 1.0 G4_AIR
#:PLACE roiFF 1 scoringPlane RM1 -0.5*$compDimX+0.5*$dosDimX -0.5*$compDimY+0.5*$dosDimY 0.
#:PLACE roiFF 2 scoringPlane RM1 -0.5*$compDimX+0.5*$dosDimX -60. 0.
#:PLACE roiFF 3 scoringPlane RM1 -20. -60. 0.
#:PLACE roiFF 4 scoringPlane RM1 -0.5*$compDimX+0.5*$dosDimX 0. 0.
#:PLACE roiFF 5 scoringPlane RM1 -40. 0. 0.
#:PLACE roiFF 6 scoringPlane RM1  0. 0. 0.
#:PLACE roiFF 7 scoringPlane RM1 -0.5*$compDimX+0.5*$dosDimX 60. 0.
#:COLOUR roiFF 0. 0. 1.

#:VOLU dosimeterROIFF BOX 0.5*$dosDimX 0.5*$dosDimY 1.0 G4_AIR
#:PLACE dosimeterROIFF 1 roiFF RM1 0. 0. 0.

#:VOLU roiFF1 BOX 0.5*$dosDimX 0.5*$dosDimY 1.0 G4_AIR
#:PLACE roiFF1 1 scoringPlane RM1 -0.5*$compDimX+0.5*$dosDimX -0.5*$compDimY+0.5*$dosDimY 0.

#:VOLU roiFF2 BOX 0.5*$dosDimX 0.5*$dosDimY 1.0 G4_AIR
#:PLACE roiFF2 1 scoringPlane RM1 -0.5*$compDimX+0.5*$dosDimX -60. 0.

#:VOLU roiFF3 BOX 0.5*$dosDimX 0.5*$dosDimY 1.0 G4_AIR
#:PLACE roiFF3 1 scoringPlane RM1 -20. -60. 0.

#:VOLU roiFF4 BOX 0.5*$dosDimX 0.5*$dosDimY 1.0 G4_AIR
#:PLACE roiFF4 1 scoringPlane RM1 -0.5*$compDimX+0.5*$dosDimX 0. 0.

#:VOLU roiFF5 BOX 0.5*$dosDimX 0.5*$dosDimY 1.0 G4_AIR
#:PLACE roiFF5 1 scoringPlane RM1 -40. 0. 0.

#:VOLU roiFF6 BOX 0.5*$dosDimX 0.5*$dosDimY 1.0 G4_AIR
#:PLACE roiFF6 1 scoringPlane RM1  0. 0. 0.

#:VOLU roiFF7 BOX 0.5*$dosDimX 0.5*$dosDimY 1.0 G4_AIR
#:PLACE roiFF7 1 scoringPlane RM1 -0.5*$compDimX+0.5*$dosDimX 60. 0.

###################################
##   SCORING ROIsPS definition   ##
###################################

#:VOLU roiPS BOX 0.5*$dosDimX 0.5*$dosDimY 0.01 "Breast_Tissue"
#:PLACE roiPS 1 scoringPlane RM1 -40. -40. 0.
#:PLACE roiPS 2 scoringPlane RM1 -20. -20. 0.
#:PLACE roiPS 3 scoringPlane1 RM1 -40. 0. 0.
#:PLACE roiPS 4 scoringPlane RM1 -20. 0. 0.
#:PLACE roiPS 5 scoringPlane RM1  0. 0. 0.
#:PLACE roiPS 6 scoringPlane RM1  -20. 20. 0.
#:PLACE roiPS 7 scoringPlane RM1 -40. 40. 0.
#:COLOUR roiPS 0. 1. 0.

#:VOLU dosimeterROIPS BOX 0.5*$dosDimX 0.5*$dosDimY 0.01 "Breast_Tissue"
#:PLACE dosimeterROIPS 1 roiPS RM1 0. 0. 0.
