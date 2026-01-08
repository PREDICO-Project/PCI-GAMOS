######################## 
##   DETECTOR Setup   ##
######################## 
#Rot Matrix 
:ROTM RMD0 0. 0. 0.

 #Params definition
:P NPixelX 1490 
:P NPixelY 2031 
:P PixelSizeX 0.085 
:P PixelSizeY 0.085 
:P DetectorDepth 0.2 

:P DetectorMarkSize 5.0 
:P DetectorMarkDepth 0.2 
 
:P DetectorSizeX $NPixelX*$PixelSizeX 
:P DetectorSizeY $NPixelY*$PixelSizeY 

######################## 
##   VOL definition   ##
######################## 

#Detector structure Placed on Z-Axis 
:VOLU detector BOX 0.5*$DetectorSizeX 0.5*$DetectorSizeY 0.5*$DetectorDepth aSe 
:PLACE detector 1 world RMD0 0. 0. $SDD+0.5*$DetectorDepth 
:COLOUR detector 0. 0. 1. 
#:VIS detector OFF 

#Detector Mark 
:VOLU DetectorMark BOX 0.5*$DetectorMarkSize 0.5*$DetectorMarkSize 0.5*$DetectorMarkDepth G4_Pb 
:PLACE DetectorMark 1 detector RMD0 0.5*($DetectorSizeX-$DetectorMarkSize) 0.5*($DetectorSizeY-$DetectorMarkSize) 0. 
:COLOUR DetectorMark 0.2 0.4 0.1 

#Stopping photons VOLUME 
:VOLU stopVOLU BOX $DetectorSizeX $DetectorSizeY 0.5*$DetectorDepth G4_Pb 
:PLACE stopVOLU 1 world  RMD0 0. 0. $SDD+$DetectorDepth 
:COLOUR stopVOLU 0. 1. 0 
